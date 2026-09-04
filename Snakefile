#ftype: python

configfile: "config.yaml"

import os
from snakemake.exceptions import WorkflowError

CANONICAL_CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

# --- Config Validation ---
if not config.get("input_vcf"):
    raise WorkflowError("Configuration error: 'input_vcf' is not specified in config.yaml.")

# --- Main Rules ---
rule all:
    input: 
        expand("results/{prefix}.sites.all.vcf.gz", prefix=config.get("output_prefix", "cohort")),
        expand("results/{prefix}.sites.pass.vcf.gz", prefix=config.get("output_prefix", "cohort")),
        expand("results/intermediate/{prefix}.full_sample_data.vcf.gz", prefix=config.get("output_prefix", "cohort"))

# --- Resources & Common Steps ---
rule index_vcf:
    input: "{filename}.vcf.gz"
    output: "{filename}.vcf.gz.tbi"
    shell:
        "tabix -pvcf {input}"

rule prepare_reference:
    input: src=lambda w: config.get("local_resources", {}).get("ref_fasta") or "resources/downloaded_hg38.fa.gz"
    output: fasta="resources/hg38.fa", fai="resources/hg38.fa.fai"
    shell:
        """
        if [[ "{input.src}" == *.gz ]]; then 
            gunzip -c {input.src} > {output.fasta}
        else 
            ln -sf $(realpath "{input.src}") {output.fasta}
        fi
        samtools faidx {output.fasta}
        """

rule download_reference_source:
    output: "resources/downloaded_hg38.fa.gz"
    params: url=config["resources"]["ref_fasta_url"]
    shell: "wget --tries=3 -O {output} {params.url}"

rule generate_chr_map:
    input:
        vcf=config["input_vcf"]
    output:
        map_file="results/temp/chr_map.txt"
    run:
        from cyvcf2 import VCF
        vcf = VCF(input.vcf)
        has_chr = any(seq.startswith("chr") for seq in vcf.seqnames)
        vcf.close()
        
        with open(output.map_file, "w") as f:
            if not has_chr:
                # Write mapping from Ensembl to UCSC (e.g. "1 chr1")
                for i in range(1, 23):
                    f.write(f"{i}\tchr{i}\n")
                f.write("X\tchrX\n")
                f.write("Y\tchrY\n")
                f.write("MT\tchrM\n")
            else:
                # No-op (empty map file means no renaming)
                pass

# --- Pipeline Rules ---
rule select_regions:
    input:
        vcf=config["input_vcf"],
        chr_map="results/temp/chr_map.txt"
    output: vcf="results/temp/01_selected.vcf.gz"
    shell:
        """
        # If the map file is not empty, rename chromosomes first
        if [ -s {input.chr_map} ]; then
            bcftools annotate --rename-chrs {input.chr_map} {input.vcf} -Ou | \
            bcftools view -t {CANONICAL_CHRS} -Ou | \
            bcftools view -e 'ALT="*"' -O z -o {output.vcf}
        else
            bcftools view -t {CANONICAL_CHRS} {input.vcf} | \
            bcftools view -e 'ALT="*"' -O z -o {output.vcf}
        fi
        """

rule normalize_and_split:
    input: vcf="results/temp/01_selected.vcf.gz", ref="resources/hg38.fa"
    output: vcf="results/temp/02_normalized.vcf.gz"
    shell:
        """
        bcftools norm --force -m -any -f {input.ref} -Ou {input.vcf} | \
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou | \
        bcftools norm --rm-dup exact -Oz -o {output.vcf}
        """


# --- QC Workflow ---

rule generate_groups_raw:
    input: samples=config["sample_info"]
    output: 
        groups="results/temp/groups_raw.txt"
    params: suffix="raw"
    shell:
        "python scripts/generate_groups.py {input.samples} {output.groups} --suffix {params.suffix}"

rule annotate_raw_vcf:
    input: 
        vcf="results/temp/02_normalized.vcf.gz", 
        groups="results/temp/groups_raw.txt"
    output: 
        vcf="results/temp/03_raw_stats.vcf.gz",
        tbi="results/temp/03_raw_stats.vcf.gz.tbi"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule genotype_masking:
    input: "results/temp/03_raw_stats.vcf.gz"
    output: 
        vcf="results/temp/04_masked.vcf.gz",
        tbi="results/temp/04_masked.vcf.gz.tbi"
    params:
        min_gq=config["qc_thresholds"]["min_gq"],
        min_dp=config["qc_thresholds"]["min_dp"],
        min_ab=config["qc_thresholds"]["ab_ratio"],
        max_ab = 1 - config["qc_thresholds"]["ab_ratio"]
    shell:
        """
        bcftools +setGT {input} -O z -o {output.vcf} -- -t q -n . -i 'FMT/GQ < {params.min_gq} | FMT/DP < {params.min_dp} | (GT="het" & (FMT/AD[*:0] + FMT/AD[*:1]) > 0 & ((FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) < {params.min_ab} | (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) > {params.max_ab}))'
        tabix -p vcf {output.vcf}
        """

rule generate_groups_final:
    input: samples=config["sample_info"]
    output: 
        groups="results/temp/groups_final.txt"
    params: suffix=""
    shell:
        "python scripts/generate_groups.py {input.samples} {output.groups}"

rule generate_ploidy_rules:
    output:
        ploidy="results/temp/ploidy_rules.txt"
    shell:
        # BCFtools returns non-zero when querying built-in ploidy definitions;
        # redirect stdout cleanly and use || true to prevent shell abort.
        "bcftools call --ploidy GRCh38? 1> {output.ploidy} 2>/dev/null || true"

rule fix_ploidy:
    input: 
        vcf="results/temp/04_masked.vcf.gz",
        samples=config["sample_info"],
        ploidy=rules.generate_ploidy_rules.output.ploidy
    output: 
        vcf="results/temp/05_ploidy_fixed.vcf.gz",
        tbi="results/temp/05_ploidy_fixed.vcf.gz.tbi",
        sex_map="results/temp/sex_map.txt"
    shell:
        """
        # Cleanly generate the standardized sex map
        python scripts/generate_sex_map.py {input.samples} {output.sex_map}
        
        bcftools +fixploidy {input.vcf} -Oz -o {output.vcf} -- -s {output.sex_map} -p {input.ploidy}
        tabix -p vcf {output.vcf}
        """

rule annotate_final_vcf:
    input: 
        vcf="results/temp/05_ploidy_fixed.vcf.gz", 
        groups="results/temp/groups_final.txt"
    output: 
        vcf="results/temp/06_final_stats.vcf.gz",
        tbi="results/temp/06_final_stats.vcf.gz.tbi"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule variant_qc_tagging:
    input: "results/temp/06_final_stats.vcf.gz"
    output: 
        vcf="results/intermediate/{prefix}.full_sample_data.vcf.gz",
        tbi="results/intermediate/{prefix}.full_sample_data.vcf.gz.tbi"
    params:
        qual=config["qc_thresholds"]["qual"],
        qd=config["qc_thresholds"]["qd"],
        mq=config["qc_thresholds"]["mq"],
        fs=config["qc_thresholds"]["fs"],
        readpos=config["qc_thresholds"]["readpos"],
        hwe=config["qc_thresholds"]["hwe"],
        maf=config["qc_thresholds"]["maf"],
        max_missing=config["qc_thresholds"]["max_missing"],
        min_dp=config["qc_thresholds"]["min_dp"],
        min_gq=config["qc_thresholds"]["min_gq"],
        ab_ratio=config["qc_thresholds"]["ab_ratio"]
    shell:
        """
        python scripts/tag_variant_qc.py {input} {output.vcf} \
            --qual {params.qual} \
            --qd {params.qd} \
            --mq {params.mq} \
            --fs {params.fs} \
            --readpos {params.readpos} \
            --hwe {params.hwe} \
            --maf {params.maf} \
            --max_missing {params.max_missing} \
            --min_dp {params.min_dp} \
            --min_gq {params.min_gq} \
            --ab_ratio {params.ab_ratio}
        tabix -p vcf {output.vcf}
        """

rule create_sites_vcf:
    input: "results/intermediate/{prefix}.full_sample_data.vcf.gz"
    output:
        all_sites="results/{prefix}.sites.all.vcf.gz",
        all_sites_tbi="results/{prefix}.sites.all.vcf.gz.tbi",
        pass_sites="results/{prefix}.sites.pass.vcf.gz",
        pass_sites_tbi="results/{prefix}.sites.pass.vcf.gz.tbi"
    shell:
        """
        bcftools view -G -O z -o {output.all_sites} {input}
        tabix -p vcf {output.all_sites}
        
        # Filter for PASS status.
        bcftools view -G -i 'INFO/QC_STATUS="PASS"' {input} -O z -o {output.pass_sites}
        tabix -p vcf {output.pass_sites}
        """
