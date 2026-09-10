#ftype: python

configfile: "config.yaml"

import os
from snakemake.exceptions import WorkflowError

CANONICAL_CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

# --- Config Validation ---
if not config.get("input_vcf"):
    raise WorkflowError("Configuration error: 'input_vcf' is not specified in config.yaml.")

# --- Intermediate Output Handling ---
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)

def intermediate(path):
    """Mark path as temp() unless keep_intermediate is set to True."""
    return path if KEEP_INTERMEDIATE else temp(path)

# --- Main Rules ---
rule all:
    input: 
        expand("results/{prefix}.sites.all.vcf.gz", prefix=config.get("output_prefix", "cohort")),
        expand("results/{prefix}.sites.pass.vcf.gz", prefix=config.get("output_prefix", "cohort")),
        expand("results/{prefix}.full_sample_data.vcf.gz", prefix=config.get("output_prefix", "cohort"))

# --- Resources & Common Steps ---
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
    output: vcf=intermediate("results/temp/01_selected.vcf.gz")
    threads: 4
    shell:
        """
        # If the map file is not empty, rename chromosomes first
        if [ -s {input.chr_map} ]; then
            bcftools annotate --threads {threads} --rename-chrs {input.chr_map} {input.vcf} -Ou | \
            bcftools view --threads {threads} -t {CANONICAL_CHRS} -Ou | \
            bcftools view --threads {threads} -e 'ALT="*"' -O z -o {output.vcf}
        else
            bcftools view --threads {threads} -t {CANONICAL_CHRS} {input.vcf} -Ou | \
            bcftools view --threads {threads} -e 'ALT="*"' -O z -o {output.vcf}
        fi
        """

rule normalize_and_split:
    input: vcf="results/temp/01_selected.vcf.gz", ref="resources/hg38.fa"
    output: vcf=intermediate("results/temp/02_normalized.vcf.gz")
    threads: 4
    shell:
        """
        bcftools norm --threads {threads} --force -m -any -f {input.ref} -Ou {input.vcf} | \
        bcftools annotate --threads {threads} -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou | \
        bcftools norm --threads {threads} --rm-dup exact -Oz -o {output.vcf}
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
        vcf=intermediate("results/temp/03_raw_stats.vcf.gz"),
        tbi=intermediate("results/temp/03_raw_stats.vcf.gz.tbi")
    threads: 4
    shell:
        """
        bcftools +fill-tags --threads {threads} {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view --threads {threads} -O z -o {output.vcf}
        bcftools index --threads {threads} -t {output.vcf}
        """

rule genotype_masking:
    input: "results/temp/03_raw_stats.vcf.gz"
    output: 
        vcf=intermediate("results/temp/04_masked.vcf.gz"),
        tbi=intermediate("results/temp/04_masked.vcf.gz.tbi")
    threads: 4
    params:
        min_gq=config["qc_thresholds"]["min_gq"],
        min_dp=config["qc_thresholds"]["min_dp"],
        min_ab=config["qc_thresholds"]["ab_ratio"],
        max_ab = 1 - config["qc_thresholds"]["ab_ratio"]
    shell:
        """
        bcftools +setGT --threads {threads} {input} -O z -o {output.vcf} -- -t q -n . -i 'FMT/GQ < {params.min_gq} | FMT/DP < {params.min_dp} | (GT="het" & (FMT/AD[*:0] + FMT/AD[*:1]) > 0 & ((FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) < {params.min_ab} | (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) > {params.max_ab}))'
        bcftools index --threads {threads} -t {output.vcf}
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
        # Node: Using bcftools call --ploidy is prone to human error.
        #       It outputs on stderr and exit code ($?) is 255
        #       Wrong invocation generates invalid files and exit code is again 255
        #       bcftools +fixploidy accepts empty and even corrupt files making these errors silent.
        # TODO: It's better to pregenerate these files or maybe hard code them in the repository it self.
        """
        bcftools call --ploidy GRCh38? 2> {output.ploidy} \
        || true # due to snakemake strict mode we have to return true (hide 255 exit code)
        """

rule fix_ploidy:
    input: 
        vcf="results/temp/04_masked.vcf.gz",
        samples=config["sample_info"],
        ploidy=rules.generate_ploidy_rules.output.ploidy
    output: 
        vcf=intermediate("results/temp/05_ploidy_fixed.vcf.gz"),
        tbi=intermediate("results/temp/05_ploidy_fixed.vcf.gz.tbi"),
        sex_map="results/temp/sex_map.txt"
    threads: 4
    shell:
        """
        # Cleanly generate the standardized sex map
        python scripts/generate_sex_map.py {input.samples} {output.sex_map}
        
        bcftools +fixploidy --threads {threads} {input.vcf} -Oz -o {output.vcf} -- -s {output.sex_map} -p {input.ploidy}
        bcftools index --threads {threads} -t {output.vcf}
        """

rule annotate_final_vcf:
    input: 
        vcf="results/temp/05_ploidy_fixed.vcf.gz", 
        groups="results/temp/groups_final.txt"
    output: 
        vcf=intermediate("results/temp/06_final_stats.vcf.gz"),
        tbi=intermediate("results/temp/06_final_stats.vcf.gz.tbi")
    threads: 4
    shell:
        """
        bcftools +fill-tags --threads {threads} {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view --threads {threads} -O z -o {output.vcf}
        bcftools index --threads {threads} -t {output.vcf}
        """

rule variant_qc_tagging:
    input: "results/temp/06_final_stats.vcf.gz"
    output: 
        vcf="results/{prefix}.full_sample_data.vcf.gz",
        tbi="results/{prefix}.full_sample_data.vcf.gz.tbi"
    threads: 4
    params:
        qual=config["qc_thresholds"]["qual"],
        qd=config["qc_thresholds"]["qd"],
        mq=config["qc_thresholds"]["mq"],
        fs=config["qc_thresholds"]["fs"],
        readpos=config["qc_thresholds"]["readpos"],
        hwe=config["qc_thresholds"]["hwe"],
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
            --max_missing {params.max_missing} \
            --min_dp {params.min_dp} \
            --min_gq {params.min_gq} \
            --ab_ratio {params.ab_ratio}
        bcftools index --threads {threads} -t {output.vcf}
        """

rule create_sites_vcf:
    input: "results/{prefix}.full_sample_data.vcf.gz"
    output:
        all_sites="results/{prefix}.sites.all.vcf.gz",
        all_sites_tbi="results/{prefix}.sites.all.vcf.gz.tbi",
        pass_sites="results/{prefix}.sites.pass.vcf.gz",
        pass_sites_tbi="results/{prefix}.sites.pass.vcf.gz.tbi"
    threads: 4
    shell:
        """
        bcftools view --threads {threads} -G -O z -o {output.all_sites} {input}
        bcftools index --threads {threads} -t {output.all_sites}
        
        # Filter for PASS status.
        bcftools view --threads {threads} -G -i 'INFO/QC_STATUS="PASS"' {input} -O z -o {output.pass_sites}
        bcftools index --threads {threads} -t {output.pass_sites}
        """
