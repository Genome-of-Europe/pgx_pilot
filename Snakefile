#ftype: python

configfile: "config.yaml"

import os, sys

# --- Config Validation ---
if not config.get("input_vcf") or not os.path.exists(config["input_vcf"]):
    sys.exit(f"ERROR: input_vcf file not found at: '{config.get('input_vcf')}'.")

if not config.get("regions_bed"):
    sys.exit("ERROR: 'regions_bed' parameter is required in config.yaml.")

# --- Helper Functions ---
def get_selection_inputs(wildcards):
    files = {"vcf": config["input_vcf"]}
    bed = config.get("regions_bed")
    if bed: 
        files["bed"] = "resources/selection_targets.bed" if bed.startswith(("http", "ftp")) else bed
    return files

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
    output: fasta="resources/hg38.fa", dct="resources/hg38.dict"
    shell:
        """
        if [[ "{input.src}" == *.gz ]]; then gunzip -c {input.src} > {output.fasta}; else ln -sf {input.src} {output.fasta}; fi
        samtools faidx {output.fasta}
        samtools dict {output.fasta} > {output.dct}
        """

rule download_reference_source:
    output: "resources/downloaded_hg38.fa.gz"
    params: url=config["resources"]["ref_fasta_url"]
    shell: "wget --tries=3 -O {output} {params.url}"

# --- Pipeline Rules ---
rule select_regions:
    input: unpack(get_selection_inputs)
    output: "results/temp/01_selected.vcf.gz"
    params: 
        bed=config.get("regions_bed", ""), 
        bed_path=lambda w: "resources/selection_targets.bed" if config.get("regions_bed", "").startswith(("http", "ftp")) else config.get("regions_bed", "")
    shell:
        """
        bcftools view -R {params.bed_path} {input.vcf} | bcftools view -e 'ALT="*"' -O z -o {output}
        tabix -p vcf {output}
        """

rule normalize_and_split:
    input: vcf="results/temp/01_selected.vcf.gz", ref="resources/hg38.fa"
    output: "results/temp/02_normalized.vcf.gz"
    shell: "bcftools norm -m -any -f {input.ref} -O z -o {output} {input.vcf} && tabix -p vcf {output}"

rule download_selection_bed:
    output: "resources/selection_targets.bed"
    params: url=config.get("regions_bed")
    shell: "wget --tries=3 -O {output} {params.url}"

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
    output: vcf="results/temp/03_raw_stats.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule genotype_masking:
    input: "results/temp/03_raw_stats.vcf.gz"
    output: "results/temp/04_masked.vcf.gz"
    params:
        min_gq=config["qc_thresholds"]["min_gq"],
        min_dp=config["qc_thresholds"]["min_dp"],
        min_ab=config["qc_thresholds"]["ab_ratio"],
        max_ab = 1 - config["qc_thresholds"]["ab_ratio"]
    shell:
        """
        bcftools +setGT {input} -O z -o {output} -- -t q -n . -i 'FMT/GQ < {params.min_gq} | FMT/DP < {params.min_dp} | (GT="het" & (FMT/AD[*:0] + FMT/AD[*:1]) > 0 & ((FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) < {params.min_ab} | (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) > {params.max_ab}))'
        tabix -p vcf {output}
        """

rule generate_groups_final:
    input: samples=config["sample_info"]
    output: 
        groups="results/temp/groups_final.txt"
    params: suffix=""
    shell:
        "python scripts/generate_groups.py {input.samples} {output.groups}"

rule annotate_final_vcf:
    input: 
        vcf="results/temp/04_masked.vcf.gz", 
        groups="results/temp/groups_final.txt"
    output: vcf="results/temp/05_final_stats.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet | \
        bcftools view -O z -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

rule variant_qc_tagging:
    input: "results/temp/05_final_stats.vcf.gz"
    output: "results/intermediate/{prefix}.full_sample_data.vcf.gz"
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
        python scripts/tag_variant_qc.py {input} {output} \
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
        tabix -p vcf {output}
        """

rule create_sites_vcf:
    input: "results/intermediate/{prefix}.full_sample_data.vcf.gz"
    output:
        all_sites="results/{prefix}.sites.all.vcf.gz",
        pass_sites="results/{prefix}.sites.pass.vcf.gz"
    shell:
        """
        bcftools view -G -O z -o {output.all_sites} {input}
        tabix -p vcf {output.all_sites}
        
        # Filter for PASS status.
        bcftools view -G -i 'INFO/QC_STATUS="PASS"' {input} -O z -o {output.pass_sites}
        tabix -p vcf {output.pass_sites}
        """
