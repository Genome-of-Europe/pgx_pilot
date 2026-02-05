#ftype: python

configfile: "config.yaml"

import os, sys, json, csv

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

def get_stats_inputs(wildcards):
    files = {"vcf": "results/temp/03_qc.vcf.gz", "samples": config["sample_info"]}
    return files

# --- Main Rules ---
rule all:
    input: 
        expand("results/{country}.sites.all.vcf.gz", country=config.get("country_code", "PT")),
        expand("results/{country}.sites.pass.vcf.gz", country=config.get("country_code", "PT")),
        expand("results/intermediate/{country}.full_sample_data.vcf.gz", country=config.get("country_code", "PT"))

rule finalize_output:
    input: 
        all_sites="results/{country}.sites.all.vcf.gz",
        pass_sites="results/{country}.sites.pass.vcf.gz"
    output:
        "results/final_annotated.sites.vcf.gz" # Optional legacy output
    shell: "cp {input.pass_sites} {output} && tabix -p vcf {output}"

# --- Resources & Common Steps ---
rule prepare_reference:
    input: src=lambda w: config.get("local_resources", {}).get("ref_fasta") or "resources/downloaded_GRCh38.fa.gz"
    output: fasta="resources/GRCh38.fa", dct="resources/GRCh38.dict"
    shell:
        """
        if [[ "{input.src}" == *.gz ]]; then gunzip -c {input.src} > {output.fasta}; else ln -sf {input.src} {output.fasta}; fi
        samtools faidx {output.fasta}
        samtools dict {output.fasta} > {output.dct}
        """

rule download_reference_source:
    output: "resources/downloaded_GRCh38.fa.gz"
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
    input: vcf="results/temp/01_selected.vcf.gz", ref="resources/GRCh38.fa"
    output: "results/temp/02_normalized.vcf.gz"
    shell: "bcftools norm -m -any -f {input.ref} -O z -o {output} {input.vcf} && tabix -p vcf {output}"

rule download_selection_bed:
    output: "resources/selection_targets.bed"
    params: url=config.get("regions_bed")
    shell: "wget --tries=3 -O {output} {params.url}"

# --- QC Workflow ---

rule annotate_vcf:
    input: vcf="results/temp/02_normalized.vcf.gz", samples=config["sample_info"]
    output: vcf="results/temp/03_raw_stats.vcf.gz"
    params: country=config.get("country_code", "PT")
    shell:
        "python scripts/annotate_vcf.py {input.vcf} {output.vcf} {input.samples} {params.country} --raw && tabix -p vcf {output.vcf}"

rule genotype_masking:
    input: "results/temp/03_raw_stats.vcf.gz"
    output: "results/temp/04_masked.vcf.gz"
    params:
        min_gq=config["qc_thresholds"]["min_gq"],
        min_dp=config["qc_thresholds"]["min_dp"],
        min_ab=config["qc_thresholds"]["ab_ratio"]
    shell:
        """
        bcftools +setGT {input} -O z -o {output} -- -t q -n . -e 'FMT/GQ < {params.min_gq} || FMT/DP < {params.min_dp} || ((FMT/AD[0:0] + FMT/AD[0:1]) > 0 && (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1])) < {params.min_ab})'
        tabix -p vcf {output}
        """

rule calc_final_stats:
    input: vcf="results/temp/04_masked.vcf.gz", samples=config["sample_info"]
    output: vcf="results/temp/05_final_stats.vcf.gz"
    params: country=config.get("country_code", "PT")
    shell:
        "python scripts/annotate_vcf.py {input.vcf} {output.vcf} {input.samples} {params.country} && tabix -p vcf {output.vcf}"

rule variant_qc_tagging:
    input: "results/temp/05_final_stats.vcf.gz"
    output: "results/intermediate/{country}.full_sample_data.vcf.gz"
    params:
        args="" 
    shell:
        "python scripts/tag_variant_qc.py {input} {output} {params.args} && tabix -p vcf {output}"

rule create_sites_vcfs:
    input: "results/intermediate/{country}.full_sample_data.vcf.gz"
    output:
        all_sites="results/{country}.sites.all.vcf.gz",
        pass_sites="results/{country}.sites.pass.vcf.gz"
    shell:
        """
        bcftools view -G -O z -o {output.all_sites} {input}
        tabix -p vcf {output.all_sites}
        
        # Filter for PASS status.
        bcftools view -G -i 'INFO/QC_STATUS="PASS"' {input} -O z -o {output.pass_sites}
        tabix -p vcf {output.pass_sites}
        """
