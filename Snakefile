#ftype: python
import os
from textwrap import dedent
from snakemake.exceptions import WorkflowError
from utils.vcf import vcf_format_from_output_type

configfile: "config.yaml"

CANONICAL_CHRS = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

# --- Config Validation ---
INPUT_VCF = config.get("input_vcf")
if not INPUT_VCF:
    raise WorkflowError("Configuration error: 'input_vcf' is not specified in config.yaml.")

QC = config.get("qc_thresholds", {})
QC_TAG_PARAMS = {
    "qual": QC["qual"],
    "qd": QC["qd"],
    "mq": QC["mq"],
    "fs": QC["fs"],
    "readpos": QC["readpos"],
    "hwe": QC["hwe"],
    "max_missing": QC["max_missing"],
    "min_dp": QC["min_dp"],
    "min_gq": QC["min_gq"],
    "ab_ratio": QC["ab_ratio"],
}

REGIONS = config.get("regions", "")
F_REGIONS = "" if not REGIONS else f"-r {REGIONS}"

# --- Compression/Decompression Threads ---
READ_THREADS = config.get("read_thread", 1)
WRITE_THREADS = config.get("write_thread", 2)
F_READ_THREADS = f"--threads {READ_THREADS}"
F_WRITE_THREADS = f"--threads {WRITE_THREADS}"

# --- Config Output Types ---
TMP_O = config.get("tmp_vcf_O", "b1")
TMP_EXT = vcf_format_from_output_type(TMP_O)
OUT_O = config.get("out_vcf_O", "z4")
OUT_EXT = vcf_format_from_output_type(OUT_O)

# --- Intermediate Output Handling ---
KEEP_INTERMEDIATE = config.get("keep_intermediate", False)
INDEX_TMP_FILES = config.get("index_tmp_files", False)
F_WRITE_TMP_INDEX = "--write-index" if KEEP_INTERMEDIATE and INDEX_TMP_FILES else ""

# --- Global Paths & Helpers ---
RESULTS_DIR = config.get("results_dir", "results")
TEMP_DIR = f"{RESULTS_DIR}/temp"
LOGS_DIR = f"{RESULTS_DIR}/logs"
SCRIPTS_DIR = f"{workflow.basedir}/scripts"

def resource(name: str) -> str:
    return f"resources/{name}"

def tmp(name: str) -> str:
    return f"{TEMP_DIR}/{name}"

def tmp_vcf(name: str) -> str:
    return tmp(f"{name}.{TMP_EXT}")

def intermediate(path: str, keep: bool = KEEP_INTERMEDIATE) -> str:
    """Mark path as temp() unless keep_intermediate is set to True."""
    return path if keep else temp(path)

def log_path(name: str) -> str:
    """Generate log file path for a given rule name."""
    return f"{LOGS_DIR}/{name}.log"

def add_log(cmd: str) -> str:
    """Wrap a shell command to stream combined stdout/stderr to console and log file."""
    return f"(\n{dedent(cmd).strip()}\n) 2>&1 | tee {{log}}"

# --- Global Bash Helper for Temporary BCF Output ---
shell.prefix(
    """
    write_tmp_vcf(){{
        bcftools view -o "$1" -O {TMP_O} {F_WRITE_THREADS} {F_WRITE_TMP_INDEX}
    }}
    """
)

# --- Main Rules ---
PREFIX = config.get("output_prefix", "cohort")

all_sites_path = expand(f"{RESULTS_DIR}/{{prefix}}.sites.all.{OUT_EXT}", prefix=PREFIX)
pass_sites_path = expand(f"{RESULTS_DIR}/{{prefix}}.sites.pass.{OUT_EXT}", prefix=PREFIX)
full_sample_data_path = expand(f"{RESULTS_DIR}/{{prefix}}.full_sample_data.{OUT_EXT}", prefix=PREFIX)

rule all:
    input:
        all_sites_path,
        pass_sites_path,
        full_sample_data_path

# --- Resources & Common Steps ---
rule prepare_reference:
    input: src=lambda w: config.get("local_resources", {}).get("ref_fasta") or resource("downloaded_hg38.fa.gz")
    output:
        fasta=resource("hg38.fa"),
        fai=resource("hg38.fa.fai")
    shell:
        """
        rm -f {output.fasta}
        if [[ "{input.src}" == *.gz ]]; then 
            gunzip -c {input.src} > {output.fasta}
        else 
            ln -sf $(realpath "{input.src}") {output.fasta}
        fi
        samtools faidx {output.fasta}
        """

rule download_reference_source:
    output: resource("downloaded_hg38.fa.gz")
    params: url=config["resources"]["ref_fasta_url"]
    shell: "wget --tries=3 -O {output} {params.url}"

rule generate_chr_rename:
    input: INPUT_VCF
    output: tmp("chr_rename_map.txt")
    shell:
        "python {SCRIPTS_DIR}/chr_rename_map.py {input} {output}"

# --- Pipeline Rules ---
rule select_regions:
    input:
        vcf=INPUT_VCF,
        chr_rename_map=tmp("chr_rename_map.txt")
    output: vcf=intermediate(tmp_vcf("01_selected"))
    log: log_path("01_selected")
    shell: add_log("""
        bcftools view {F_READ_THREADS} {F_REGIONS} {input.vcf} -Ou |
        {{
            # If the map file is not empty, rename chromosomes
            if [ -s {input.chr_rename_map} ]; then
                bcftools annotate --rename-chrs {input.chr_rename_map} -Ou
            else
                # If empty, pass the stream through unmodified
                cat
            fi
        }} |
        bcftools view -t {CANONICAL_CHRS} -Ou |
        write_tmp_vcf {output.vcf}
    """)

rule normalize_and_split:
    input:
        vcf=tmp_vcf("01_selected"),
        ref=resource("hg38.fa")
    output: vcf=intermediate(tmp_vcf("02_normalized"))
    log: log_path("02_normalized")
    shell: add_log("""
        bcftools norm {F_READ_THREADS} --force -m -any --check-ref eC -f {input.ref} -Ou {input.vcf} |
        bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -Ou |
        bcftools view {F_READ_THREADS} -e 'ALT="*"' -Ou |
        bcftools norm --rm-dup exact -Ou |
        write_tmp_vcf {output.vcf}
    """)

# --- QC Workflow ---

rule generate_groups_raw:
    input: samples=config["sample_info"]
    output:
        groups=tmp("groups_raw.txt")
    params: suffix="raw"
    shell:
        "python {SCRIPTS_DIR}/generate_groups.py {input.samples} {output.groups} --suffix {params.suffix}"

rule annotate_raw_vcf:
    input:
        vcf=tmp_vcf("02_normalized"),
        groups=tmp("groups_raw.txt")
    output:
        vcf=intermediate(tmp_vcf("03_raw_stats"))
    log: log_path("03_raw_stats")
    shell: add_log("""
        bcftools view {F_READ_THREADS} {input.vcf} -Ou |
        bcftools +fill-tags -Ou -- -S {input.groups} -t AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF |
        bcftools +fill-tags -Ou -- -S {input.groups} -t HWE,NS,F_MISSING,ExcHet |
        write_tmp_vcf {output.vcf}
    """)

rule genotype_masking:
    input: tmp_vcf("03_raw_stats")
    output:
        vcf=intermediate(tmp_vcf("04_masked"))
    log: log_path("04_masked")
    params:
        min_gq=QC["min_gq"],
        min_dp=QC["min_dp"],
        min_ab=QC["ab_ratio"],
        max_ab=1 - QC["ab_ratio"]
    shell: add_log("""
        bcftools +setGT {F_READ_THREADS} {input} -Ou -- -t q -n . -i 'FMT/GQ < {params.min_gq} | FMT/DP < {params.min_dp} | (GT="het" & (FMT/AD[*:0] + FMT/AD[*:1]) > 0 & ((FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) < {params.min_ab} | (FMT/AD[*:1])/(FMT/AD[*:0]+FMT/AD[*:1]) > {params.max_ab}))' |
        write_tmp_vcf {output.vcf}
    """)

rule generate_groups_final:
    input: samples=config["sample_info"]
    output:
        groups=tmp("groups_final.txt")
    params: suffix=""
    shell:
        "python {SCRIPTS_DIR}/generate_groups.py {input.samples} {output.groups}"

rule generate_sex_map:
    input: samples=config["sample_info"]
    output: sex_map=tmp("sex_map.txt")
    shell:
        "python {SCRIPTS_DIR}/generate_sex_map.py {input.samples} {output.sex_map}"

rule fix_ploidy:
    input:
        vcf=tmp_vcf("04_masked"),
        samples=config["sample_info"],
        ploidy=resource("ploidy_rules.txt"),
        sex_map=tmp("sex_map.txt")
    output:
        vcf=intermediate(tmp_vcf("05_ploidy_fixed"))
    log: log_path("05_ploidy_fixed")
    shell: add_log("""
        bcftools +fixploidy {F_READ_THREADS} {input.vcf} -Ou -- -s {input.sex_map} -p {input.ploidy} |
        write_tmp_vcf {output.vcf}
    """)

use rule annotate_raw_vcf as annotate_final_vcf with:
    input:
        vcf=tmp_vcf("05_ploidy_fixed"),
        groups=tmp("groups_final.txt")
    output:
        vcf=intermediate(tmp_vcf("06_final_stats"))
    log: log_path("06_final_stats")

rule variant_qc_tagging:
    input: tmp_vcf("06_final_stats")
    output:
        vcf=full_sample_data_path
    log: expand(log_path("{prefix}_07_variant_qc_tagging"), prefix=PREFIX)
    params:
        **QC_TAG_PARAMS
    shell: add_log("""
        python {SCRIPTS_DIR}/tag_variant_qc.py --in-threads {READ_THREADS} {input} -Ou \
            --qual {params.qual} \
            --qd {params.qd} \
            --mq {params.mq} \
            --fs {params.fs} \
            --readpos {params.readpos} \
            --hwe {params.hwe} \
            --max_missing {params.max_missing} \
            --min_dp {params.min_dp} \
            --min_gq {params.min_gq} \
            --ab_ratio {params.ab_ratio} |
        bcftools view {F_WRITE_THREADS} -O{OUT_O} -o {output.vcf} --write-index
    """)

rule create_sites_vcf:
    input:
        in_vcf=full_sample_data_path
    output:
        all_sites=all_sites_path,
        pass_sites=pass_sites_path
    log: expand(log_path("{prefix}_08_sites"), prefix=PREFIX)
    shell: add_log("""
        bcftools view {F_READ_THREADS} {input.in_vcf} -G -Ou |
        tee >(bcftools view {F_WRITE_THREADS} -O{OUT_O} -o {output.all_sites} --write-index) |
        bcftools view {F_WRITE_THREADS} -O{OUT_O} -o {output.pass_sites} -i 'INFO/QC_STATUS="PASS"' --write-index
    """)

onsuccess:
    from pathlib import Path
    if os.path.exists(LOGS_DIR):
        for f in Path(LOGS_DIR).glob("*.log"):
            if f.is_file() and f.stat().st_size == 0:
                f.unlink()

onerror:
    print(f"\\nWorkflow encountered an error. Check detailed step logs in: {LOGS_DIR}/\\n")
