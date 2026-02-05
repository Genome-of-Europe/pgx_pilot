# GDI Variant Harmonization Pipeline

**Project:** WP6 GoE  
**Author:** Hugo Martiniano

## 1. Overview
This pipeline aims to standardize genomic variant processing across multiple sites. It automates:

1. Filtering variants based on specified genomic regions.
2. Splitting multiallelic sites and left-aligning indels.
3. Genotype masking and variant filtering.
4. Calculating Allele Frequencies (AF) stratified by Sex and Filter Status.
5. Generating a sites-only VCF output.

The goal is to allow partners to run the exact same workflow on their local data and produce standardized, comparable outputs.

## 2. Quick Start
The pipeline can be run as Docker container, or directly from the code in the GitHub repository.

**Prerequisites:** [Docker](https://docs.docker.com/get-docker/) or [Apptainer/Singularity](https://apptainer.org/).

### Pull the Image (NOT Working yet - use the Dockerfile or run snakemake locally)
```bash
docker pull goe/gdi-pipeline:latest
```

### Build Instructions
#### Using conda
1. Install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install).
2. Create and activate conda env
```bash
conda create -n pgx_pilot -f env.yml
conda activate pgx_pilot
```

#### Using docker
```bash
git clone https://github.com/your-repo/gdi_pilot.git
cd gdi_pilot
docker build -t goe/gdi-pipeline:latest .
```

## 3. Setup Your Local Environment
Repository structure:

```text
/analysis_directory/
├── config.yaml          # Default configuration for LOCAL execution
├── config_docker.yaml   # Configuration for Docker/Apptainer execution
├── data/
│   ├── cohort.vcf.gz    # Your input VCF (filename is set in the config file)
│   └── samples.tsv      # Metadata (SampleID, Sex)
├── resources/           # Folder for resources (BED targets file, PGx genes and other external resources)
└── results/             # Output directory
```

## 4. Configuration
The pipeline uses two configuration files:
*   **`config.yaml`**: The default configuration for **local execution**, which uses relative paths (e.g., `data/...`).
*   **`config_docker.yaml`**: The configuration for **Docker/Apptainer execution**, which uses absolute paths inside the container (e.g., `/data/...`).

When running in a container, you must mount your `config_docker.yaml` to `/pipeline/config.yaml` to override the default.

**Important Settings in `config.yaml`:**
* **`country_code`**: Country code for stratified AF (default: `PT`, oubviously change for your ).
*  **`input_vcf`:** Path to your input VCF file (default: `data/cohort.vcf.gz`).*  ** `samples_metadata`:** Path to your sample metadata TSV (default: `data/samples.tsv`).

**Other settings:**
*  **`genome_build`:** Specify your input build (default: `GRCh38`).
*  **Target Selection:** Specify `regions_bed` (default: `resources/targets.bed`).

## 5. Pipeline Workflows

This repository contains three specialized workflows.

### A. GoE Pipeline (`Snakefile`)
The primary pipeline for generating standardized frequency data and high-quality internal datasets.

**Architecture:**
(See `docs/graphs/main_workflow.dot` for visual representation)
1.  **Select & Normalization:** Filters variants to target regions and normalizes indels/multiallelics.
2.  **Raw Stats:** Calculates initial population statistics.
3.  **Genotype QC:** Masks genotypes based on quality thresholds:
    *   `GQ < 20`: Low Genotype Quality.
    *   `DP < 10`: Low Depth.
    *   `AB < 0.2`: Allele Balance deviation (Alt / Total Reads < 0.2).
4.  **Final Stats:** Recalculates statistics on the masked data.
5.  **Variant QC:** Tags variants based on site-level metrics:
    *   `QD < 2.0`: Quality by Depth.
    *   `FS > 60.0`: Fisher Strand Bias.
    *   `MQ < 40.0`: RMS Mapping Quality.
    *   `ReadPosRankSum < -8.0`: Read Position Bias.
6.  **Sites Output:** Generates sites-only VCFs (All & PASS only).

**Run Command:**
```bash
snakemake -s Snakefile -j 8
```

### B. Beacon Pipeline (`Snakefile.AF_bcftools_pipeline`)
A wrapper around the EGA/Beacon Nextflow pipeline.

**Run Command:**
```bash
snakemake -s Snakefile.AF_bcftools_pipeline -j 8
```

### C. PyPGX Pipeline (`Snakefile.pypgx`)
A specialized workflow for Pharmacogenomics (PGx) calling using the `PyPGX` tool suite.

**Workflow:**
1.  **Engine:** Runs the `run-ngs-pipeline` on the input VCF.
2.  **Calculation:** Computes star-alleles and phenotypes for target genes.
3.  **Merging:** Aggregates results into a single CSV report.

**Run Command:**
```bash
snakemake -s Snakefile.pypgx -j 8
```

## 6. Outputs
Results are written to your local `results/` folder. The exact outputs depend on which pipeline was run.

### Main Output Files
| File | Pipeline | Description |
|------|----------|-------------|
| `final_annotated.sites.vcf.gz` | **GoE & Beacon** | **Submission File.** A sites-only VCF containing the calculated allele frequencies (**AF**). This is the primary output for data sharing. |
| `final_annotated.vcf.gz` | **GoE only** | The full VCF with all genotypes and the calculated allele frequencies (**AF**). This file should be kept private. |
| `results/pypgx/` | **PyPGX only** | Directory containing PGx star-allele calls and phenotype predictions. |

### Calculated Statistics (GoE Pipeline)
The following statistics are added to the VCF INFO field:

**Metrics:**
*   `AC`: Allele Count
*   `AN`: Allele Number
*   `AF`: Allele Frequency
*   `AC_Het`: Count of Heterozygous genotypes
*   `AC_Homo`: Count of Homozygous genotypes
*   `AC_Hemi`: Count of Hemizygous genotypes

**Stratification:**
Statistics are stratified by Country (e.g., `_PT`) and Sex (`_M`, `_F`).
*   *Example:* `AF_PT_M` (Allele Frequency for Males in the Portuguse cohort).


