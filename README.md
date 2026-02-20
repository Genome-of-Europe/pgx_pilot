# GoE PGx/GDI MAP Stage 1 pipelines

## 1. Overview
This repository contains snakemake workflow files for the analysis of GoE legacy data as part of the GoE PGx pilot and the GDI MAP Stage 1.

These perform the following: 

1. Filtering variants based on specified genomic regions (selected pharmacogenes).
2. Splitting multiallelic sites and left-aligning indels.
3. Genotype masking and variant filtering.
4. Calculating Allele Frequencies (AF) stratified by Country and Sex.
5. Generating a sites-only VCF, compatible with the format expected for GDI MAP Stage 1.
6. Additionally, generating star allele and phenotype (drug response) frequencies, stratified in the same groups as the above.

The goal is to allow partners to run the exact same workflow on their local data and produce standardized, comparable outputs.

## 2. Quick Start
The pipeline can be run as Docker container, or directly from the code in the GitHub repository.

### Build Instructions
#### Using conda
1. Install [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install).
2. Create and activate conda env
```bash
conda env create -f env.yml 
conda activate pgx_pilot
```

#### Using docker
```bash
git clone https://github.com/Genome-of-Europe/pgx_pilot.git
cd pgx_pilot
docker build -t goe/pgx-pipeline:latest .
```

## 3. Setup Your Local Environment
Repository structure:

```text
/analysis_directory/
├── config.yaml           # Pipeline configuration
├── data/
│   ├── raw_cohort.vcf.gz # Your input VCF (path set in config.yaml)
│   └── samples.tsv       # Metadata (SampleID, Sex, CountryCode)
├── resources/            # Folder for resources
│   ├── targets.bed       # BED file with PGx genes
│   └── pypgx_genes.txt   # Required for PyPGX workflow
└── results/              # Output directory
```

## 4. Configuration
The pipeline uses a single configuration file: `config.yaml`.

**Important Settings in `config.yaml`:**
* **`output_prefix`**: Prefix for final VCF files (e.g., `pgx_pilot`).
* **`input_vcf`**: Path to your input VCF file.
* **`sample_info`**: Path to your sample metadata TSV.

The sample metadata file must be a three column TSV file, with or without header, where the columns are:
* **sample\_id**: Sample ID (same as the sample ID in the VCF file)
* **sex**: one of (M, F, 1, 2, Male, Female, XY, XX). Case-insensitive.
* **country\_code**: the two letter country code (ISO 3166-1 alpha-2 code). Case-insensitive.

**QC Thresholds:**
The pipeline uses standardized thresholds, as defined in the `qc_thresholds` section of `config.yaml`.

## 5. Pipeline Workflows

This repository contains two specialized workflows.

### A. AF Pipeline (`Snakefile`)
The primary pipeline for generating standardized frequency data and high-quality internal datasets.

**Run Command:**
using conda
```bash
snakemake -s Snakefile -j 8
```
using docker
```bash
docker run --rm -v $(pwd):/pipeline goe/pgx-pipeline:latest snakemake -s Snakefile -j 8
```

### B. PyPGX Pipeline (`Snakefile.pypgx`)
A specialized workflow for Pharmacogenomics (PGx) calling using the `PyPGX` tool suite.

**Run Command:**
using conda
```bash
snakemake -s Snakefile.pypgx -j 8
```

using docker
```bash
docker run --rm -v $(pwd):/pipeline goe/pgx-pipeline:latest snakemake -s Snakefile.pypgx -j 8
```

## 6. Outputs
Results are written to the `results/` folder.

### Main Output Files
*   **GoE Pipeline Outputs**
    *   `results/{output_prefix}.sites.pass.vcf.gz`: **Submission File for GDI MAP Stage 1.** PASS variants only.
    *   `results/{output_prefix}.sites.all.vcf.gz`: **For GoE PGx pilot.** All variants including those that failed QC.
    *   `results/intermediate/{output_prefix}.full_sample_data.vcf.gz`: The full VCF with all sample genotypes. **Keep private.**

*   **PyPGX Pipeline Outputs**
    *   `results/pgx/merged_alleles.csv`: Aggregated star-allele calls across genes, stratified by Group.
    *   `results/pgx/merged_phenotypes.csv`: Aggregated phenotype predictions, stratified by Group.
    *   `results/pgx/merged_genotypes.csv`: Aggregated genotype calls and frequencies, stratified by Group.

### Calculated Statistics
Statistics (AC, AN, AF, etc.) are stratified by Country and Sex based on the input metadata.
*   *Example:* `AF_PT_M` (Allele Frequency for Males in the PT cohort).
