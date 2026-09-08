# GoE PGx/GDI MAP Stage 1 pipelines

## 1. Overview
This repository contains snakemake workflow files for the analysis of GoE legacy data as part of the GoE PGx pilot and the GDI MAP Stage 2.

These perform the following: 

1. Filtering variants based on specified genomic regions (selected pharmacogenes).
2. Splitting multiallelic sites and left-aligning indels.
3. Genotype masking and variant filtering.
4. Calculating Allele Frequencies (AF) stratified by Country and Sex.
5. Generating a sites-only VCF, compatible with the format expected for upload to Beacon.
6. Additionally, generating star allele and phenotype (drug response) frequencies, stratified in the same groups as the above.

The goal is to allow partners to run the exact same workflow on their local data and produce standardized, comparable outputs.
The workflow requires input data aligned to the **GRCh38 (hg38)** reference genome (both `chr1` and `1` chromosome naming formats are automatically detected and supported).

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
docker build -t goe/pgx_pilot:latest .
```

#### Using Apptainer / Singularity
Pull or build (as above) the Docker image:
```bash
docker pull ghcr.io/genome-of-europe/pgx_pilot:latest
docker tag ghcr.io/genome-of-europe/pgx_pilot:latest goe/pgx_pilot:latest
```


Convert the docker image to a Singularity image file (`.sif`):
```bash
apptainer build pgx_pilot.sif docker-daemon://goe/pgx_pilot:latest
```

## 3. Setup Your Local Environment
Repository structure:

```text
/analysis_directory/
├── config.yaml           # Pipeline configuration
├── data/
│   ├── <your_cohort>.vcf.gz     # Your input VCF (bgzipped)
│   ├── <your_cohort>.vcf.gz.tbi # Tabix index file for input VCF
│   └── samples.tsv              # Metadata (sample_id, sex, country_code)
├── resources/            # Folder for resources
│   └── target_genes.txt  # Target pharmacogenes list
└── results/              # Output directory
```

## 4. Configuration
The pipeline uses a single configuration file: `config.yaml`.

**Important Settings in `config.yaml`:**
* **`output_prefix`**: Prefix for final VCF files (e.g., `pgx_pilot`).
* **`input_vcf`**: Path to your input VCF file (e.g., `data/<your_cohort>.vcf.gz`). The corresponding tabix index (`<your_cohort>.vcf.gz.tbi`) must be located in the same directory.
* **`sample_info`**: Path to your sample metadata TSV.
* **`keep_intermediate`**: Set to `false` (default) to automatically remove intermediate VCFs once downstream rules complete, saving disk space. Set to `true` to retain all intermediate VCF files.

The sample metadata file must be a three column TSV file, with or without header, where the columns are:
* **sample_id**: Sample ID (same as the sample ID in the VCF file)
* **sex**: one of (M, F, 1, 2, Male, Female, XY, XX). Case-insensitive.
* **country_code**: the two letter country code (ISO 3166-1 alpha-2 code). Case-insensitive.

**QC Thresholds:**
The pipeline uses standardized thresholds defined in the `qc_thresholds` section of `config.yaml`:

## 5. Pipeline Workflows

This repository contains two specialized workflows.

### A. AF Pipeline (`Snakefile`)
The primary pipeline for generating standardized allele frequency data.

**Run Command:**
using conda
```bash
snakemake -s Snakefile --cores 4
```
using docker
```bash
docker run --rm -v $(pwd):/pipeline goe/pgx_pilot:latest snakemake -s Snakefile --cores 4
```
using Apptainer
```bash
apptainer exec --bind $(pwd):/pipeline pgx_pilot.sif snakemake -s Snakefile --cores 4
```

### B. PyPGX Pipeline (`Snakefile.pypgx`)
A specialized workflow for Pharmacogenomics (PGx) calling using the `PyPGX` tool suite.

**Run Command:**
using conda
```bash
snakemake -s Snakefile.pypgx --cores 8
```
using docker
```bash
docker run --rm -v $(pwd):/pipeline goe/pgx_pilot:latest snakemake -s Snakefile.pypgx --cores 8
```
using Apptainer
```bash
apptainer exec --bind $(pwd):/pipeline pgx_pilot.sif snakemake -s Snakefile.pypgx --cores 8
```

## 6. Outputs
Results are written to the `results/` folder.

### Main Output Files
*   **GoE Pipeline Outputs**
    *   `results/{output_prefix}.sites.pass.vcf.gz`: **Submission File for beacon.** PASS variants only.
    *   `results/{output_prefix}.sites.all.vcf.gz`: **For GoE PGx pilot (Internal).** All variants including those that failed QC.
    *   `results/{output_prefix}.full_sample_data.vcf.gz`: **The full VCF with all sample genotypes.** All variants including those that failed QC.

*   **PyPGX Pipeline Outputs (Shareable Summaries)**
    *   `results/pgx/merged_alleles.csv`: Aggregated star-allele calls across genes, stratified by Group.
    *   `results/pgx/merged_phenotypes.csv`: Aggregated phenotype predictions, stratified by Group.
    *   `results/pgx/merged_genotypes.csv`: Aggregated genotype calls and frequencies, stratified by Group.

### Calculated Statistics
Statistics (AC, AN, AF, etc.) are stratified by Country and Sex based on the input metadata.
*   *Example:* `AF_PT_M` (Allele Frequency for Males in the PT cohort).
