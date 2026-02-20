# Beacon Allele Frequency Pipeline

This workflow is a wrapper around the EGA/Beacon Nextflow pipeline for calculating stratified allele frequencies.

## 1. Overview
The Beacon pipeline calculates Allele Frequencies (AF) stratified by Country and Sex, following the standards of the EGA/Beacon project.

## 2. Setup
Ensure you have the repository structure as described in the main `README.md`.

## 3. Configuration
The Beacon pipeline uses its own configuration file: `config.beacon.yaml`.

**Important Settings in `config.beacon.yaml`:**
* **`output_prefix`**: Prefix for final VCF files (e.g., `beacon_pilot`).
* **`input_vcf`**: Path to your input VCF file.
* **`sample_info`**: Path to your sample metadata TSV.
* **`country_code`**: ISO 3166-1 alpha-2 code for your cohort.

## 4. Running the Pipeline
You can run the pipeline using snakemake:

```bash
snakemake -s Snakefile.AF_bcftools_pipeline -j 8
```

Using docker:
```bash
docker run --rm -v $(pwd):/pipeline goe/pgx-pipeline:latest snakemake -s Snakefile.AF_bcftools_pipeline -j 8
```

## 5. Outputs
Results are written to the `results/` folder with the prefix specified in `config.beacon.yaml`.
* `results/{output_prefix}.sites.pass.vcf.gz`: Stratified frequency data in Beacon format.
