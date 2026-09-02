# PGX-Pilot Pipeline Rules & Precedence

This document outlines the execution flow and dependencies of the main `Snakefile` for the PGx-Pilot pipeline.

## Main Execution Flow

The pipeline follows a linear progression from initial data ingestion to QC-tagged site-specific VCF generation.

### 1. Target Rule (`all`)
Defines the final outputs:
- `results/{prefix}.sites.all.vcf.gz`
- `results/{prefix}.sites.pass.vcf.gz`
- `results/intermediate/{prefix}.full_sample_data.vcf.gz`

### 2. Site VCF Generation (`create_sites_vcf`)
Generates site-only VCFs (all sites and PASS-only sites).
- **Input:** `results/intermediate/{prefix}.full_sample_data.vcf.gz`

### 3. QC Tagging (`variant_qc_tagging`)
Applies QC filters and tags variants using `scripts/tag_variant_qc.py`.
- **Input:** `results/temp/06_final_stats.vcf.gz`
- **Precedence:** Depends on `annotate_final_vcf`.

### 4. Final Annotation (`annotate_final_vcf`)
Fills variant tags (AC, AF, etc.) using `bcftools +fill-tags`.
- **Input:** `results/temp/05_ploidy_fixed.vcf.gz`, `results/temp/groups_final.txt`
- **Output:** `results/temp/06_final_stats.vcf.gz`
- **Precedence:** Depends on `fix_ploidy` and `generate_groups_final`.

### 5. Ploidy Correction (`fix_ploidy`)
Adjusts male sex-chromosome non-PAR ploidy (hemizygous) and masks spurious female ChrY calls.
- **Input:** `results/temp/04_masked.vcf.gz`, `data/samples.tsv`
- **Output:** `results/temp/05_ploidy_fixed.vcf.gz`
- **Precedence:** Depends on `genotype_masking`.

### 6. Genotype Masking (`genotype_masking`)
Sets low-quality genotypes to missing based on GQ, DP, and AB ratio thresholds.
- **Input:** `results/temp/03_raw_stats.vcf.gz`
- **Output:** `results/temp/04_masked.vcf.gz`
- **Precedence:** Depends on `annotate_raw_vcf`.

### 6. Raw Annotation (`annotate_raw_vcf`)
Calculates initial population statistics on the normalized VCF.
- **Input:** `results/temp/02_normalized.vcf.gz`, `results/temp/groups_raw.txt`
- **Precedence:** Depends on `normalize_and_split` and `generate_groups_raw`.

### 7. Normalization (`normalize_and_split`)
Normalizes variants and splits multi-allelic sites using the reference genome.
- **Input:** `results/temp/01_selected.vcf.gz`, `resources/GRCh38.fa`
- **Precedence:** Depends on `select_regions` and `prepare_reference`.

### 8. Region Selection (`select_regions`)
Subsets the input VCF to specific regions and excludes symbolic alleles.
- **Input:** Configured `input_vcf` and `regions_bed`.
- **Precedence:** Depends on `download_selection_bed` if the BED is a URL.

---

## Resource & Support Rules

- **`index_vcf`**: Generic rule to index any `.vcf.gz` file with `tabix`.
- **`prepare_reference`**: Ensures the GRCh38 reference FASTA and dictionary are available.
- **`download_reference_source`**: Downloads the reference if not found locally.
- **`download_selection_bed`**: Downloads the target regions BED file if provided as a URL.
- **`generate_groups_raw` / `generate_groups_final`**: Creates group sample files for `bcftools +fill-tags`.

## Simplified Precedence Summary
`select_regions` → `normalize_and_split` → `annotate_raw_vcf` → `genotype_masking` → `fix_ploidy` → `annotate_final_vcf` → `variant_qc_tagging` → `create_sites_vcf` → `all`
