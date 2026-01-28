#!/bin/bash
set -e

# --- Configuration ---
BASE_DIR="data/test"
NUM_SAMPLES=5
OUTPUT_DIR_GRCH37="$BASE_DIR/grch37"
OUTPUT_DIR_GRCH38="$BASE_DIR/grch38"
PROCESS_GRCH37=false
NUM_JOBS=4

# 1000 Genomes URLs
# GRCh37: Phase 3
GRCH37_BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
# GRCh38: 30x High Coverage (20201028 release)
GRCH38_BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"
# Pedigree info
PEDIGREE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt"

# --- Usage ---
usage() {
    echo "Usage: $0 [-j JOBS] [--include-37]"
    echo "  -j, --jobs JOBS       Number of parallel download jobs (default: 4)"
    echo "  --include-37          Also process GRCh37 (default: GRCh38 only)"
    echo "  -h, --help            Show this help message"
    exit 1
}

# --- Argument Parsing ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -j|--jobs) NUM_JOBS="$2"; shift ;; 
        --include-37) PROCESS_GRCH37=true ;; 
        -h|--help) usage ;; 
        *) echo "Unknown parameter: $1"; usage ;; 
    esac
    shift
done

mkdir -p "$OUTPUT_DIR_GRCH38"
if [ "$PROCESS_GRCH37" = true ]; then
    mkdir -p "$OUTPUT_DIR_GRCH37"
fi

# --- Dependency Check ---
for cmd in bcftools wget xargs; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found."
        exit 1
    fi
done

# --- Helper Functions ---
# Export variables needed by subshells in xargs
export BASE_DIR GRCH37_BASE_URL GRCH38_BASE_URL OUTPUT_DIR_GRCH37 OUTPUT_DIR_GRCH38

process_grch38_chr() {
    local CHR=$1
    echo "      [38] Processing Chr$CHR..."

    # GRCh38 30x naming pattern (recalibrated variants)
    local VCF_FILENAME="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.recalibrated_variants.vcf.gz"
    local URL="$GRCH38_BASE_URL/$VCF_FILENAME"
    local OUT_PART="$OUTPUT_DIR_GRCH38/chr${CHR}.part.vcf.gz"

    # Download, subset samples, filter, and output
    bcftools view "$URL" \
        --samples-file "$BASE_DIR/sample_list.txt" \
        --force-samples \
        --min-ac 1 \
        -O z -o "$OUT_PART"

    bcftools index "$OUT_PART"
}
export -f process_grch38_chr

process_grch37_chr() {
    local CHR=$1
    echo "      [37] Processing Chr$CHR..."

    local VCF_FILENAME=""
    if [ "$CHR" == "X" ]; then
        VCF_FILENAME="ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
    elif [ "$CHR" == "Y" ]; then
        VCF_FILENAME="ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
    elif [ "$CHR" == "MT" ]; then
        VCF_FILENAME="ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
    else
        VCF_FILENAME="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    fi

    local URL="$GRCH37_BASE_URL/$VCF_FILENAME"
    local OUT_PART="$OUTPUT_DIR_GRCH37/chr${CHR}.part.vcf.gz"

    bcftools view "$URL" \
        --samples-file "$BASE_DIR/sample_list.txt" \
        --force-samples \
        --min-ac 1 \
        -O z -o "$OUT_PART"

    bcftools index "$OUT_PART"
}
export -f process_grch37_chr

# --- Step 1: Metadata ---
echo "[1/3] Preparing Metadata ($NUM_SAMPLES samples)..."
rm -f "$BASE_DIR/samples.tsv" "$BASE_DIR/sample_list.txt"
wget -q -O "$BASE_DIR/1kGP.3202_samples.pedigree_info.txt" "$PEDIGREE_URL"

# Parse PED file:
# Column 2: Sample ID
# Column 5: Sex (1=Male, 2=Female) -> Convert to M/F/U
awk 'NR>1 { if ($4==1) s="M"; else if ($4==2) s="F"; else s="U"; print $1"\t"s }' "$BASE_DIR/1kGP.3202_samples.pedigree_info.txt" | head -n "$NUM_SAMPLES" > "$BASE_DIR/samples.tsv"
cut -f1 "$BASE_DIR/samples.tsv" > "$BASE_DIR/sample_list.txt"
echo "      Created metadata."

# --- Step 2: GRCh38 Whole Genome ---
echo "[2/3] Processing GRCh38 (30x) - Whole Genome ($NUM_JOBS parallel jobs)..."
CHROMS_38=$(seq 1 22)
CHROMS_38="$CHROMS_38 X Y"

# Use printf + xargs for safe parallelism handling of the list
printf "%s\n" $CHROMS_38 | xargs -P "$NUM_JOBS" -I {} bash -c 'process_grch38_chr "$@"' _ {}

echo "      Concatenating GRCh38..."
FILES_TO_CONCAT_38=""
for CHR in $CHROMS_38; do FILES_TO_CONCAT_38="$FILES_TO_CONCAT_38 $OUTPUT_DIR_GRCH38/chr${CHR}.part.vcf.gz"; done

bcftools concat $FILES_TO_CONCAT_38 -O z -o "$OUTPUT_DIR_GRCH38/cohort.vcf.gz"
bcftools index -t "$OUTPUT_DIR_GRCH38/cohort.vcf.gz"
rm $FILES_TO_CONCAT_38 "$OUTPUT_DIR_GRCH38"/*.part.vcf.gz.csi
echo "      Saved $OUTPUT_DIR_GRCH38/cohort.vcf.gz"

# --- Step 3: GRCh37 Whole Genome (Optional) ---
if [ "$PROCESS_GRCH37" = true ]; then
    echo "[3/3] Processing GRCh37 (Phase 3) - Whole Genome ($NUM_JOBS parallel jobs)..."
    CHROMS_37=$(seq 1 22)
    CHROMS_37="$CHROMS_37 X Y MT"

    printf "%s\n" $CHROMS_37 | xargs -P "$NUM_JOBS" -I {} bash -c 'process_grch37_chr "$@"' _ {}

    echo "      Concatenating GRCh37..."
    FILES_TO_CONCAT_37=""
    for CHR in $CHROMS_37; do FILES_TO_CONCAT_37="$FILES_TO_CONCAT_37 $OUTPUT_DIR_GRCH37/chr${CHR}.part.vcf.gz"; done

    bcftools concat $FILES_TO_CONCAT_37 -O z -o "$OUTPUT_DIR_GRCH37/cohort.vcf.gz"
    bcftools index -t "$OUTPUT_DIR_GRCH37/cohort.vcf.gz"
    rm $FILES_TO_CONCAT_37 "$OUTPUT_DIR_GRCH37"/*.part.vcf.gz.csi
    echo "      Saved $OUTPUT_DIR_GRCH37/cohort.vcf.gz"
else
    echo "[3/3] Skipping GRCh37 (use --include-37 to enable)."
fi

echo ""
echo "Done! Full genome test VCFs created."
