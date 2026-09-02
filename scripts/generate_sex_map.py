"""Generate bcftools-compatible sex map file for ploidy fixing."""

import argparse
import sys
import pandas as pd


def generate_sex_map(input_tsv: str, output_map: str) -> None:
    """Generate a sample_id to sex (M/F) mapping file for bcftools +fixploidy.

    Parameters
    ----------
    input_tsv : str
        Path to input sample info TSV file.
    output_map : str
        Path to output sex map file.
    """
    try:
        # Robust header detection: inspect exact keywords of the first columns
        with open(input_tsv, "r") as f:
            first_line = f.readline()
        first_fields = [f.strip().lower() for f in first_line.split("\t")]
        has_header = (
            (len(first_fields) > 0 and first_fields[0] in {"sample_id", "sample", "id", "sampleid", "samples"})
            or (len(first_fields) > 1 and first_fields[1] in {"sex", "gender"})
        )

        if has_header:
            df = pd.read_csv(input_tsv, sep="\t")
            df.columns = [c.lower().strip() for c in df.columns]
        else:
            # Position-based fallback: col 0 is sample_id, col 1 is sex
            df = pd.read_csv(input_tsv, sep="\t", header=None)
            if len(df.columns) >= 2:
                df.columns = ["sample_id", "sex"] + list(df.columns[2:])
            else:
                raise ValueError("Input file must have at least 2 tab-separated columns.")
            
        # Standardize sex mapping
        sex_map = {
            'M': 'M', 'MALE': 'M', '1': 'M', 'XY': 'M',
            'F': 'F', 'FEMALE': 'F', '2': 'F', 'XX': 'F'
        }
        df['sex_clean'] = df['sex'].astype(str).str.upper().str.strip().map(sex_map)
        
        # Verify all samples are mapped successfully
        if df['sex_clean'].isna().any():
            invalid = df[df['sex_clean'].isna()]['sample_id'].tolist()
            print(f"Warning: The following samples could not be mapped to M or F: {invalid}")
            
        df[["sample_id", "sex_clean"]].to_csv(output_map, sep="\t", index=False, header=False)
        print(f"Successfully generated sex map: {output_map}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


def main() -> None:
    """Parse CLI arguments and generate sex map."""
    parser = argparse.ArgumentParser(description="Generate bcftools-compatible sex map file.")
    parser.add_argument("input", help="Input sample info TSV file")
    parser.add_argument("output", help="Output sex map file")
    args = parser.parse_args()

    generate_sex_map(args.input, args.output)


if __name__ == "__main__":
    main()
