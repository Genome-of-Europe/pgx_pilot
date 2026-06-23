import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Generate bcftools-compatible sex map file.")
    parser.add_argument("input", help="Input sample info TSV file")
    parser.add_argument("output", help="Output sex map file")
    args = parser.parse_args()
    
    try:
        # Simple header check: does the first column header equal "sample_id"?
        with open(args.input, 'r') as f:
            first_line = f.readline()
        has_header = first_line.split('\t')[0].lower().strip() == "sample_id"
        
        if has_header:
            df = pd.read_csv(args.input, sep='\t')
            df.columns = [c.lower().strip() for c in df.columns]
        else:
            # Position-based fallback: col 0 is sample_id, col 1 is sex
            df = pd.read_csv(args.input, sep='\t', header=None)
            if len(df.columns) >= 2:
                df.columns = ['sample_id', 'sex'] + list(df.columns[2:])
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
            
        df[['sample_id', 'sex_clean']].to_csv(args.output, sep='\t', index=False, header=False)
        print(f"Successfully generated sex map: {args.output}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
