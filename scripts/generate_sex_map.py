import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Generate bcftools-compatible sex map file.")
    parser.add_argument("input", help="Input sample info TSV file")
    parser.add_argument("output", help="Output sex map file")
    args = parser.parse_args()
    
    try:
        # Detect header
        with open(args.input, 'r') as f:
            first_line = f.readline()
            has_header = any(keyword in first_line.lower() for keyword in ["sample", "sex", "id", "country"])
        
        if has_header:
            df = pd.read_csv(args.input, sep='\t')
            df.columns = [c.lower().strip() for c in df.columns]
        else:
            df = pd.read_csv(args.input, sep='\t', header=None, names=['sample_id', 'sex', 'country_code'])
            
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
