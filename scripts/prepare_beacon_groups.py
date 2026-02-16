import pandas as pd
import argparse
import sys

def generate_nextflow_groups(input_tsv, output_csv, country_code="Unknown"):
    try:
        with open(input_tsv, 'r') as f:
            first_line = f.readline()
            has_header = any(keyword in first_line.lower() for keyword in ["sample", "sex", "id", "country"])
        
        if has_header:
            df = pd.read_csv(input_tsv, sep='	')
            df.columns = [c.lower().strip() for c in df.columns]
            sid_col = [c for c in df.columns if 'id' in c or 'sample' in c][0]
            sex_col = [c for c in df.columns if 'sex' in c][0]
        else:
            df = pd.read_csv(input_tsv, sep='	', header=None, names=['sample_id', 'sex', 'country_code'])
            sid_col, sex_col = 'sample_id', 'sex'
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Nextflow pipeline expects 1 for Male, 2 for Female
    sex_map = {'M':'1','F':'2','m':'1','f':'2','1':'1','2':'2', 'MALE':'1', 'FEMALE':'2'}
    
    out_df = pd.DataFrame()
    out_df['SAMPLE'] = df[sid_col].astype(str).str.strip()
    out_df['SEX'] = df[sex_col].astype(str).str.strip().map(sex_map).fillna('0')
    
    # Use country_code from the input file if available, otherwise use the provided country_code argument
    if 'country_code' in df.columns:
        out_df['ANCESTRY'] = df['country_code'].astype(str).str.strip()
    else:
        out_df['ANCESTRY'] = country_code
    
    out_df.to_csv(output_csv, index=False)
    print(f"Created Nextflow groups file: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Nextflow groups for Beacon pipeline.")
    parser.add_argument("input", help="Input TSV file (sample_id, sex, country_code)")
    parser.add_argument("output", help="Output CSV file")
    parser.add_argument("--country", default="Unknown", help="Country code for groups")
    
    args = parser.parse_args()
    generate_nextflow_groups(args.input, args.output, args.country)
