import pandas as pd
import argparse
import sys

def generate_genomics_metadata(input_tsv, groups_output, suffix=""):
    try:
        # Read the first few lines to detect header
        with open(input_tsv, 'r') as f:
            first_line = f.readline()
            has_header = any(keyword in first_line.lower() for keyword in ["sample", "sex", "id", "country"])
        
        if has_header:
            df = pd.read_csv(input_tsv, sep='\t')
            df.columns = [c.lower().strip() for c in df.columns]
        else:
            # Assume order: sample_id, sex, country_code
            df = pd.read_csv(input_tsv, sep='\t', header=None, names=['sample_id', 'sex', 'country_code'])
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # 1. Generate Group Mapping File (for +fill-tags -S)
    if suffix and not suffix.startswith('_'):
        suffix = f"_{suffix}"

    group_data = []
    for _, row in df.iterrows():
        sample = str(row['sample_id']).strip()
        raw_sex = str(row['sex']).upper().strip()
        
        sex = ""
        if raw_sex in ['M', 'MALE', '1', 'XY']: sex = 'M'
        elif raw_sex in ['F', 'FEMALE', '2', 'XX']: sex = 'F'
        
        country = str(row['country_code']).upper().strip()
        
        # Build comma-separated groups: Country, Sex, Country_Sex
        groups = []
        if country:
            groups.append(f"{country}{suffix}")
        if sex:
            groups.append(f"{sex}{suffix}")
        if country and sex:
            groups.append(f"{country}_{sex}{suffix}")
            
        if groups:
            group_data.append([sample, ",".join(groups)])

    groups_df = pd.DataFrame(group_data)
    groups_df.to_csv(groups_output, sep='\t', index=False, header=False)
    print(f"Created group mapping file: {groups_output} with suffix '{suffix}'")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate bcftools group mapping files.")
    parser.add_argument("input", help="Input TSV file (sample_id, sex, country_code)")
    parser.add_argument("groups", help="Output group mapping file")
    parser.add_argument("--suffix", default="", help="Suffix for group tags (e.g., raw)")
    
    args = parser.parse_args()
    generate_genomics_metadata(args.input, args.groups, args.suffix)
