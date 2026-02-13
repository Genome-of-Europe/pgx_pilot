import pandas as pd
import zipfile
import io
import csv
import argparse
import pypgx
import os

def calculate_stats(zip_path, groups_path, gene, out_alleles, out_genotypes, out_phenotypes):
    # 1. Load sample to group mapping
    sample_groups = {}
    all_groups = set(["ALL"])
    with open(groups_path, 'r') as f:
        for line in f:
            parts = line.strip().split('	')
            if len(parts) < 2: continue
            sid, groups_str = parts[0], parts[1]
            groups = groups_str.split(',')
            sample_groups[sid] = groups + ["ALL"]
            for g in groups:
                all_groups.add(g)
    
    all_groups = sorted(list(all_groups))

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        tsv_path = next((f for f in zip_ref.namelist() if f.endswith('data.tsv')), None)
        if not tsv_path:
            print(f"Error: No data.tsv found in {zip_path}")
            return

        with zip_ref.open(tsv_path) as tsv_file:
            df = pd.read_csv(io.StringIO(tsv_file.read().decode('utf-8')), sep="\t")

        # Rename the first column to 'Sample' if it's unnamed (PyPGX format)
        if 'Sample' not in df.columns and 'Unnamed: 0' in df.columns:
            df.rename(columns={'Unnamed: 0': 'Sample'}, inplace=True)
        elif 'Sample' not in df.columns:
            # Fallback: assume the first column is the sample ID
            df.rename(columns={df.columns[0]: 'Sample'}, inplace=True)
        
        # Prepare output files
        with (
            open(out_alleles, "w", newline="") as f_alleles,
            open(out_genotypes, "w", newline="") as f_genotypes,
            open(out_phenotypes, "w", newline="") as f_phenotypes,
        ):
            
            w_alleles = csv.writer(f_alleles)
            w_genotypes = csv.writer(f_genotypes)
            w_phenotypes = csv.writer(f_phenotypes)
            
            w_alleles.writerow(["Gene", "Group", "Allele", "Function", "AC", "AN", "AF"])
            w_genotypes.writerow(["Gene", "Group", "Genotype", "Phenotype", "Count", "Frequency"])
            w_phenotypes.writerow(["Gene", "Group", "Phenotype", "Count", "Frequency"])

            for group in all_groups:
                # Subset samples for this group
                if group == "ALL":
                    subset_df = df
                else:
                    group_samples = [sid for sid, groups in sample_groups.items() if group in groups]
                    subset_df = df[df['Sample'].isin(group_samples)]
                
                if subset_df.empty:
                    continue
                
                total_samples = len(subset_df)
                an = total_samples * 2

                # --- Phenotypes ---
                pheno_counts = subset_df["Phenotype"].value_counts().reset_index()
                # Handle pandas version differences in value_counts output
                if 'index' in pheno_counts.columns:
                    pheno_counts.columns = ['Phenotype', 'count']
                
                for _, row in pheno_counts.iterrows():
                    w_phenotypes.writerow([gene, group, row['Phenotype'], row['count'], row['count']/total_samples])

                # --- Genotypes ---
                gt_counts = subset_df.groupby(["Genotype", "Phenotype"]).size().reset_index(name='count')
                for _, row in gt_counts.iterrows():
                    w_genotypes.writerow([gene, group, row['Genotype'], row['Phenotype'], row['count'], row['count']/total_samples])

                # --- Alleles ---
                all_alleles = []
                for genotype in subset_df["Genotype"]:
                    all_alleles.extend(genotype.replace('|', '/').split('/'))
                
                al_counts = pd.Series(all_alleles).value_counts().reset_index()
                if 'index' in al_counts.columns:
                    al_counts.columns = ['Allele', 'AC']
                else:
                    al_counts.columns = ['Allele', 'count'] # Support newer pandas
                    al_counts.rename(columns={'count': 'AC'}, inplace=True)

                for _, row in al_counts.iterrows():
                    allele = row['Allele']
                    ac = row['AC']
                    try:
                        function = pypgx.get_function(gene, allele)
                    except:
                        function = "Unknown"
                    w_alleles.writerow([gene, group, allele, function, ac, an, ac/an])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out-alleles", required=True)
    parser.add_argument("--out-genotypes", required=True)
    parser.add_argument("--out-phenotypes", required=True)
    args = parser.parse_args()
    
    calculate_stats(args.zip, args.groups, args.gene, args.out_alleles, args.out_genotypes, args.out_phenotypes)
