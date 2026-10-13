"""Extract sample-level pharmacogenomic calls from PyPGX results."""

import argparse
import io
import os
import zipfile
from typing import Tuple

import pandas as pd
import pypgx


def extract_sample_dfs(
    zip_path: str,
    gene: str,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Extract sample-level DataFrames from a PyPGX results.zip archive.

    Parameters
    ----------
    zip_path : str
        Path to PyPGX results.zip archive.
    gene : str
        Target gene symbol.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Tuple of (df_alleles, df_genotypes, df_phenotypes).

    Raises
    ------
    FileNotFoundError
        If data.tsv is missing from the zip archive.

    """
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        tsv_path = next((f for f in zip_ref.namelist() if f.endswith("data.tsv")), None)
        if not tsv_path:
            raise FileNotFoundError(f"No data.tsv found in '{zip_path}'.")

        with zip_ref.open(tsv_path) as tsv_file:
            df = pd.read_csv(io.StringIO(tsv_file.read().decode("utf-8")), sep="\t")

    # Standardize sample column name
    if "Sample" not in df.columns and "Unnamed: 0" in df.columns:
        df.rename(columns={"Unnamed: 0": "Sample"}, inplace=True)
    elif "Sample" not in df.columns:
        df.rename(columns={df.columns[0]: "Sample"}, inplace=True)

    df["Sample"] = df["Sample"].astype(str)

    # 1. Sample Alleles (long/tidy format: Sample, Gene, Allele, Function)
    allele_rows = []
    for _, row in df.iterrows():
        sample = str(row["Sample"])
        genotype = row.get("Genotype")
        if pd.isna(genotype) or not str(genotype).strip():
            continue

        raw_alleles = str(genotype).replace("|", "/").split("/")
        for allele in raw_alleles:
            allele_clean = allele.strip()
            if not allele_clean:
                continue
            try:
                function = pypgx.get_function(gene, allele_clean)
            except Exception:
                function = "Unknown"
            allele_rows.append({
                "Sample": sample,
                "Gene": gene,
                "Allele": allele_clean,
                "Function": function,
            })

    df_alleles = pd.DataFrame(
        allele_rows,
        columns=["Sample", "Gene", "Allele", "Function"],
    )

    # 2. Sample Genotypes (Sample, Gene, Genotype, Phenotype)
    genotype_rows = []
    for _, row in df.iterrows():
        sample = str(row["Sample"])
        gt = row.get("Genotype")
        genotype = None if pd.isna(gt) else str(gt)
        ph = row.get("Phenotype")
        phenotype = None if pd.isna(ph) else str(ph)
        genotype_rows.append({
            "Sample": sample,
            "Gene": gene,
            "Genotype": genotype,
            "Phenotype": phenotype,
        })

    df_genotypes = pd.DataFrame(
        genotype_rows,
        columns=["Sample", "Gene", "Genotype", "Phenotype"],
    )

    # 3. Sample Phenotypes (Sample, Gene, Phenotype)
    phenotype_rows = []
    for _, row in df.iterrows():
        sample = str(row["Sample"])
        ph = row.get("Phenotype")
        phenotype = None if pd.isna(ph) else str(ph)
        phenotype_rows.append({
            "Sample": sample,
            "Gene": gene,
            "Phenotype": phenotype,
        })

    df_phenotypes = pd.DataFrame(
        phenotype_rows,
        columns=["Sample", "Gene", "Phenotype"],
    )

    return df_alleles, df_genotypes, df_phenotypes


def extract_sample_data(
    zip_path: str,
    gene: str,
    out_alleles: str,
    out_genotypes: str,
    out_phenotypes: str,
) -> None:
    """Extract sample-level pharmacogenomic calls and write to CSV files.

    Parameters
    ----------
    zip_path : str
        Path to PyPGX results.zip archive.
    gene : str
        Target gene symbol.
    out_alleles : str
        Output CSV path for sample star alleles.
    out_genotypes : str
        Output CSV path for sample genotypes.
    out_phenotypes : str
        Output CSV path for sample phenotypes.

    """
    df_alleles, df_genotypes, df_phenotypes = extract_sample_dfs(zip_path, gene)

    # Ensure parent directories exist
    for out_path in (out_alleles, out_genotypes, out_phenotypes):
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    df_alleles.to_csv(out_alleles, index=False)
    df_genotypes.to_csv(out_genotypes, index=False)
    df_phenotypes.to_csv(out_phenotypes, index=False)


def main() -> None:
    """Parse CLI arguments and extract sample-level PGx data."""
    parser = argparse.ArgumentParser(
        description="Extract sample star alleles, genotypes, and phenotypes from PyPGX."
    )
    parser.add_argument("--zip", required=True, help="Path to PyPGX results.zip")
    parser.add_argument("--gene", required=True, help="Target gene symbol")
    parser.add_argument(
        "--out-alleles", required=True, help="Output sample alleles CSV"
    )
    parser.add_argument(
        "--out-genotypes", required=True, help="Output sample genotypes CSV"
    )
    parser.add_argument(
        "--out-phenotypes", required=True, help="Output sample phenotypes CSV"
    )
    args = parser.parse_args()

    extract_sample_data(
        args.zip,
        args.gene,
        args.out_alleles,
        args.out_genotypes,
        args.out_phenotypes,
    )


if __name__ == "__main__":
    main()
