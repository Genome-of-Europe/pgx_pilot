"""Calculate stratified pharmacogenomic allele, genotype, and phenotype frequencies."""

import argparse
import csv
import os
from typing import Dict, List, Set, Tuple

import pandas as pd
import pypgx

try:
    from scripts.extract_pgx_samples import extract_sample_dfs
except ModuleNotFoundError:
    from extract_pgx_samples import extract_sample_dfs


def load_sample_groups(groups_path: str) -> Tuple[Dict[str, List[str]], List[str]]:
    """Load sample to group mappings from TSV file.

    Parameters
    ----------
    groups_path : str
        Path to group mapping TSV file.

    Returns
    -------
    Tuple[Dict[str, List[str]], List[str]]
        Mapping of sample_id -> list of groups, and sorted list of all
        unique group names.

    """
    sample_groups: Dict[str, List[str]] = {}
    all_groups: Set[str] = {"ALL"}

    with open(groups_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sid, groups_str = str(parts[0]).strip(), parts[1].strip()
            groups = [g.strip() for g in groups_str.split(",") if g.strip()]
            sample_groups[sid] = groups + ["ALL"]
            for g in groups:
                all_groups.add(g)

    return sample_groups, sorted(list(all_groups))


def aggregate_pgx_stats(
    df_alleles: pd.DataFrame,
    df_genotypes: pd.DataFrame,
    df_phenotypes: pd.DataFrame,
    groups_path: str,
    gene: str,
    out_alleles: str,
    out_genotypes: str,
    out_phenotypes: str,
) -> None:
    """Calculate and export stratified frequencies from sample-level DataFrames.

    Parameters
    ----------
    df_alleles : pd.DataFrame
        Sample-level star alleles DataFrame (Sample, Gene, Allele, Function).
    df_genotypes : pd.DataFrame
        Sample-level genotypes DataFrame (Sample, Gene, Genotype, Phenotype).
    df_phenotypes : pd.DataFrame
        Sample-level phenotypes DataFrame (Sample, Gene, Phenotype).
    groups_path : str
        Path to group mapping TSV file.
    gene : str
        Target gene symbol.
    out_alleles : str
        Output CSV path for allele frequencies.
    out_genotypes : str
        Output CSV path for genotype frequencies.
    out_phenotypes : str
        Output CSV path for phenotype frequencies.

    """
    sample_groups, all_groups = load_sample_groups(groups_path)

    # Ensure sample IDs are string type
    df_alleles = df_alleles.copy()
    df_genotypes = df_genotypes.copy()
    df_phenotypes = df_phenotypes.copy()
    df_alleles["Sample"] = df_alleles["Sample"].astype(str)
    df_genotypes["Sample"] = df_genotypes["Sample"].astype(str)
    df_phenotypes["Sample"] = df_phenotypes["Sample"].astype(str)

    # Ensure parent output directories exist
    for out_path in (out_alleles, out_genotypes, out_phenotypes):
        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    with (
        open(out_alleles, "w", newline="") as f_alleles,
        open(out_genotypes, "w", newline="") as f_genotypes,
        open(out_phenotypes, "w", newline="") as f_phenotypes,
    ):
        w_alleles = csv.writer(f_alleles)
        w_genotypes = csv.writer(f_genotypes)
        w_phenotypes = csv.writer(f_phenotypes)

        w_alleles.writerow(
            ["Gene", "Group", "Allele", "Function", "AC", "AN", "AF"]
        )
        w_genotypes.writerow(
            ["Gene", "Group", "Genotype", "Phenotype", "Count", "Frequency"]
        )
        w_phenotypes.writerow(
            ["Gene", "Group", "Phenotype", "Count", "Frequency"]
        )

        for group in all_groups:
            # Subset samples for this group based on sample roster
            if group == "ALL":
                sub_gt = df_genotypes
                sub_ph = df_phenotypes
                sub_al = df_alleles
            else:
                group_samples = {
                    sid for sid, groups in sample_groups.items() if group in groups
                }
                sub_gt = df_genotypes[df_genotypes["Sample"].isin(group_samples)]
                sub_ph = df_phenotypes[df_phenotypes["Sample"].isin(group_samples)]
                sub_al = df_alleles[df_alleles["Sample"].isin(group_samples)]

            if sub_gt.empty and sub_ph.empty:
                continue

            total_samples = len(sub_gt)

            # --- Phenotypes ---
            pheno_counts = sub_ph["Phenotype"].value_counts().reset_index()
            if "index" in pheno_counts.columns:
                pheno_counts.columns = ["Phenotype", "count"]

            for _, row in pheno_counts.iterrows():
                ph = row["Phenotype"]
                cnt = row["count"]
                freq = (cnt / total_samples) if total_samples > 0 else 0.0
                w_phenotypes.writerow([gene, group, ph, cnt, freq])

            # --- Genotypes ---
            gt_counts = (
                sub_gt.groupby(["Genotype", "Phenotype"])
                .size()
                .reset_index(name="count")
            )
            for _, row in gt_counts.iterrows():
                gt = row["Genotype"]
                ph = row["Phenotype"]
                cnt = row["count"]
                freq = (cnt / total_samples) if total_samples > 0 else 0.0
                w_genotypes.writerow([gene, group, gt, ph, cnt, freq])

            # --- Alleles ---
            an = len(sub_al)
            al_counts = sub_al["Allele"].value_counts().reset_index()
            if "index" in al_counts.columns:
                al_counts.columns = ["Allele", "AC"]
            else:
                al_counts.columns = ["Allele", "count"]
                al_counts.rename(columns={"count": "AC"}, inplace=True)

            for _, row in al_counts.iterrows():
                allele = row["Allele"]
                ac = row["AC"]
                try:
                    function = pypgx.get_function(gene, allele)
                except Exception:
                    function = "Unknown"
                af = (ac / an) if an > 0 else 0.0
                w_alleles.writerow([gene, group, allele, function, ac, an, af])


def calculate_stats_from_samples(
    sample_alleles_path: str,
    sample_genotypes_path: str,
    sample_phenotypes_path: str,
    groups_path: str,
    gene: str,
    out_alleles: str,
    out_genotypes: str,
    out_phenotypes: str,
) -> None:
    """Calculate and export stratified frequencies from sample-level CSV files.

    Parameters
    ----------
    sample_alleles_path : str
        Path to sample star alleles CSV file.
    sample_genotypes_path : str
        Path to sample genotypes CSV file.
    sample_phenotypes_path : str
        Path to sample phenotypes CSV file.
    groups_path : str
        Path to group mapping TSV file.
    gene : str
        Target gene symbol.
    out_alleles : str
        Output CSV path for allele frequencies.
    out_genotypes : str
        Output CSV path for genotype frequencies.
    out_phenotypes : str
        Output CSV path for phenotype frequencies.

    """
    df_alleles = pd.read_csv(
        sample_alleles_path,
        dtype={"Sample": str, "Gene": str, "Allele": str, "Function": str},
    )
    df_genotypes = pd.read_csv(
        sample_genotypes_path,
        dtype={"Sample": str, "Gene": str, "Genotype": str, "Phenotype": str},
    )
    df_phenotypes = pd.read_csv(
        sample_phenotypes_path,
        dtype={"Sample": str, "Gene": str, "Phenotype": str},
    )

    aggregate_pgx_stats(
        df_alleles=df_alleles,
        df_genotypes=df_genotypes,
        df_phenotypes=df_phenotypes,
        groups_path=groups_path,
        gene=gene,
        out_alleles=out_alleles,
        out_genotypes=out_genotypes,
        out_phenotypes=out_phenotypes,
    )


def calculate_stats_from_zip(
    zip_path: str,
    groups_path: str,
    gene: str,
    out_alleles: str,
    out_genotypes: str,
    out_phenotypes: str,
) -> None:
    """Calculate and export stratified frequencies directly from PyPGX results.zip.

    Parameters
    ----------
    zip_path : str
        Path to PyPGX results.zip archive.
    groups_path : str
        Path to group mapping TSV file.
    gene : str
        Target gene symbol.
    out_alleles : str
        Output CSV path for allele frequencies.
    out_genotypes : str
        Output CSV path for genotype frequencies.
    out_phenotypes : str
        Output CSV path for phenotype frequencies.

    """
    df_alleles, df_genotypes, df_phenotypes = extract_sample_dfs(zip_path, gene)

    aggregate_pgx_stats(
        df_alleles=df_alleles,
        df_genotypes=df_genotypes,
        df_phenotypes=df_phenotypes,
        groups_path=groups_path,
        gene=gene,
        out_alleles=out_alleles,
        out_genotypes=out_genotypes,
        out_phenotypes=out_phenotypes,
    )


def calculate_stats(
    zip_path: str,
    groups_path: str,
    gene: str,
    out_alleles: str,
    out_genotypes: str,
    out_phenotypes: str,
) -> None:
    """Backwards-compatible wrapper to calculate PGx stats from PyPGX results.zip.

    Parameters
    ----------
    zip_path : str
        Path to PyPGX results.zip archive.
    groups_path : str
        Path to group mapping TSV file.
    gene : str
        Target gene symbol.
    out_alleles : str
        Output CSV path for allele frequencies.
    out_genotypes : str
        Output CSV path for genotype frequencies.
    out_phenotypes : str
        Output CSV path for phenotype frequencies.

    """
    calculate_stats_from_zip(
        zip_path=zip_path,
        groups_path=groups_path,
        gene=gene,
        out_alleles=out_alleles,
        out_genotypes=out_genotypes,
        out_phenotypes=out_phenotypes,
    )


def main() -> None:
    """Parse CLI arguments and compute PGx statistics."""
    parser = argparse.ArgumentParser(
        description="Calculate stratified pharmacogenomic frequencies."
    )
    parser.add_argument(
        "--zip", required=False, help="Path to PyPGX results.zip (standalone mode)"
    )
    parser.add_argument(
        "--sample-alleles", required=False, help="Path to sample alleles CSV"
    )
    parser.add_argument(
        "--sample-genotypes", required=False, help="Path to sample genotypes CSV"
    )
    parser.add_argument(
        "--sample-phenotypes", required=False, help="Path to sample phenotypes CSV"
    )
    parser.add_argument("--groups", required=True, help="Path to groups mapping file")
    parser.add_argument("--gene", required=True, help="Target gene symbol")
    parser.add_argument("--out-alleles", required=True, help="Output alleles CSV")
    parser.add_argument("--out-genotypes", required=True, help="Output genotypes CSV")
    parser.add_argument("--out-phenotypes", required=True, help="Output phenotypes CSV")
    args = parser.parse_args()

    has_samples = (
        args.sample_alleles and args.sample_genotypes and args.sample_phenotypes
    )
    if not args.zip and not has_samples:
        parser.error(
            "Either --zip OR (--sample-alleles, --sample-genotypes, "
            "--sample-phenotypes) must be provided."
        )

    if has_samples:
        calculate_stats_from_samples(
            sample_alleles_path=args.sample_alleles,
            sample_genotypes_path=args.sample_genotypes,
            sample_phenotypes_path=args.sample_phenotypes,
            groups_path=args.groups,
            gene=args.gene,
            out_alleles=args.out_alleles,
            out_genotypes=args.out_genotypes,
            out_phenotypes=args.out_phenotypes,
        )
    else:
        calculate_stats_from_zip(
            zip_path=args.zip,
            groups_path=args.groups,
            gene=args.gene,
            out_alleles=args.out_alleles,
            out_genotypes=args.out_genotypes,
            out_phenotypes=args.out_phenotypes,
        )


if __name__ == "__main__":
    main()
