"""Generate a chromosome rename mapping file for bcftools annotate --rename-chrs."""

import argparse

from cyvcf2 import VCF


def check_for_chr(vcf_path: str) -> bool:
    """Check if any contig in the VCF header starts with 'chr'.

    Parameters
    ----------
    vcf_path : str
        Path to input VCF/BCF file.

    Returns
    -------
    bool
        True if at least one contig name starts with 'chr', False otherwise.

    """
    with VCF(vcf_path) as vcf:
        return any(seq.startswith("chr") for seq in vcf.seqnames)


def main() -> None:
    """Parse arguments and write chromosome rename map."""
    parser = argparse.ArgumentParser(
        description="Generate chromosome rename map from Ensembl to UCSC."
    )
    parser.add_argument("input", help="Input VCF path")
    parser.add_argument("output", help="Output map file path")
    args = parser.parse_args()

    with open(args.output, "w") as f:
        has_chr = check_for_chr(args.input)
        if not has_chr:
            # Map Ensembl numerals to UCSC 'chr' prefix
            for i in range(1, 23):
                f.write(f"{i}\tchr{i}\n")
            f.write("X\tchrX\n")
            f.write("Y\tchrY\n")
            f.write("MT\tchrM\n")
        else:
            # No-op (empty map file signals bcftools no renaming is necessary)
            pass


if __name__ == "__main__":
    main()
