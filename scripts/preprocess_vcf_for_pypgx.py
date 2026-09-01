"""Preprocess VCF files for compatibility with PyPGX.

PyPGX run-ngs-pipeline crashes on VCF files when it encounters comma-separated
missing values in format fields (e.g., AD=".,." for masked/missing genotypes).
This script normalizes any all-missing AD vector to a single "." using HTSlib
vector-end encoding.
"""

import argparse
import subprocess
from cyvcf2 import VCF, Writer
import numpy as np

BCF_INT32_MISSING = -2147483648
BCF_INT32_VECTOR_END = -2147483647


def preprocess_vcf_for_pypgx(input_vcf: str, output_vcf: str) -> None:
    """Normalize missing AD format values for PyPGX compatibility.

    Parameters
    ----------
    input_vcf : str
        Path to the input VCF/BCF file.
    output_vcf : str
        Path to the output VCF file (should end with .vcf.gz).
    """
    vcf = VCF(input_vcf)
    writer = Writer(output_vcf, vcf)

    for variant in vcf:
        ad = variant.format("AD")
        if ad is not None:
            # Detect rows where every value is missing or vector-end
            fully_missing = np.all(
                (ad == BCF_INT32_MISSING) | (ad == BCF_INT32_VECTOR_END), axis=1
            )
            if fully_missing.any():
                ad[fully_missing, 0] = BCF_INT32_MISSING
                ad[fully_missing, 1:] = BCF_INT32_VECTOR_END
                variant.set_format("AD", ad)
        writer.write_record(variant)

    writer.close()
    vcf.close()

    # Index the output VCF file with tabix
    subprocess.run(["tabix", "-f", "-p", "vcf", output_vcf], check=True)


def main() -> None:
    """Parse arguments and run VCF preprocessing."""
    parser = argparse.ArgumentParser(
        description="Preprocess VCF format fields for PyPGX compatibility."
    )
    parser.add_argument("input_vcf", help="Input VCF file path")
    parser.add_argument("output_vcf", help="Output VCF file path (.vcf.gz)")
    args = parser.parse_args()

    preprocess_vcf_for_pypgx(args.input_vcf, args.output_vcf)


if __name__ == "__main__":
    main()
