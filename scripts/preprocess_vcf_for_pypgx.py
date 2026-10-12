"""Preprocess VCF files for compatibility with PyPGX.

PyPGX run-ngs-pipeline crashes on VCF files when it encounters comma-separated
missing values in format fields (e.g., AD=".,." for masked/missing genotypes).
This script normalizes any all-missing AD vector to a single "." using HTSlib
vector-end encoding.
"""

import argparse

import numpy as np
from cyvcf2 import VCF, Writer

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

    has_ad = "AD" in vcf

    for variant in vcf:
        if has_ad:
            try:
                ad = variant.format("AD")
            except KeyError:
                ad = None

            if ad is not None:
                if ad.ndim == 1:
                    ad = ad.reshape(1, -1)
                # Detect rows where any element is missing. PyPGX's underlying fuc
                # parser attempts int(x) on comma-separated values and crashes with
                # ValueError on '.' (e.g., AD='.,.' or AD='.,12'). Sanitizing any
                # vector with missing values to a single '.' ensures seamless execution.
                has_missing = np.any(ad == BCF_INT32_MISSING, axis=1)
                if has_missing.any():
                    ad[has_missing, 0] = BCF_INT32_MISSING
                    ad[has_missing, 1:] = BCF_INT32_VECTOR_END
                    variant.set_format("AD", ad)
        writer.write_record(variant)

    writer.close()
    vcf.close()


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
