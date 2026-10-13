"""Check that REF/ALT alleles in a VCF are accepted by Beacon."""

import argparse
import re
import sys
from typing import List, Tuple

from cyvcf2 import VCF

# Pattern Beacon v2 applies to referenceBases and alternateBases (LegacyVariation in
# EGA-archive/beacon2-ri-tools-v2, validators/ga4gh/vcf.py). Uppercase IUPAC codes only.
BEACON_ALLELE_PATTERN = re.compile(r"^([ACGTUNRYSWKMBDHV\-\.]*)$")


def find_invalid_alleles(vcf_path: str, max_examples: int = 10) -> Tuple[int, List[str]]:
    """Find records whose REF or any ALT allele does not match the Beacon pattern.

    Parameters
    ----------
    vcf_path : str
        Path to input VCF file.
    max_examples : int
        Maximum number of offending records to return.

    Returns
    -------
    Tuple[int, List[str]]
        Number of offending records, and up to `max_examples` of them as
        CHROM:POS:REF:ALT strings.

    """
    n_invalid = 0
    examples = []
    vcf = VCF(vcf_path)
    for record in vcf:
        alleles = [record.REF] + record.ALT
        if not all(BEACON_ALLELE_PATTERN.match(allele) for allele in alleles):
            n_invalid += 1
            if len(examples) < max_examples:
                examples.append(f"{record.CHROM}:{record.POS}:{record.REF}:{','.join(record.ALT)}")
    vcf.close()
    return n_invalid, examples


def main() -> None:
    """Parse CLI arguments and exit with an error if any allele is rejected."""
    parser = argparse.ArgumentParser(description="Check REF/ALT alleles against the Beacon allele pattern.")
    parser.add_argument("vcf", help="VCF to check, e.g. the sites.pass output")
    args = parser.parse_args()

    n_invalid, examples = find_invalid_alleles(args.vcf)
    if n_invalid:
        print(
            f"Error: {n_invalid} record(s) in {args.vcf} have REF/ALT alleles that Beacon rejects "
            f"(allowed: {BEACON_ALLELE_PATTERN.pattern}). First records:",
            file=sys.stderr,
        )
        for example in examples:
            print(f"  {example}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
