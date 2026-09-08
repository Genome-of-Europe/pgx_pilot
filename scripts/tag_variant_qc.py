#!/usr/bin/env python

"""Tag variants with quality control flags in a VCF file using cyvcf2."""

import argparse
import shlex
import sys
from typing import Any

from cyvcf2 import VCF, Writer, Variant

def _extract_scalar(val: Any) -> Any:
    """Extract first element if value is a tuple or list, otherwise return as-is.

    Parameters
    ----------
    val : Any
        Value from VCF record INFO field.

    Returns
    -------
    Any
        Scalar float/int or None.
    """
    if isinstance(val, (list, tuple)):
        return val[0] if len(val) > 0 else None
    return val


def tag_variant_qc(in_vcf: VCF, out_vcf: Writer, thresholds: argparse.Namespace) -> None:
    """Tag variant records with QC_STATUS info field based on thresholds.

    Parameters
    ----------
    in_vcf : VCF
        cyfcf2.VCF file
    out_vcf : Writer
        cyfcf2.Writer file
    thresholds : argparse.Namespace
        QC threshold object with cutoffs.
    """

    qual_thresh = thresholds.qual
    qd_thresh = thresholds.qd
    mq_thresh = thresholds.mq
    fs_thresh = thresholds.fs
    readpos_thresh = thresholds.readpos
    min_dp_thresh = thresholds.min_dp
    hwe_thresh = thresholds.hwe
    max_missing_thresh = thresholds.max_missing

    record: Variant
    for record in in_vcf:
        reasons = []

        if record.QUAL is not None and 0<=record.QUAL<qual_thresh:
            reasons.append("FAIL_QUAL")

        info = record.INFO

        qd = _extract_scalar(info.get("QD"))
        if qd is not None and qd < qd_thresh:
            reasons.append("FAIL_QD")

        mq = _extract_scalar(info.get("MQ"))
        if mq is not None and mq < mq_thresh:
            reasons.append("FAIL_MQ")

        fs = _extract_scalar(info.get("FS"))
        if fs is not None and fs > fs_thresh:
            reasons.append("FAIL_FS")

        readpos = _extract_scalar(info.get("ReadPosRankSum"))
        if readpos is not None and readpos < readpos_thresh:
            reasons.append("FAIL_ReadPosRankSum")

        dp = _extract_scalar(info.get("DP"))
        if dp is not None and dp < min_dp_thresh:
            reasons.append("FAIL_DP")

        hwe = _extract_scalar(info.get("HWE"))
        if hwe is not None and hwe < hwe_thresh:
            reasons.append("FAIL_HWE")

        f_missing = _extract_scalar(info.get("F_MISSING"))
        if f_missing is not None and f_missing > max_missing_thresh:
            reasons.append("FAIL_MISSING")

        status = ",".join(reasons) if reasons else "PASS"
        record.INFO["QC_STATUS"] = status
        out_vcf.write_record(record)


def modify_vcf_header(vcf: VCF) -> None:
    """Ensure all the header fields are present."""

    # Reconstruct the command (simple version), safely handles string and quotes via shelx.
    command_args=" ".join(shlex.quote(arg) for arg in sys.argv)
    vcf.add_to_header(f"##commandline={command_args}")

    # TODO: we should use official FILTER
    if "QC_STATUS" not in vcf:
        vcf.add_info_to_header({
            "ID":"QC_STATUS",
            "Number":"1",
            "Type":"String",
            "Description":"Variant QC status: PASS or comma-separated reasons for failure",
        })

    # Ensure required metrics are present in header if missing
    for tag, (num, typ, desc) in {
        "HWE": ("A", "Float", "HWE p-value"),
        "MAF": ("1", "Float", "Minor Allele Frequency"),
        "F_MISSING": ("1", "Float", "Fraction of missing genotypes"),
    }.items():
        if tag not in vcf:
            info_row = { "ID": tag, "Number": num, "Type": typ, "Description": desc, }
            vcf.add_info_to_header(info_row)

def main() -> None:
    """Parse CLI arguments and run variant QC tagging."""
    parser = argparse.ArgumentParser(description="Tag variant QC in VCF.")
    parser.add_argument("input_vcf", nargs="?", default="-", help="Input VCF path")
    parser.add_argument("-o", "--output", default="/dev/stdout", help="Output VCF path")

    parser.add_argument("--in-threads", help="Same as bcftools --threads but only for de-compression [0]", type=int, default=0)
    parser.add_argument("-O", "--output-type", help="u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]", default="v")

    parser.add_argument("--qual", type=float, default=30.0)
    parser.add_argument("--qd", type=float, default=2.0)
    parser.add_argument("--mq", type=float, default=40.0)
    parser.add_argument("--fs", type=float, default=60.0)
    parser.add_argument("--readpos", type=float, default=-8.0)

    parser.add_argument("--min_dp", type=int, default=10)
    parser.add_argument("--min_gq", type=int, default=20)
    parser.add_argument("--ab_ratio", type=float, default=0.2)
    parser.add_argument("--hwe", type=float, default=1e-6)
    parser.add_argument("--max_missing", type=float, default=0.1)
    args = parser.parse_args()

    input_path = "/dev/stdin" if args.input_vcf == "-" else args.input_vcf

    ot = args.output_type
    if ot.startswith("u"): # cyvcf2 thinks this is z, likely a bug
        ot = "b0"

    threads = args.in_threads or None # This is a quirk of cyvcf2 not handling threads=0

    with VCF(input_path, threads=threads) as in_vcf:
        # header can't be changed on out_vcf so modify in_vcf as it serves as template
        modify_vcf_header(in_vcf)
        with Writer(args.output, in_vcf, mode="w"+ot) as out_vcf:
            tag_variant_qc(in_vcf, out_vcf, args)
            pass


if __name__ == "__main__":
    main()
