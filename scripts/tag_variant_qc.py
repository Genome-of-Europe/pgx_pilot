#!/usr/bin/env python3
"""Tag variants with quality control flags in a VCF/BCF file using cyvcf2."""

import argparse
import shlex
import sys
from typing import Any, Dict, Optional, Union

from cyvcf2 import VCF, Variant, Writer


def _extract_scalar(val: Any) -> Optional[float]:
    """Extract first numeric float value from a scalar, sequence, or string.

    Parameters
    ----------
    val : Any
        Value from VCF record INFO field.

    Returns
    -------
    Optional[float]
        Numeric float value, or None if missing or non-numeric.

    """
    if isinstance(val, (list, tuple)):
        val = val[0] if len(val) > 0 else None
    if val is None or val == "." or val == "":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def modify_vcf_header(vcf: VCF) -> None:
    """Ensure all required QC header fields and command line provenance are present.

    Parameters
    ----------
    vcf : VCF
        cyvcf2.VCF instance whose header will be updated before creating Writer.

    """
    command_args = " ".join(shlex.quote(arg) for arg in sys.argv)
    vcf.add_to_header(f"##commandline={command_args}")

    if "QC_STATUS" not in vcf:
        vcf.add_info_to_header(
            {
                "ID": "QC_STATUS",
                "Number": "1",
                "Type": "String",
                "Description": (
                    "Variant QC status: PASS or comma-separated reasons for failure"
                ),
            }
        )

    for tag, (num, typ, desc) in {
        "HWE": ("A", "Float", "HWE p-value"),
        "MAF": ("1", "Float", "Minor Allele Frequency"),
        "F_MISSING": ("1", "Float", "Fraction of missing genotypes"),
    }.items():
        if tag not in vcf:
            vcf.add_info_to_header(
                {
                    "ID": tag,
                    "Number": num,
                    "Type": typ,
                    "Description": desc,
                }
            )


def _get_threshold(
    thresholds: Union[Dict[str, Any], argparse.Namespace], key: str, default: Any = None
) -> Any:
    """Retrieve threshold cutoff from dict or Namespace object."""
    if isinstance(thresholds, dict):
        return thresholds.get(key, default)
    return getattr(thresholds, key, default)


def _tag_variant_records(
    vcf_in: VCF,
    vcf_out: Writer,
    thresholds: Union[Dict[str, Any], argparse.Namespace],
) -> None:
    """Iterate over records in vcf_in, compute QC_STATUS, and write to vcf_out."""
    qual_thresh = _get_threshold(thresholds, "qual", 30.0)
    qd_thresh = _get_threshold(thresholds, "qd", 2.0)
    mq_thresh = _get_threshold(thresholds, "mq", 40.0)
    fs_thresh = _get_threshold(thresholds, "fs", 60.0)
    readpos_thresh = _get_threshold(thresholds, "readpos", -8.0)
    min_dp_thresh = _get_threshold(thresholds, "min_dp", 10)
    hwe_thresh = _get_threshold(thresholds, "hwe", 1e-6)
    max_missing_thresh = _get_threshold(thresholds, "max_missing", 0.1)

    record: Variant
    for record in vcf_in:
        reasons = []

        if record.QUAL is not None and record.QUAL >= 0 and record.QUAL < qual_thresh:
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
        vcf_out.write_record(record)


def tag_variant_qc(
    input_vcf: Union[str, VCF],
    output_vcf: Union[str, Writer],
    thresholds: Union[Dict[str, Any], argparse.Namespace],
) -> None:
    """Tag variant records with QC_STATUS info field based on thresholds.

    Supports both file paths and active cyvcf2 VCF/Writer objects.

    Parameters
    ----------
    input_vcf : str or VCF
        Input VCF file path or open cyvcf2.VCF instance.
    output_vcf : str or Writer
        Output VCF file path or open cyvcf2.Writer instance.
    thresholds : dict or argparse.Namespace
        QC threshold cutoffs.

    """
    close_in = False
    close_out = False

    if isinstance(input_vcf, str):
        vcf_in = VCF(input_vcf)
        close_in = True
    else:
        vcf_in = input_vcf

    modify_vcf_header(vcf_in)

    if isinstance(output_vcf, str):
        vcf_out = Writer(output_vcf, vcf_in)
        close_out = True
    else:
        vcf_out = output_vcf

    try:
        _tag_variant_records(vcf_in, vcf_out, thresholds)
    finally:
        if close_in:
            vcf_in.close()
        if close_out:
            vcf_out.close()


def main() -> None:
    """Parse CLI arguments and run variant QC tagging."""
    parser = argparse.ArgumentParser(description="Tag variant QC in VCF.")
    parser.add_argument(
        "input_vcf", nargs="?", default="-", help="Input VCF path (default: stdin)"
    )
    parser.add_argument(
        "output_vcf",
        nargs="?",
        default=None,
        help="Output VCF path (optional positional)",
    )
    parser.add_argument(
        "-o", "--output", default=None, help="Output VCF path (default: stdout)"
    )
    parser.add_argument(
        "--in-threads",
        type=int,
        default=0,
        help="Decompression threads for input VCF [0]",
    )
    parser.add_argument(
        "-O",
        "--output-type",
        default="v",
        help="Output format: u/b: un/compressed BCF, v/z: un/compressed VCF [v]",
    )

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

    input_path = "/dev/stdin" if args.input_vcf in ("-", "", None) else args.input_vcf
    output_path = args.output or args.output_vcf or "/dev/stdout"

    ot = args.output_type
    # In cyvcf2, mode="wu" is treated as compressed VCF ("wz").
    # For uncompressed BCF streaming, mode="wb0" must be used.
    if ot.startswith("u"):
        ot = "b0"

    threads = args.in_threads or None

    with VCF(input_path, threads=threads) as in_vcf:
        modify_vcf_header(in_vcf)
        with Writer(output_path, in_vcf, mode="w" + ot) as out_vcf:
            tag_variant_qc(in_vcf, out_vcf, args)


if __name__ == "__main__":
    main()
