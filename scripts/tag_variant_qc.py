"""Tag variants with quality control flags in a VCF file using cyvcf2."""

import argparse
from typing import Any, Dict
from cyvcf2 import VCF, Writer


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


def tag_variant_qc(input_vcf: str, output_vcf: str, thresholds: Dict[str, Any]) -> None:
    """Tag variant records with QC_STATUS info field based on thresholds.

    Parameters
    ----------
    input_vcf : str
        Path to input VCF file.
    output_vcf : str
        Path to output VCF file (.vcf.gz).
    thresholds : dict
        QC threshold dictionary with cutoffs.
    """
    vcf_in = VCF(input_vcf)

    if "QC_STATUS" not in vcf_in:
        vcf_in.add_info_to_header({
            "ID": "QC_STATUS",
            "Number": "1",
            "Type": "String",
            "Description": "Variant QC status: PASS or comma-separated reasons for failure",
        })

    # Ensure required metrics are present in header if missing
    for tag, (num, typ, desc) in {
        "HWE": ("A", "Float", "HWE p-value"),
        "MAF": ("1", "Float", "Minor Allele Frequency"),
        "F_MISSING": ("1", "Float", "Fraction of missing genotypes"),
    }.items():
        if tag not in vcf_in:
            vcf_in.add_info_to_header({
                "ID": tag,
                "Number": num,
                "Type": typ,
                "Description": desc,
            })

    vcf_out = Writer(output_vcf, vcf_in)

    qual_thresh = thresholds["qual"]
    qd_thresh = thresholds["qd"]
    mq_thresh = thresholds["mq"]
    fs_thresh = thresholds["fs"]
    readpos_thresh = thresholds["readpos"]
    min_dp_thresh = thresholds["min_dp"]
    hwe_thresh = thresholds["hwe"]
    maf_thresh = thresholds["maf"]
    max_missing_thresh = thresholds["max_missing"]

    for record in vcf_in:
        reasons = []

        if record.QUAL is not None and record.QUAL < qual_thresh:
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

        maf = _extract_scalar(info.get("MAF"))
        if maf is not None and maf < maf_thresh:
            reasons.append("FAIL_MAF")

        f_missing = _extract_scalar(info.get("F_MISSING"))
        if f_missing is not None and f_missing > max_missing_thresh:
            reasons.append("FAIL_MISSING")

        status = ",".join(reasons) if reasons else "PASS"
        record.INFO["QC_STATUS"] = status
        vcf_out.write_record(record)

    vcf_in.close()
    vcf_out.close()


def main() -> None:
    """Parse CLI arguments and run variant QC tagging."""
    parser = argparse.ArgumentParser(description="Tag variant QC in VCF.")
    parser.add_argument("input_vcf", help="Input VCF path")
    parser.add_argument("output_vcf", help="Output VCF path")

    parser.add_argument("--qual", type=float, default=30.0)
    parser.add_argument("--qd", type=float, default=2.0)
    parser.add_argument("--mq", type=float, default=40.0)
    parser.add_argument("--fs", type=float, default=60.0)
    parser.add_argument("--readpos", type=float, default=-8.0)

    parser.add_argument("--min_dp", type=int, default=10)
    parser.add_argument("--min_gq", type=int, default=20)
    parser.add_argument("--ab_ratio", type=float, default=0.2)
    parser.add_argument("--hwe", type=float, default=1e-6)
    parser.add_argument("--maf", type=float, default=0.0)
    parser.add_argument("--max_missing", type=float, default=0.1)

    args = parser.parse_args()

    thresholds = {
        "qual": args.qual,
        "qd": args.qd,
        "mq": args.mq,
        "fs": args.fs,
        "readpos": args.readpos,
        "min_dp": args.min_dp,
        "min_gq": args.min_gq,
        "ab_ratio": args.ab_ratio,
        "hwe": args.hwe,
        "maf": args.maf,
        "max_missing": args.max_missing,
    }

    tag_variant_qc(args.input_vcf, args.output_vcf, thresholds)


if __name__ == "__main__":
    main()
