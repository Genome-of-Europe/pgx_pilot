"""VCF and BCF format helpers for pipeline output configuration."""


def vcf_format_from_output_type(ot: str) -> str:
    """Determine the file extension from bcftools output type (b|u|z|v[0-9]).

    Parameters
    ----------
    ot : str
        bcftools output type code, e.g. 'b1', 'u', 'z4', 'v'.

    Returns
    -------
    str
        File extension corresponding to the format: 'bcf', 'vcf', or 'vcf.gz'.

    Raises
    ------
    ValueError
        If the output type is unsupported.

    """
    if ot.startswith("v"):
        return "vcf"
    if ot.startswith("z"):
        return "vcf.gz"
    if ot.startswith("b"):
        return "bcf"
    if ot.startswith("u"):
        if len(ot) > 2:
            raise ValueError("Compression not supported for uncompressed format")
        return "bcf"
    raise ValueError(f"Unsupported output type {ot=}")
