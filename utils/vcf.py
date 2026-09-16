
def vcf_format_from_output_type(ot: str) -> str:
    """Determine the format of the VCF file based on output type b|u|z|v[0-9]

    Output type represents: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v).
    The compression level of the compressed formats (b and z) can be set by appending a number between 0-9.
    """
    if ot.startswith("v"): return "vcf"
    if ot.startswith("z"): return "vcf.gz"
    if ot.startswith("b"): return "bcf"
    if ot.startswith("u"):
        if len(ot) > 2:
            raise ValueError(f"Compression not supported  for uncompressed VCF")
        return "bcf"
    raise ValueError(f"Unsupported output type {ot=}")


# TODO: This is WIP. Actually bcftools annotate is dificult as we have to rename and assign (duplicate).
def vcf_tag_suffix(tags: str,suffix: str,old_suffix: str= "",where="INFO") -> str:
    assignments = []
    for tag in tags.split(","):
        assignments.append( f"{where}/{tag}{suffix}:={where}/{tag}{old_suffix}")
    return ",".join(assignments)

RAW_ASSIGNMENTS = ",".join([
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_raw_BE"),
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_raw_M_BE",old_suffix="_M"),
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_raw_F_BE",old_suffix="_F"),
])

FINAL_ASSIGNMENTS = ",".join([
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_BE"),
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_M_BE",old_suffix="_M"),
    vcf_tag_suffix("AC,AN,AF,AC_Het,AC_Hom,AC_Hemi,MAF,HWE,NS,F_MISSING,ExcHet","_F_BE",old_suffix="_F"),
])


