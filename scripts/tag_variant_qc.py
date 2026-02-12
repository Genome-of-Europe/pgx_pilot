import pysam
import argparse

def tag_variant_qc(input_vcf, output_vcf, thresholds):
    vcf_in = pysam.VariantFile(input_vcf)
    header = vcf_in.header
    
    if 'QC_STATUS' not in header.info:
        header.info.add('QC_STATUS', number='1', type='String', description='Variant QC status: PASS or semicolon-separated reasons for failure')

    for tag, (num, typ, desc) in {
        'HWE': ('A', 'Float', 'HWE p-value'),
        'MAF': ('1', 'Float', 'Minor Allele Frequency'),
        'F_MISSING': ('1', 'Float', 'Fraction of missing genotypes')
    }.items():
        if tag not in header.info:
            header.info.add(tag, number=num, type=typ, description=desc)

    vcf_out = pysam.VariantFile(output_vcf, 'w', header=header)

    for record in vcf_in:
        reasons = []
        
        if record.qual is not None and record.qual < thresholds['qual']:
            reasons.append(f"FAIL_QUAL")

        info = record.info
        
        if 'QD' in info and info['QD'] is not None and info['QD'] < thresholds['qd']:
            reasons.append("FAIL_QD")
        
        if 'MQ' in info and info['MQ'] is not None and info['MQ'] < thresholds['mq']:
            reasons.append("FAIL_MQ")
            
        if 'FS' in info and info['FS'] is not None and info['FS'] > thresholds['fs']:
            reasons.append("FAIL_FS")
            
        if 'ReadPosRankSum' in info and info['ReadPosRankSum'] is not None and info['ReadPosRankSum'] < thresholds['readpos']:
            reasons.append("FAIL_ReadPosRankSum")

        if 'DP' in info and info['DP'] is not None and info['DP'] < thresholds['min_dp']:
            reasons.append("FAIL_DP")

        if 'HWE' in info:
            hwe_val = info['HWE']
            if isinstance(hwe_val, (list, tuple)):
                hwe_val = hwe_val[0]
            if hwe_val is not None and hwe_val < thresholds['hwe']:
                reasons.append("FAIL_HWE")

        if 'MAF' in info:
            maf_val = info['MAF']
            if isinstance(maf_val, (list, tuple)):
                maf_val = maf_val[0]
            if maf_val is not None and maf_val < thresholds['maf']:
                reasons.append("FAIL_MAF")

        if 'F_MISSING' in info:
            f_missing = info['F_MISSING']
            if isinstance(f_missing, (list, tuple)):
                f_missing = f_missing[0]
            if f_missing is not None and f_missing > thresholds['max_missing']:
                reasons.append("FAIL_MISSING")

        status = "PASS"
        if reasons:
            status = ",".join(reasons)
        
        record.info['QC_STATUS'] = status
        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_vcf")
    parser.add_argument("output_vcf")
    
    parser.add_argument("--qual", type=float, default=30.0)
    parser.add_argument("--qd", type=float, default=2.0)
    parser.add_argument("--mq", type=float, default=40.0)
    parser.add_argument("--fs", type=float, default=60.0)
    parser.add_argument("--readpos", type=float, default=-8.0)
    
    parser.add_argument("--min_dp", type=int, default=10)
    parser.add_argument("--min_gq", type=int, default=20)
    parser.add_argument("--ab_ratio", type=float, default=0.2)
    parser.add_argument("--hwe", type=float, default=1e-6)
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--max_missing", type=float, default=0.1)
    
    args = parser.parse_args()
    
    thresholds = {
        'qual': args.qual,
        'qd': args.qd,
        'mq': args.mq,
        'fs': args.fs,
        'readpos': args.readpos,
        'min_dp': args.min_dp,
        'min_gq': args.min_gq,
        'ab_ratio': args.ab_ratio,
        'hwe': args.hwe,
        'maf': args.maf,
        'max_missing': args.max_missing
    }
    
    tag_variant_qc(args.input_vcf, args.output_vcf, thresholds)
