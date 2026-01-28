import pysam
import argparse

def tag_variant_qc(input_vcf, output_vcf, thresholds):
    vcf_in = pysam.VariantFile(input_vcf)
    header = vcf_in.header
    
    # Add QC_STATUS tag to header
    if 'QC_STATUS' not in header.info:
        header.info.add('QC_STATUS', number='1', type='String', description='Variant QC status: PASS or semicolon-separated reasons for failure')

    vcf_out = pysam.VariantFile(output_vcf, 'w', header=header)

    for record in vcf_in:
        reasons = []
        
        # QUAL
        if record.qual is not None and record.qual < thresholds['qual']:
            reasons.append(f"FAIL_QUAL")

        # INFO fields
        info = record.info
        
        # QD
        if 'QD' in info and info['QD'] < thresholds['qd']:
            reasons.append("FAIL_QD")
        
        # DP (Site depth)
        if 'DP' in info and info['DP'] < thresholds['dp']:
            reasons.append("FAIL_DP")
            
        # MQ
        if 'MQ' in info and info['MQ'] < thresholds['mq']:
            reasons.append("FAIL_MQ")
            
        # FS
        if 'FS' in info and info['FS'] > thresholds['fs']:
            reasons.append("FAIL_FS")
            
        # ReadPosRankSum
        if 'ReadPosRankSum' in info and info['ReadPosRankSum'] < thresholds['readpos']:
            reasons.append("FAIL_ReadPosRankSum")

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
    # Add threshold args with defaults
    parser.add_argument("--qual", type=float, default=30.0)
    parser.add_argument("--qd", type=float, default=2.0)
    parser.add_argument("--dp", type=int, default=10)
    parser.add_argument("--mq", type=float, default=40.0)
    parser.add_argument("--fs", type=float, default=60.0)
    parser.add_argument("--readpos", type=float, default=-8.0)
    
    args = parser.parse_args()
    
    thresh = {
        'qual': args.qual,
        'qd': args.qd,
        'dp': args.dp,
        'mq': args.mq,
        'fs': args.fs,
        'readpos': args.readpos
    }
    
    tag_variant_qc(args.input_vcf, args.output_vcf, thresh)
