"""Generate a map for renaming chromosomes, to be fed to bcftools."""
import argparse


def check_for_chr(vcf_path):
    from cyvcf2 import VCF
    with VCF(vcf_path) as vcf:
        return any(seq.startswith("chr") for seq in vcf.seqnames)


def main() -> None:
    parser=argparse.ArgumentParser(description="Generate chromosome rename file from Ensembl to USCS (e.g 1 to chr1)")
    parser.add_argument("input",help="Input VCF path")
    parser.add_argument("output",help="Output map file")
    args=parser.parse_args()

    print("write", args.output)
    with open(args.output,"w") as f:
        has_chr = check_for_chr(args.input)

        if not has_chr:
            for i in range(1,23):
                f.write(f"{i}\tchr{i}\n")
            f.write("X\tchrX\n")
            f.write("Y\tchrY\n")
            f.write("MT\tchrM\n")
        else:
            # No-op (empty map file means no renaming)
            pass

if __name__=="__main__":
    main()
