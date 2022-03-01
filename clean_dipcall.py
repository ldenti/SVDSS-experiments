import sys
from pysam import VariantFile

def main():
    vcf_path = sys.argv[1]

    vcf = VariantFile(vcf_path)
    print(vcf.header, end='')
    for rec in vcf.fetch():
        if len(rec.alts) == 1:
            print(rec, end='')
        else:
            alts = rec.alts
            h1, h2 = rec.samples[0]["GT"]
            if h1 != 0:
                rec.alts = (alts[0],)
                rec.samples[0]["GT"] = (1,0)
                print(rec, end='')
            if h2 != 0:
                rec.alts = (alts[1],)
                rec.samples[0]["GT"] = (0,1)
                print(rec, end='')

if __name__ == "__main__":
    main()