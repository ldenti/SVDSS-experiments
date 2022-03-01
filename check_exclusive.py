import sys
from pysam import VariantFile

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_vcf(vcf_path):
    vcf = VariantFile(vcf_path)
    Vs = {}
    for rec in vcf.fetch():
        if rec.chrom not in Vs:
            Vs[rec.chrom] = []
        Vs[rec.chrom].append(rec.pos)
    return Vs


def main():
    exclusive_vcf_path = sys.argv[1]
    dipcall_vcf_path = sys.argv[2]

    exclusive_vars = parse_vcf(exclusive_vcf_path)
    dipcall_vars = parse_vcf(dipcall_vcf_path)

    for chrom, exVs in exclusive_vars.items():
        for exV in exVs:
            min_d = float("inf")
            count = 0
            for dipV in dipcall_vars[chrom]:
                if abs(dipV - exV) < min_d:
                    min_d = abs(dipV - exV)
                    count = 1
                elif abs(dipV - exV) == min_d:
                    count += 1
            if min_d == 0 and count == 1:
                min_d = float("inf")
                count = 0
                for dipV in dipcall_vars[chrom]:
                    if abs(dipV - exV) < min_d and abs(dipV - exV) != 0:
                        min_d = abs(dipV - exV)
                        count = 1
                    elif abs(dipV - exV) == min_d:
                        count += 1
            print(chrom, exV, min_d, count)


def bars():
    txt_path = sys.argv[1]

    Ds = []
    for line in open(txt_path):
        d = int(line.strip("\n").split(" ")[2])
        Ds.append(d)

    max_d = 250
    bars = {}
    for d in Ds:
        if d == 0:
            key = "0bp"
        elif d > max_d:
            key = f">{max_d}bp"
        else:
            min_split, max_split = int(d / 50) * 50 + 1, int(d / 50) * 50 + 50
            key = f"{min_split}-{max_split}bp"
        if key not in bars:
            bars[key] = 0
        bars[key] += 1

    keys = [
        "0bp",
        "1-50bp",
        "51-100bp",
        "101-150bp",
        "151-200bp",
        "201-250bp",
        ">250bp",
    ]
    values = [bars[k] for k in keys]
    bars_dict = {"Distance": keys, "Count": values}
    df = pd.DataFrame.from_dict(bars_dict)

    print(df)
    sns.barplot(x="Distance", y="Count", data=df)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("dist-bars.png")


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "get":
        main()
    elif mode == "plot":
        bars()
