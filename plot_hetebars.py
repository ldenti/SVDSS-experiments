import sys, os, glob
from collections import Counter
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_vcf(vcf_path, tag=""):
    vcf = VariantFile(vcf_path)
    Vs = {}
    for rec in vcf.fetch():
        if list(rec.filter) != ["PASS"] and list(rec.filter) != []:
            continue
        if "SVLEN" in rec.info:
            l = rec.info["SVLEN"]
        else:
            l = len(rec.alts[0]) - len(rec.ref)
        if type(l) == type(()):
            l = l[0]
        if abs(l) < 50:
            continue

        a1, a2 = rec.samples[0]["GT"]
        a1 = 0 if a1 == None else a1
        a2 = 0 if a2 == None else a2
        gt = f"{min(a1,a2)}/{max(a1,a2)}"
        key = (rec.contig, rec.pos)
        if key in Vs:
            Vs[key] = "1/2"
        else:
            Vs[key] = gt
    return Vs


def main():
    dipcall_path = sys.argv[1]
    giab_path = sys.argv[2]
    out_dir = sys.argv[3]

    df_dip = pd.DataFrame(columns=["VCF", "Type", "Count"])
    df_giab = pd.DataFrame(columns=["VCF", "Type", "Count"])

    dip = parse_vcf(dipcall_path)
    dip = dict(Counter(dip.values()))
    for k, v in dip.items():
        df_dip = df_dip.append(
            {"VCF": "Truth", "Type": k, "Count": v},
            ignore_index=True,
        )

    giab = parse_vcf(giab_path)
    giab = dict(Counter(giab.values()))
    giab["1/2"] = 0
    for k, v in giab.items():
        df_giab = df_giab.append(
            {"VCF": "Truth", "Type": k, "Count": v},
            ignore_index=True,
        )

    for vcf_path in glob.glob(
        os.path.join(out_dir, "truvari", "*", "full", "tp-base.vcf")
    ):
        tool = vcf_path.split("/")[-3]
        if tool == "pp":
            tool = "SVDSS"
        elif tool == "cutesv":
            tool = "cuteSV"
        elif tool == "svim":
            tool == "SVIM"
        Vs = parse_vcf(vcf_path)
        Vs = dict(Counter(Vs.values()))
        for k, v in Vs.items():
            df_dip = df_dip.append(
                {"VCF": tool, "Type": k, "Count": v},
                ignore_index=True,
            )

    for vcf_path in glob.glob(
        os.path.join(out_dir, "truvari-giab", "*", "full", "tp-base.vcf")
    ):
        tool = vcf_path.split("/")[-3]
        if tool == "pp":
            tool = "SVDSS"
        elif tool == "cutesv":
            tool = "cuteSV"
        elif tool == "svim":
            tool == "SVIM"
        Vs = parse_vcf(vcf_path)
        Vs = dict(Counter(Vs.values()))
        Vs["1/2"] = 0
        for k, v in Vs.items():
            df_giab = df_giab.append(
                {"VCF": tool, "Type": k, "Count": v},
                ignore_index=True,
            )

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(7, 7))

    print(df_giab)
    print(df_dip)

    sns.barplot(
        x="VCF",
        y="Count",
        hue="Type",
        hue_order=["0/1", "1/1", "1/2"],
        orient="v",
        data=df_dip,
        ax=ax1,
    )
    sns.barplot(
        x="VCF",
        y="Count",
        hue="Type",
        hue_order=["0/1", "1/1", "1/2"],
        orient="v",
        data=df_giab,
        ax=ax2,
    )

    ax1.set_title("Truth = dipcall")
    ax2.set_title("Truth = GIAB")
    ax1.legend_.remove()
    ax2.set_ylabel("")
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=60)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=60)

    plt.tight_layout()
    plt.savefig("hetbar.png")


if __name__ == "__main__":
    main()
