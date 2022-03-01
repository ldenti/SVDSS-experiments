import sys
from pysam import VariantFile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({"font.size": 17})
sns.set_style("whitegrid")

colors = [sns.color_palette("bright")[1]] + [
    sns.color_palette("dark")[i] for i in [0, 2, 6, 7, 9]
]


def parse_vcf(vcf_path, tag=""):
    vcf = VariantFile(vcf_path)
    Vs = []
    for rec in vcf.fetch():
        if list(rec.filter) != ["PASS"] and list(rec.filter) != []:
            continue
        if "SVLEN" in rec.info:
            l = rec.info["SVLEN"]
        else:
            l = len(rec.alts[0]) - len(rec.ref)
        if type(l) == type(()):  # fix for pbsv
            l = l[0]
        if rec.info["SVTYPE"] == "DEL" and l > 0:  # fix for debreak
            l = -l
        if abs(l) < 30:
            continue
        if l > 0:
            Vs.append([tag, l, "Ins"])
        else:
            Vs.append([tag, -l, "Del"])
    return Vs


def main():
    df = pd.DataFrame(columns=["VCF", "Length", "Type"], data=[])
    for name, vcf_path in zip(sys.argv[1::2], sys.argv[2::2]):
        Vs = parse_vcf(vcf_path, name)
        print(name, len(Vs))
        print(name, "del", len([v for v in Vs if v[2] == "Del"]))
        print(name, "ins", len([v for v in Vs if v[2] == "Ins"]))
        new_df = pd.DataFrame(columns=["VCF", "Length", "Type"], data=Vs)
        df = pd.concat([df, new_df], ignore_index=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(7, 7))

    sns.kdeplot(
        data=df.loc[df["Type"] == "Del"],
        x="Length",
        log_scale=True,
        hue="VCF",
        palette=colors,
        hue_order=["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM", "debreak"],
        legend=False,
        ax=ax1,
        linewidth=1.3,
    )
    sns.kdeplot(
        data=df.loc[df["Type"] == "Ins"],
        x="Length",
        log_scale=True,
        hue="VCF",
        palette=colors,
        hue_order=["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM", "debreak"],
        ax=ax2,
        linewidth=1.3,
    )

    ax1.set_title("Deletions")
    ax2.set_title("Insertions")
    ax1.set_xlabel("")
    ax2.set_xlabel("")
    ax1.set_ylabel("Density")

    plt.xscale("symlog")
    ax1.invert_xaxis()
    fig.supxlabel("Length")

    plt.subplots_adjust(wspace=0, hspace=0)

    ax2.set_xticks(ax2.get_xticks()[2::2])
    ax1.set_xticks(ax2.get_xticks())

    # plt.suptitle("(a) SVs size distribution")
    plt.savefig("var-hist.png", dpi=300)
    plt.savefig("var-hist.svg", dpi=300)


if __name__ == "__main__":
    main()
