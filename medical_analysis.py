import sys, os, glob
from pysam import VariantFile
import pandas as pd
from matplotlib.pyplot import figure
import matplotlib.pyplot as plt
import seaborn as sns
from venn import venn
from supervenn import supervenn

# plt.rcParams.update({"font.size": 27})
sns.set_style("whitegrid")

colors = [sns.color_palette("bright")[1]] + [sns.color_palette("dark")[i] for i in [0,2,6,7,9]]


def parse_vcf(vcf_path):
    vcf = VariantFile(vcf_path)
    Vs = set()
    for rec in vcf.fetch():
        l = len(rec.alts[0]) - len(rec.ref)
        idx = f"{rec.chrom}-{rec.pos}-{rec.ref}-{rec.alts[0]}"
        Vs.add(idx)
    return Vs


def parse_bed(bed_path):
    genes = {}
    for line in open(bed_path):
        chrom, s, e, gidx = line.strip("\n").split("\t")
        genes[(chrom, s, e)] = gidx
    return genes


def build_gene_tp_map(tps, pos2gene):
    tp_map = {gidx: 0 for gidx in pos2gene.values()}
    for v in tps:
        chrom, p = v.split("-")[0], v.split("-")[1]
        for (gchrom, s, e), gidx in pos2gene.items():
            if chrom != gchrom:
                continue
            if s <= p and p <= e:
                tp_map[gidx] += 1
    return tp_map


def plot_venn():
    truv_dir = sys.argv[1]
    data = {}
    for vcf_path in list(glob.glob(os.path.join(truv_dir, "*", "tp-base.vcf"))):
        tool = vcf_path.split("/")[-2]
        if tool == "pp":
            tool = "SVDSS"
        elif tool == "svim":
            tool = "SVIM"
        elif tool == "cutesv":
            tool = "cuteSV"
        elif tool == "sniffles2":
            continue
        tp = parse_vcf(vcf_path)
        data[tool] = tp

    venn_labels = ["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM"]  # , "debreak"]
    venn_dict = {k: data[k] for k in venn_labels}

    supervenn_labels = venn_labels + ["debreak"]
    supervenn_sets = [data[k] for k in supervenn_labels]

    plt.figure(figsize=(7, 7), dpi=300)

    venn(venn_dict, cmap=colors, fontsize=17)
    plt.tight_layout()
    plt.savefig("medical-venn.png", dpi=300)
    plt.savefig("medical-venn.svg", dpi=300)

    plt.clf()
    supervenn(
        supervenn_sets,
        supervenn_labels,
        chunks_ordering="occurrence",
        widths_minmax_ratio=0.05,
        color_cycle=colors,
        fontsize=17,
    )
    plt.xlabel("True positive calls from CMRG callset")
    plt.ylabel("Tools")
    plt.tight_layout()
    plt.savefig("medical-supervenn.png", dpi=300)
    plt.savefig("medical-supervenn.svg", dpi=300)


def main():
    truth_vcf_path = sys.argv[1]
    genes_bed_path = sys.argv[2]
    truv_dir = sys.argv[3]

    pos2gene = parse_bed(genes_bed_path)

    # Getting genes effectively exhibithing at least one SV
    truth = parse_vcf(truth_vcf_path)
    truth_map = build_gene_tp_map(truth, pos2gene)
    genes = set([gidx for gidx in truth_map if truth_map[gidx] != 0])
    print(f"Considering {len(genes)} instead of {len(pos2gene.values())}")

    # All genes (even if not covered by a SV)
    # genes = set(pos2gene.values())

    matrix = []
    venn_dict = {}
    for vcf_path in list(glob.glob(os.path.join(truv_dir, "*", "tp-base.vcf"))):
        tool = vcf_path.split("/")[-2]
        if tool == "pp":
            tool = "SVDSS"
        elif tool == "svim":
            tool = "SVIM"
        elif tool == "cutesv":
            tool = "cuteSV"
        tp = parse_vcf(vcf_path)
        venn_dict[tool] = tp
        tp_map = build_gene_tp_map(tp, pos2gene)
        for gidx in genes:
            matrix.append([gidx, tool, tp_map[gidx]])

    df = pd.DataFrame(columns=["Gene", "Tool", "SVs"], data=matrix)
    df_heatmap = df.pivot("Gene", "Tool", "SVs")

    # figure(figsize=(8, 24), dpi=150)
    # sns.heatmap(df_heatmap)
    # plt.tight_layout()
    # plt.savefig("heatmap.png")

    venn_dict = {
        k: venn_dict[k] for k in ["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM"]
    }
    # plt.clf()
    plt.figure(figsize=(7, 7), dpi=300)
    v = venn(venn_dict, cmap=colors, fontsize=17)

    # plt.title("(d) TPs distribution against CMRG")
    plt.tight_layout()
    plt.savefig("medical-venn.png", dpi=300)
    plt.savefig("medical-venn.svg", dpi=300)

    print("Found by SVDSS:")
    for v in venn_dict["SVDSS"]:
        gidx = ""
        for (chrom, s, e), gidx in pos2gene.items():
            if s <= v.split("-")[1] and v.split("-")[1] <= e:
                break
        print("", v.split("-")[0], v.split("-")[1], gidx)
    print("")
    print("Found by SVDSS only:")
    for v in (
        venn_dict["SVDSS"]
        - venn_dict["pbsv"]
        - venn_dict["cuteSV"]
        - venn_dict["SVIM"]
        - venn_dict["sniffles"]
    ):
        gidx = ""
        for (chrom, s, e), gidx in pos2gene.items():
            if s <= v.split("-")[1] and v.split("-")[1] <= e:
                break
        print("", v.split("-")[0], v.split("-")[1], gidx)
    print("")
    print("Not found by SVDSS:")
    for v in venn_dict["pbsv"] - venn_dict["SVDSS"]:
        gidx = ""
        for (chrom, s, e), gidx in pos2gene.items():
            if s <= v.split("-")[1] and v.split("-")[1] <= e:
                break
        print("", v.split("-")[0], v.split("-")[1], gidx)


if __name__ == "__main__":
    plot_venn()
