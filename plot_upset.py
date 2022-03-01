import sys
import os, glob
from itertools import combinations
from pysam import VariantFile
from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot as plt

# import seaborn as sns

plt.rcParams.update({"font.size": 17})
# sns.set_style("whitegrid")

# colors = [sns.color_palette("bright")[1]] + [
#     sns.color_palette("dark")[i] for i in [0, 2, 6, 7, 9]
# ]


class Variant:
    def __init__(self, chrom, pos, ref, alt, l):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.l = l
        self.end = self.pos if self.l > 0 else self.pos + (-self.l)
        self.type = "DEL" if self.l < 0 else "INS"
        self.idx = "."

    def __hash__(self):
        return hash((self.chrom, self.pos, self.ref, self.alt, self.l))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return NotImplemented
        return (
            self.chrom == other.chrom
            and self.pos == other.pos
            and self.ref == other.ref
            and self.alt == other.alt
            and self.l == other.l
        )

    def __repr__(self):
        return (
            self.chrom
            + "\t"
            + f"{self.pos}"
            + "\t"
            + self.idx
            + "\t"
            + self.ref
            + "\t"
            + self.alt
            + "\t"
            + "."
            + "\t"
            + "PASS"
            + "\t"
            + f"VARTYPE=SV;SVTYPE={self.type};SVLEN={self.l};END={self.end}"
            + "\t"
            + "GT"
            + "\t"
            + "0/1"
        )


def parse_vcf(vcf_path, print_header=False):
    vcf = VariantFile(vcf_path)
    if print_header:
        print(vcf.header, end="")
    Vs = {}
    for rec in vcf.fetch():
        if rec.chrom not in Vs:
            Vs[rec.chrom] = set()
        if "SVLEN" in rec.info:
            l = rec.info["SVLEN"]
        else:
            l = len(rec.alts[0]) - len(rec.ref)
        if type(l) == type(()):
            l = l[0]

        # if abs(l) < 50:
        #     continue

        Vs[rec.chrom].add(Variant(rec.chrom, rec.pos, rec.ref, rec.alts[0], l))
    return Vs


def plot_upset():
    truv_dir = sys.argv[1]
    mode = sys.argv[2]  # e.g., "tier2conf"

    truth = {}
    data = {}
    first = True
    VCFs = list(
        glob.glob(os.path.join(truv_dir, "*", mode, "tp-base.vcf"))
    )  # + [os.path.join(truv_dir, "..", "..", "dip-vs-giab", "tier2conf", "tp-base.vcf")]
    for vcf_path in VCFs:
        tool = vcf_path.split("/")[-3]
        if tool == "pp":
            tool = "SVDSS"
        elif tool == "svim":
            tool = "SVIM"
        elif tool == "cutesv":
            tool = "cuteSV"
        elif tool == "sniffles2":
            continue
        data[tool] = parse_vcf(vcf_path, tool == "SVDSS")
        if first:
            truth = parse_vcf(vcf_path)
            FNs = parse_vcf("/".join(vcf_path.split("/")[:-1] + ["fn.vcf"]))
            for chrom, Vs in FNs.items():
                if chrom not in truth:
                    truth[chrom] = []
                truth[chrom] |= Vs
            first = False

    for tool in data:
        N = 0
        for chrom in data[tool]:
            N += len(data[tool][chrom])
        print(tool, N, file=sys.stderr)

    N = 0
    for chrom in truth:
        N += len(truth[chrom])
    print("truth", N, file=sys.stderr)

    df = {}
    for n in range(len(data) + 1):
        for comb in combinations(data, n):
            df["-".join(sorted(comb))] = 0

    for chrom, Vs in truth.items():
        for v in Vs:
            Ts = []
            for tool in data:
                if chrom not in data[tool]:
                    continue
                if v in data[tool][chrom]:
                    Ts.append(tool)
            if Ts == ["SVDSS"]:
                print(v)
            df["-".join(sorted(Ts))] += 1

    # print(df)
    keys = [k.split("-") if len(k) > 0 else [] for k in sorted(list(df.keys()))]
    df_membership = from_memberships(
        keys, data=[df[k] for k in sorted(list(df.keys()))]
    )
    plot(df_membership)

    plt.savefig("upset.png", dpi=300)
    plt.savefig("upset.svg", dpi=300)


# def plot_supervenn():
#     truv_dir = sys.argv[1]
#     mode = sys.argv[2]  # e.g., "tier2conf"

#     truth = {}
#     data = {}
#     first = True
#     VCFs = list(glob.glob(os.path.join(truv_dir, "*", mode, "tp-base.vcf")))
#     for vcf_path in VCFs:
#         tool = vcf_path.split("/")[-3]
#         if tool == "pp":
#             tool = "SVDSS"
#         elif tool == "svim":
#             tool = "SVIM"
#         elif tool == "cutesv":
#             tool = "cuteSV"
#         elif tool == "sniffles":
#             continue
#         data[tool] = parse_vcf(vcf_path, tool == "SVDSS")
#         if first:
#             truth = parse_vcf(vcf_path)
#             FNs = parse_vcf("/".join(vcf_path.split("/")[:-1] + ["fn.vcf"]))
#             for chrom, Vs in FNs.items():
#                 if chrom not in truth:
#                     truth[chrom] = []
#                 truth[chrom] |= Vs
#             first = False

#     df = {}
#     for tool in data:
#         N = 0
#         df[tool] = set()
#         for chrom in data[tool]:
#             df[tool] |= set(data[tool][chrom])
#             N += len(data[tool][chrom])
#         print(tool, N, file=sys.stderr)
#     df["dipcall (truthset)"] = set()
#     for chrom in truth:
#         df["dipcall (truthset)"] |= set(truth[chrom])
#         N += len(truth[chrom])
#     print("dipcall (truthset)", N, file=sys.stderr)

#     labels = [
#         "dipcall (truthset)",
#         "SVDSS",
#         "cuteSV",
#         "pbsv",
#         "sniffles2",
#         "SVIM",
#         "debreak",
#     ]
#     supervenn_sets = [df[k] for k in labels]
#     plt.figure(figsize=(21, 7), dpi=300)

#     supervenn(
#         supervenn_sets,
#         labels,
#         chunks_ordering="occurrence",
#         widths_minmax_ratio=0.5,
#         color_cycle=colors,
#     )
#     plt.xlabel("Calls on HG007 (Ext. Tier 2)")
#     plt.ylabel("Tools")
#     plt.tight_layout()
#     plt.savefig("upset-supervenn.png", dpi=300)


if __name__ == "__main__":
    plot_upset()
