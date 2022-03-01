import sys
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

plt.rcParams.update({"font.size": 17})
sns.set_style("whitegrid")


colors = [sns.color_palette("bright")[1]] + [
    sns.color_palette("dark")[i] for i in [0, 2, 6, 7, 9]
]


def parse_res(fpath, label=""):
    data = []
    for line in open(fpath):
        if line.startswith("fullconf"):
            line = line.strip("\n").split(",")
            tool, P, R, F1 = line[1], float(line[5]), float(line[6]), float(line[7])
            if tool == "pp":
                tool = "SVDSS"
            elif tool == "cutesv":
                tool = "cuteSV"
            elif tool == "svim":
                tool = "SVIM"
            elif tool == "sniffles2":
                continue
            data.append([label, tool, P, R, F1])
    return data


def main():
    inres_5x = sys.argv[1]
    inres_10x = sys.argv[2]
    inres_15x = sys.argv[3]

    res_5x = parse_res(inres_5x, "5x")
    res_10x = parse_res(inres_10x, "10x")
    res_15x = parse_res(inres_15x, "15x")
    res = res_5x + res_10x + res_15x
    df = pd.DataFrame(res, columns=["Coverage", "Tool", "Precision", "Recall", "F1"])
    print(df)

    plt.figure(figsize=(7, 7))

    sns.lineplot(
        data=df,
        x="Precision",
        y="Recall",
        hue="Tool",
        palette=colors,
        hue_order=["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM", "debreak"],
        legend=None,
    )
    sns.lineplot(
        data=df,
        x="Precision",
        y="Recall",
        hue="Tool",
        style="Coverage",
        palette=colors,
        hue_order=["SVDSS", "cuteSV", "pbsv", "sniffles", "SVIM", "debreak"],
        markers=True,
        dashes=False,
        markersize=11,
        legend=None,
    )

    plt.xlim([50, 100])
    plt.ylim([20, 100])

    # plt.legend(framealpha=1)
    # plt.title("(b) Coverage titration on HG007")
    plt.tight_layout()
    plt.savefig("coverage.png", dpi=300)
    plt.savefig("coverage.svg", dpi=300)


if __name__ == "__main__":
    main()
