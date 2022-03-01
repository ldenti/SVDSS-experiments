import sys
import os
import glob


def parse_summary(fpath):
    TP, FP, FN = 0, 0, 0
    for line in open(fpath):
        line = line.strip("\n \t\"")
        if line.startswith("TP-base\""):
            TP = int(line.split(" ")[1][:-1])
        elif line.startswith("FP\""):
            FP = int(line.split(" ")[1][:-1])
        elif line.startswith("FN\""):
            FN = int(line.split(" ")[1][:-1])
    return TP, FP, FN


def main():
    indir = sys.argv[1]

    print("Mode,Tool,TP,FP,FN,P,R,F1")

    for f in glob.glob(os.path.join(indir, "*", "*", "summary.txt")):
        tool = f.split("/")[-3]
        mode = f.split("/")[-2]
        tp, fp, fn = parse_summary(f)
        P = round(tp/(tp+fp)*100, 1)
        R = round(tp/(tp+fn)*100, 1)
        F = round(2*(P*R)/(P+R) if P+R != 0 else 0.0, 1)
        print(mode, tool, tp, fp, fn, P, R, F, sep = ",")


if __name__ == "__main__":
    main()
