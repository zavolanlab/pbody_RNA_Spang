"""
Validation: global fraction F per cohort from several fracSanity runs vs the imaging P-body intensity
(pb/cyto, mean with bootstrap CI), reproducing the right panel of globalF_vs_pBodyIntensity in
pbody_RNA_Spang.current.ipynb (F +- 1.96 sqrt(var_F), clipped to [0, 1]; Spearman and Pearson r).

Usage: compare_globalF_with_imaging.py OUT_PREFIX LABEL=FITS_DIR [LABEL=FITS_DIR ...]
FITS_DIR holds one pickle per cohort (fit.<time>_<cond>.pkl or <time>__<cond>.pkl).
"""

import glob
import os
import pickle
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import bootstrap

from common import CONDITION_COLORS, TIME_MARKERS

INTENSITY = "/scicore/home/zavolan/GROUP/pbody_RNA_spang_zavolan_v2/materials/pbody_intensity.tsv"
KEYS = ["F", "F_CI_lower", "F_CI_upper", "var_F", "F_identifiable", "pi_outlier", "tau", "a", "b"]


def load_fits(fits_dir):
    rows = []
    for path in sorted(glob.glob(os.path.join(fits_dir, "*.pkl"))):
        m = re.search(r"(10min|60min)_+(\w+)\.pkl$", os.path.basename(path))
        if not m:
            continue
        with open(path, "rb") as fh:
            fit = pickle.load(fh)
        rows.append({"time_point": m.group(1), "condition": m.group(2), **{k: fit.get(k) for k in KEYS}})
    return pd.DataFrame(rows)


def intensity_summary():
    df = pd.read_csv(INTENSITY, sep="\t")
    rows = []
    for (tp, cond), x in df.groupby(["time_point", "condition"])["pb/cyto"]:
        x = x.dropna().to_numpy()
        ci = bootstrap((x,), np.mean, n_resamples=10000, confidence_level=0.95, method="percentile",
                       rng=np.random.default_rng(0)).confidence_interval
        rows.append({"time_point": tp, "condition": cond, "pb/cyto": x.mean(), "y_CI_lower": ci.low, "y_CI_upper": ci.high})
    return pd.DataFrame(rows)


def main():
    out_prefix, runs = sys.argv[1], [a.split("=", 1) for a in sys.argv[2:]]
    inten = intensity_summary()
    fig, axes = plt.subplots(1, len(runs), figsize=(4.2 * len(runs), 3.8), squeeze=False)
    tables = []
    for ax, (label, fits_dir) in zip(axes[0], runs):
        data = load_fits(fits_dir).merge(inten, on=["time_point", "condition"], how="inner")
        se = np.sqrt(data["var_F"].astype(float))
        lo, hi = np.clip(data["F"] - 1.96 * se, 0, 1), np.clip(data["F"] + 1.96 * se, 0, 1)
        ax.errorbar(data["F"], data["pb/cyto"], xerr=[data["F"] - lo, hi - data["F"]],
                    yerr=[data["pb/cyto"] - data["y_CI_lower"], data["y_CI_upper"] - data["pb/cyto"]],
                    fmt="none", ecolor="black", capsize=0.1, alpha=0.7)
        for _, r in data.iterrows():
            ax.scatter(r["F"], r["pb/cyto"], color=CONDITION_COLORS[r["condition"]], marker=TIME_MARKERS[r["time_point"]],
                       s=80, edgecolor="black", zorder=3)
        rho, p_rho = stats.spearmanr(data["F"], data["pb/cyto"])
        r, p_r = stats.pearsonr(data["F"], data["pb/cyto"])
        ax.set_title(f"{label}\nSpearman ρ = {rho:.2f} (p = {p_rho:.2g})\nPearson r = {r:.2f} (p = {p_r:.2g}), n = {len(data)}",
                     fontsize=9)
        ax.set(ylim=(0, 7), xlabel=r"Estimated Global Fraction ($\hat{F}$)")
        tables.append(data.assign(run=label))
    axes[0, 0].set_ylabel("normalized mean gray value\n(PB / cytoplasm)")
    handles = [plt.Line2D([], [], marker="o", ls="", color=CONDITION_COLORS[c], mec="black", label=c)
               for c in ["minusD", "DTT", "H2O2", "CCCP", "Tun"]]
    handles += [plt.Line2D([], [], marker=m, ls="", color="white", mec="black", label=t) for t, m in TIME_MARKERS.items()]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    fig.savefig(out_prefix + ".pdf", bbox_inches="tight", dpi=300)
    fig.savefig(out_prefix + ".png", bbox_inches="tight", dpi=100)
    pd.concat(tables).to_csv(out_prefix + ".tsv", sep="\t", index=False)
    print(pd.concat(tables).pivot_table(index=["time_point", "condition"], columns="run", values="F").round(4))


if __name__ == "__main__":
    main()
