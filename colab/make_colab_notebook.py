"""Generates colab/pbody_explorer.ipynb, the Colab notebook shared with the Spang lab."""
import sys
import nbformat as nbf

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


md(r"""
# P-body recruitment explorer (yeast PD vs Total, fracSanity)

Interactive view of the fracSanity results of the P-body pull-down (**PD**, Dcp2-HBH cCLAP) vs **Total** RNA-seq, for UT, DMSO, glucose starvation (minusD), DTT, H2O2, CCCP and tunicamycin (Tun) at 10 and 60 min.

- $\alpha_g$: fraction of the mRNA of gene $g$ that is in the P-body pull-down, with a 95% credible interval.
- $\Delta\log_2\alpha_g$: change of recruitment of a stress relative to its control (UT for minusD, DTT, H2O2 and DMSO; DMSO for CCCP and Tun). `lfsr` is the local false sign rate (probability that the sign of the change is wrong).
- $F$: fraction of all mRNA molecules in the pull-down.

**How to use.** Run the cells from top to bottom (*Runtime → Run all*). Change the gene list in section 2 and re-run the cells below it. The data are read from a SWITCHdrive folder that the Zavolan lab keeps in sync with the analysis on sciCORE; the first run downloads ~100 MB into this Colab session.
""")

md(r"""
## 1. Data access

Add two **secrets** (key icon in the left sidebar of Colab, enable *Notebook access*):
- `PBODY_SHARE_TOKEN`: the last part of the SWITCHdrive share link `https://drive.switch.ch/index.php/s/<token>`;
- `PBODY_SHARE_PASSWORD`: the password of the share link.
""")

code(r"""
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

try:
    from google.colab import userdata
    TOKEN = userdata.get("PBODY_SHARE_TOKEN")
    PASSWORD = userdata.get("PBODY_SHARE_PASSWORD")
except ImportError:  # outside Colab: environment variables
    TOKEN = os.environ["PBODY_SHARE_TOKEN"]
    PASSWORD = os.environ.get("PBODY_SHARE_PASSWORD", "")

# Nextcloud public-share WebDAV endpoints (new and legacy); the first that answers is used
ENDPOINTS = [
    (f"https://drive.switch.ch/public.php/dav/files/{TOKEN}/", ("anonymous", PASSWORD)),
    ("https://drive.switch.ch/public.php/webdav/", (TOKEN, PASSWORD)),
]
DATA_DIR = "/content/pbody_data" if os.path.isdir("/content") else "pbody_data"


def fetch(relpath, refresh=False):
    # download one shared file into DATA_DIR (cached) and return its local path
    local = os.path.join(DATA_DIR, relpath)
    if os.path.exists(local) and not refresh:
        return local
    os.makedirs(os.path.dirname(local), exist_ok=True)
    for base, auth in ENDPOINTS:
        r = requests.get(base + relpath, auth=auth, timeout=300)
        if r.status_code == 200:
            with open(local, "wb") as f:
                f.write(r.content)
            return local
    raise RuntimeError(f"could not download {relpath} (HTTP {r.status_code}); check token and password")


def read(relpath, **kw):
    return pd.read_csv(fetch(relpath), sep="\t", **kw)


alpha = read("fracSanity/fracSanity_alpha_per_cohort.tsv")
diff = read("fracSanity/fracSanity_differential_recruitment.stress_vs_control.tsv")
global_F = read("fracSanity/fracSanity_global_F.tsv")
features = read("gene_features/gene_features.tsv")
set_diff = read("gene_set_activities/differential_recruitment.tsv")
gene_sets = read("gene_set_activities/gene_sets.tsv")
print(f"{alpha['gene_id'].nunique()} genes, {alpha['cohort'].nunique()} cohorts, "
      f"{diff['comparison'].nunique()} stress-vs-control comparisons")
""")

code(r"""
CONDITION_ORDER = ["UT", "DMSO", "minusD", "DTT", "H2O2", "CCCP", "Tun"]
CONDITION_COLORS = dict(zip(CONDITION_ORDER, ["grey", "skyblue", "orangered", "lime", "royalblue", "yellow", "fuchsia"]))
TIME_POINTS = ["10min", "60min"]
TIME_MARKERS = {"10min": "o", "60min": "s"}
TAU = 0.5  # threshold on |delta log2 alpha| for "more / less recruited"

# gene names and systematic names both work
name_of = dict(zip(alpha["gene_id"], alpha["gene_name"].fillna(alpha["gene_id"])))
id_of = {**{g: g for g in name_of}, **{n: g for g, n in name_of.items()}}


def to_ids(genes):
    missing = [g for g in genes if g not in id_of]
    if missing:
        print("not in the gene universe:", missing)
    return [id_of[g] for g in genes if g in id_of]
""")

md(r"""
## 2. Genes of interest
""")

code(r"""
GENES = ["BSC1", "TPI1", "RLM1", "ATP11", "ILM1", "MRPL38", "AIM2", "MET10", "ACT1", "PGK1",
         "RPS9A", "RPS9B", "RPS6A", "RPL19B"]
gene_ids = to_ids(GENES)
""")

md(r"""
## 3. Fraction in P-bodies, $\alpha_g$, per condition and time point
""")

code(r"""
def plot_alpha(gene_ids, ncols=4, log=True):
    nrows = int(np.ceil(len(gene_ids) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows), squeeze=False, sharey=log)
    for ax, g in zip(axes.flat, gene_ids):
        d = alpha.loc[alpha["gene_id"] == g]
        for i, cond in enumerate(CONDITION_ORDER):
            for j, tp in enumerate(TIME_POINTS):
                r = d.loc[(d["condition"] == cond) & (d["time_point"] == tp)]
                if r.empty:
                    continue
                r = r.iloc[0]
                x = i + (j - 0.5) * 0.3
                ax.errorbar(x, r["alpha"], yerr=[[r["alpha"] - r["alpha_CI_lower"]], [r["alpha_CI_upper"] - r["alpha"]]],
                            fmt=TIME_MARKERS[tp], color=CONDITION_COLORS[cond], mec="black", ms=6, ecolor="black", capsize=2)
        ax.set_xticks(range(len(CONDITION_ORDER)), CONDITION_ORDER, rotation=60)
        ax.set_title(name_of[g])
        if log:
            ax.set_yscale("log")
    for ax in axes.flat[len(gene_ids):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\alpha_g$ (95% CI)")
    handles = [plt.Line2D([], [], marker=TIME_MARKERS[tp], ls="", color="white", mec="black", label=tp) for tp in TIME_POINTS]
    fig.legend(handles=handles, loc="upper right", frameon=False)
    fig.tight_layout()
    return fig


fig = plot_alpha(gene_ids)
""")

md(r"""
## 4. Change of recruitment, stress vs control

Heatmap of $\Delta\log_2\alpha_g$; a dot marks lfsr < 0.05 and $|\Delta\log_2\alpha_g| > \tau$.
""")

code(r"""
def plot_delta(gene_ids):
    d = diff.loc[diff["gene_id"].isin(gene_ids)]
    mat = d.pivot(index="gene_id", columns="comparison", values="delta_log2_alpha").reindex(gene_ids)
    sig = d.assign(s=(d["lfsr"] < 0.05) & (d["delta_log2_alpha"].abs() > TAU)) \
           .pivot(index="gene_id", columns="comparison", values="s").reindex(index=mat.index, columns=mat.columns)
    lim = np.nanmax(np.abs(mat.values))
    fig, ax = plt.subplots(figsize=(0.55 * mat.shape[1] + 2, 0.35 * mat.shape[0] + 1.5))
    sns.heatmap(mat.rename(index=name_of), cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                cbar_kws={"label": r"$\Delta\log_2\alpha_g$"})
    yy, xx = np.where(sig.fillna(False).values)
    ax.scatter(xx + 0.5, yy + 0.5, s=12, color="black")
    ax.set(xlabel="", ylabel="")
    fig.tight_layout()
    return fig, d


fig, delta_table = plot_delta(gene_ids)
delta_table[["gene_name", "comparison", "delta_log2_alpha", "delta_log2_alpha_CI_lower",
             "delta_log2_alpha_CI_upper", "lfsr"]].round(3)
""")

md(r"""
## 5. Global fraction $F$ of mRNA in the pull-down
""")

code(r"""
fig, ax = plt.subplots(figsize=(5, 3.2))
for j, tp in enumerate(TIME_POINTS):
    d = global_F.loc[global_F["time_point"] == tp].set_index("condition").reindex(CONDITION_ORDER)
    x = np.arange(len(CONDITION_ORDER)) + (j - 0.5) * 0.3
    ax.errorbar(x, d["F"], yerr=[d["F"] - d["F_CI_lower"], d["F_CI_upper"] - d["F"]], fmt="none", ecolor="black", capsize=2)
    ax.scatter(x, d["F"], marker=TIME_MARKERS[tp], c=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
               edgecolor="black", zorder=3, label=tp)
ax.set_xticks(range(len(CONDITION_ORDER)), CONDITION_ORDER, rotation=60)
ax.set(ylabel="$F$ (95% CI)")
ax.legend(frameon=False)
fig.tight_layout()
""")

md(r"""
## 6. Genes of a gene set, and gene sets of a gene

Gene-set activities: ISMARA-style fit of $\log_2(\mathrm{PD}/\mathrm{Total})$ on GO, Reactome, KEGG and WikiPathways sets; `delta_activity` is the change of the shared recruitment of the set's genes, stress vs control.
""")

code(r"""
QUERY = "cytoplasmic translation"   # part of a set name (case-insensitive)

hits = set_diff.loc[set_diff["name"].str.contains(QUERY, case=False, na=False)]
hits.pivot_table(index=["set", "name", "size"], columns="comparison", values="delta_activity").round(3)
""")

code(r"""
# genes of the first matching set, ready to be used as GENES in section 2
members = gene_sets.loc[gene_sets["name"].str.contains(QUERY, case=False, na=False)]
if len(members):
    print(members.iloc[0]["name"], "|", members.iloc[0]["genes"])
""")

code(r"""
GENE = "RPS9A"   # all gene sets that contain this gene
set_genes = gene_sets["genes"].str.split(",")
gene_sets.loc[set_genes.apply(lambda g: name_of.get(id_of.get(GENE), GENE) in g), ["set_id", "database", "name", "size"]]
""")

md(r"""
## 7. Sequence features of the selected genes
""")

code(r"""
features.set_index("gene_id").loc[[g for g in gene_ids if g in set(features["gene_id"])]]
""")

md(r"""
## 8. Save a figure or table

Files written to `/content` disappear when the session ends; download them from the file browser (folder icon), or save to your Google Drive with `from google.colab import drive; drive.mount('/content/drive')`.
""")

code(r"""
# fig.savefig("my_genes.pdf", bbox_inches="tight")
# delta_table.to_csv("my_genes.delta_log2_alpha.tsv", sep="\t", index=False)
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                  "language_info": {"name": "python"}}
nbf.write(nb, sys.argv[1])
print("written", sys.argv[1], len(cells), "cells")
