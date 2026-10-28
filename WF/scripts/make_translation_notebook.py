"""Generates pbody_RNA_Spang.translation.ipynb: translation efficiency (public Ribo-seq) vs P-body recruitment."""
import sys
import nbformat as nbf

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


md(r"""
# Translation efficiency under stress vs P-body recruitment

Are mRNAs whose translation drops under a stress the ones recruited to P-bodies? Public yeast ribosome profiling (Ribo-seq) with matched RNA-seq, processed with `WF/Snakefile-Riboseq-public`, is compared with the P-body pull-down of this project (`pbody_RNA_Spang.fracSanity.ipynb`, variant `Sanity_utils_allTPM`).

| public dataset | stress, control | strain | our comparison |
|---|---|---|---|
| Zid & O'Shea, Nature 2014 (PRJNA223423) | glucose starvation (−GLU, ~15 min) vs +GLU, 3 replicates | BY4741 (WT) | minusD vs UT, 10 and 60 min |
| van Dalfsen et al., Dev Cell 2018 (PRJNA474715) | DTT, tunicamycin vs untreated, 2 replicates | strain "1362" (assumed wild type; check the GEO records) | DTT vs UT, Tun vs DMSO (60 min) |

**Processing.** Poly(A)-tailing libraries: 3' poly(A) trimmed, unique end-to-end STAR alignment to R64-1-1. Ribo-seq: P-site offsets per read length calibrated on annotated start codons (lengths 25–35 with a clear offset), P-sites counted in the CDS without its first 15 and last 5 codons; P-site bigWigs per strand (RPM) for IGV. RNA-seq: fragment midpoints in the same trimmed CDS.

**Translation efficiency with fracSanity.** Sanity normalizes the Ribo and RNA samples separately (log2 TPM over the trimmed CDS length), replicates are averaged per condition, and fracSanity fits $D_g = \langle A^{Ribo}_g\rangle - \langle A^{RNA}_g\rangle = \log_2\alpha_g - \log_2 F + \varepsilon_g$, i.e. Ribo-seq takes the place of the pull-down and RNA-seq that of Total. $\alpha_g$ is the translation efficiency up to the global factor $F$, which (as for the pull-down) is weakly identified. Comparisons between genes and conditions therefore use the **relative** log2 TE (centered $D_g$) and its change between stress and control, which do not depend on $F$; Spearman correlations are invariant to global shifts anyway.
""")

code(r"""
import warnings
warnings.simplefilter('ignore')

import pandas as pd
import numpy as np
import os
import sys
import glob
import csv

from dotenv import load_dotenv
from scipy import stats

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from IPython.display import display
from adjustText import adjust_text
%matplotlib inline
""")

code(r"""
load_dotenv("pbody_RNA_Spang.scicore.env", override=True)
subdirs = {k: os.getenv(v) for k, v in {"main_project_dir": "MAIN_PROJECT_DIR", "wf_dir": "WF_DIR", "figures_dir": "FIGURES_DIR",
                                       "tables_dir": "TABLES_DIR", "wf_runs_dir": "WF_RUNS_DIR"}.items()}

RIBO_OUT = os.path.join(subdirs["wf_runs_dir"], "riboseq_public", "output")
SU_RUN = "Sanity_utils_allTPM"
SU_TABLES = os.path.join(subdirs["tables_dir"], "gene_expression_analysis", "bulkRNAseq", SU_RUN)
TE_TABLES = os.path.join(subdirs["tables_dir"], "translation", "riboseq_public")
TE_FIGS = os.path.join(subdirs["figures_dir"], "translation", "riboseq_public")
for d in (TE_TABLES, TE_FIGS):
    os.makedirs(d, exist_ok=True)

from IPython.display import HTML

def figure_caption(path, formats=("png", "pdf")):
    # file name of a saved figure, shown with the plot; the links point to the full-size files relative to the
    # HTML report (WF/scripts/make_html_reports.py copies them to <report>/figures/, same subfolders as FIGURES_DIR)
    path = os.path.splitext(path)[0] if path.endswith((".png", ".pdf")) else path
    rel = os.path.relpath(path, subdirs["figures_dir"])
    links = " | ".join(f'<a href="figures/{rel}.{f}" target="_blank">{f}</a>' for f in formats)
    display(HTML(f'<p><b>Figure:</b> <code>{os.path.basename(path)}</code> &nbsp;({links})</p>'))

def savefig(fig, name, formats=("png", "pdf")):
    for fmt in formats:
        fig.savefig(os.path.join(TE_FIGS, f"{name}.{fmt}"), bbox_inches="tight", dpi=300)
    figure_caption(os.path.join(TE_FIGS, name), formats)

def save_tsv(df, name):
    df.to_csv(os.path.join(TE_TABLES, name), sep="\t", index=False, quoting=csv.QUOTE_NONE)

cat_color_dict = dict(zip(["UT", "DMSO", "minusD", "DTT", "H2O2", "CCCP", "Tun"],
                          ["grey", "skyblue", "orangered", "lime", "royalblue", "yellow", "fuchsia"]))
samples = pd.read_csv(os.path.join(subdirs["wf_runs_dir"], "riboseq_public", "samples.tsv"), sep="\t")
genes_df = pd.read_csv(os.path.join(subdirs["wf_runs_dir"], "v2", "output", SU_RUN, "inputs", "genes.tsv"), sep="\t", index_col=0)
gene_names = genes_df["gene_name"].fillna(pd.Series(genes_df.index, index=genes_df.index))
rp_genes = set(gene_names.index[gene_names.str.match(r"^RP[LS]\d+[AB]?$")])
display(samples.groupby(["dataset", "assay", "condition"]).size().unstack())
""")

md("# Quality control of the Ribo-seq data")

code(r"""
# P-site offsets per read length (used: clear offset), and the fraction of P-sites in frame 0 of CDSs
off = pd.concat([pd.read_csv(f, sep="\t").assign(sample=os.path.basename(f).split(".")[0])
                 for f in glob.glob(os.path.join(RIBO_OUT, "psites", "*.psite_offsets.tsv"))])
off = off.merge(samples[["sample", "dataset", "condition", "replicate"]], on="sample")
display(off.loc[off["used"]].pivot_table(index="length", columns=["dataset", "condition", "replicate"], values="offset"))

# mapping (STAR) and P-sites per sample
rows = []
for s in samples["sample"]:
    f = os.path.join(RIBO_OUT, "map", f"{s}.Log.final.out")
    if os.path.exists(f):
        st = dict(l.strip().split(" |\t") for l in open(f) if " |\t" in l)
        rows.append({"sample": s, "input reads, M": int(st["Number of input reads"]) / 1e6,
                     "% uniquely mapped": float(st["Uniquely mapped reads %"].rstrip("%"))})
mapping = samples.merge(pd.DataFrame(rows), on="sample")
cds = {s: pd.read_csv(os.path.join(RIBO_OUT, "psites", f"{s}.cds_counts.tsv"), sep="\t")["count"].sum() / 1e6
       for s in samples["sample"] if os.path.exists(os.path.join(RIBO_OUT, "psites", f"{s}.cds_counts.tsv"))}
mapping["CDS counts, M"] = mapping["sample"].map(cds)
display(mapping.drop(columns=["url", "md5"], errors="ignore").round(2))

sns.set(font_scale=1, style="ticks")
fig, ax = plt.subplots(figsize=(6, 3.5))
sns.pointplot(data=off.loc[off["used"]], x="length", y="frame0_frac", hue="dataset", errorbar=None, ax=ax)
ax.axhline(1 / 3, color="grey", ls="--", lw=0.8)
ax.set(ylabel="fraction of P-sites in frame 0", xlabel="footprint length (nt)", ylim=(0, 1))
savefig(fig, "QC.psite_frame0_by_length")
plt.show()

print("P-site bigWigs for IGV (RPM, one file per strand):")
for f in sorted(glob.glob(os.path.join(RIBO_OUT, "psites", "*.psites.*.rpm.bw"))):
    print(" ", f)
""")

md(r"""
# Translation efficiency per condition and its change under stress

`te_global_F.tsv`: $F$ per condition (weakly identified, see the caveats of the fracSanity notebook). Relative log2 TE: $D_g$ minus its median over genes in that condition.
""")

code(r"""
te = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(RIBO_OUT, "te", "*", "te_per_condition.tsv"))])
teF = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(RIBO_OUT, "te", "*", "te_global_F.tsv"))])
ted = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(os.path.join(RIBO_OUT, "te", "*", "te_differential.tsv"))])
display(teF.round(4))
te["RP"] = te["gene_id"].isin(rp_genes)
ted["RP"] = ted["gene_id"].isin(rp_genes)

sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
sns.violinplot(data=te, x="condition", y="log2_TE_relative", hue="RP", split=True, cut=0, inner="quart",
               ax=axes[0], palette={True: "orangered", False: "lightgrey"})
axes[0].set(title="relative log2 TE per condition", xlabel="")
d = ted.assign(label=ted["dataset"].str.split("_").str[0] + ": " + ted["comparison"].str.replace("_vs_", " vs "))
sns.violinplot(data=d, x="label", y="delta_log2_TE_relative", hue="RP", split=True, cut=0, inner="quart",
               ax=axes[1], palette={True: "orangered", False: "lightgrey"})
axes[1].axhline(0, color="black", lw=0.8)
axes[1].set(title="change of relative log2 TE (stress − control)", xlabel="")
axes[1].tick_params(axis="x", rotation=20)
savefig(fig, "TE_relative.per_condition_and_change")
plt.show()
save_tsv(ted, "TE_differential.all.tsv")
""")

md(r"""
# Translation efficiency change vs P-body recruitment

For each matched pair, per gene: change of relative log2 TE under the stress (public Ribo-seq) vs change of P-body recruitment in our data, as the change of $\log_2\alpha_g$ (fracSanity, `delta_log2_alpha`) and as the F-independent change of the centered $D_g = \log_2(\mathrm{PD}/\mathrm{Total})$. A negative correlation means that mRNAs whose translation drops are those recruited to P-bodies.
""")

code(r"""
alpha = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_alpha_per_cohort.tsv"), sep="\t")
alpha["D_centered"] = alpha["log2_PD_over_Total"] - alpha.groupby("cohort")["log2_PD_over_Total"].transform("median")
dalpha = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_differential_recruitment.stress_vs_control.tsv"), sep="\t")

pairs = [  # (dataset, TE comparison, our comparison)
    ("Zid2014_glucose", "minusD_vs_UT", "10min;minusD_vs_UT"),
    ("Zid2014_glucose", "minusD_vs_UT", "60min;minusD_vs_UT"),
    ("vanDalfsen2018_ER", "DTT_vs_UT", "10min;DTT_vs_UT"),
    ("vanDalfsen2018_ER", "DTT_vs_UT", "60min;DTT_vs_UT"),
    ("vanDalfsen2018_ER", "Tun_vs_UT", "10min;Tun_vs_DMSO"),
    ("vanDalfsen2018_ER", "Tun_vs_UT", "60min;Tun_vs_DMSO"),
]

def recruitment_change(comp):
    test, ref = comp.split("_vs_")
    tp = test.split(";")[0]
    D = alpha.pivot(index="gene_id", columns="cohort", values="D_centered")
    out = pd.DataFrame({"delta_D_centered": D[test] - D[f"{tp};{ref}"]})
    dd = dalpha.loc[dalpha["comparison"] == comp].set_index("gene_id")
    out["delta_log2_alpha"] = dd["delta_log2_alpha"]
    out["lfsr_alpha"] = dd["lfsr"]
    return out

merged, summary = [], []
for ds, te_comp, comp in pairs:
    t = ted.loc[(ted["dataset"] == ds) & (ted["comparison"] == te_comp)].set_index("gene_id")
    m = recruitment_change(comp).join(t[["delta_log2_TE_relative", "delta_log2_TE_relative_SE", "padj"]], how="inner").dropna(
        subset=["delta_log2_TE_relative", "delta_D_centered"])
    m["gene_name"] = m.index.map(gene_names)
    m["RP"] = m.index.isin(rp_genes)
    m["pair"] = f"{comp}  ~  {ds.split('_')[0]} {te_comp}"
    merged.append(m.reset_index())
    rho_D, p_D = stats.spearmanr(m["delta_log2_TE_relative"], m["delta_D_centered"])
    rho_a, p_a = stats.spearmanr(m["delta_log2_TE_relative"], m["delta_log2_alpha"], nan_policy="omit")
    summary.append({"our comparison": comp, "Ribo-seq": f"{ds} {te_comp}", "n genes": len(m),
                    "Spearman (ΔTE, ΔD centered)": rho_D, "p": p_D, "Spearman (ΔTE, Δlog2α)": rho_a,
                    "RP: median ΔTE": m.loc[m["RP"], "delta_log2_TE_relative"].median(),
                    "RP: median ΔD centered": m.loc[m["RP"], "delta_D_centered"].median()})
merged = pd.concat(merged)
summary = pd.DataFrame(summary)
save_tsv(merged, "TE_change_vs_recruitment_change.per_gene.tsv")
save_tsv(summary, "TE_change_vs_recruitment_change.summary.tsv")
display(summary.round(3))
""")

code(r"""
genes_of_interest_list = ["ACT1", "PGK1", "BSC1", "DCP2", "TPI1", "RLM1", "ATP11", "MET10", "RPS6A", "RPS9A", "RPS9B", "RPL19B"]
sns.set(font_scale=1, style="ticks")
n = len(pairs)
fig, axes = plt.subplots(2, int(np.ceil(n / 2)), figsize=(5 * int(np.ceil(n / 2)), 9.5), squeeze=False)
for ax, (pair, m) in zip(axes.flat, merged.groupby("pair", sort=False)):
    ax.hexbin(m["delta_log2_TE_relative"], m["delta_D_centered"], gridsize=60, bins="log", cmap="Greys", mincnt=1)
    rp = m.loc[m["RP"]]
    ax.scatter(rp["delta_log2_TE_relative"], rp["delta_D_centered"], s=8, color="orangered", label="RP mRNAs")
    # running median of the recruitment change in bins of the TE change
    b = pd.qcut(m["delta_log2_TE_relative"], 20, duplicates="drop")
    med = m.groupby(b, observed=True)[["delta_log2_TE_relative", "delta_D_centered"]].median()
    ax.plot(med["delta_log2_TE_relative"], med["delta_D_centered"], color="royalblue", lw=2, label="binned median")
    g = m.loc[m["gene_name"].isin(genes_of_interest_list)]
    texts = [ax.text(r["delta_log2_TE_relative"], r["delta_D_centered"], r["gene_name"], fontsize=7) for _, r in g.iterrows()]
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.6)
    rho, p = stats.spearmanr(m["delta_log2_TE_relative"], m["delta_D_centered"])
    ax.set(title=f"{pair}\nSpearman ρ = {rho:.2f} (p = {p:.1g}), n = {len(m)}",
           xlabel="Δ relative log2 TE (stress − control)", ylabel="Δ P-body recruitment\n(centered log2 PD/Total)")
    ax.title.set_fontsize(9)
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5))
axes.flat[0].legend(loc="upper right", fontsize=8)
for ax in axes.flat[n:]:
    ax.axis("off")
fig.tight_layout()
savefig(fig, "TE_change_vs_recruitment_change.scatter")
plt.show()
""")

code(r"""
# quadrants: genes with a strong change of relative TE (|Δ| > 1 and |Δ| > 2 SE; F-independent) and their
# recruitment change; Fisher's exact test: are TE-down genes more often among the most recruited (top 20% of ΔD)?
rows = []
for pair, m in merged.groupby("pair", sort=False):
    z = m["delta_log2_TE_relative"] / m["delta_log2_TE_relative_SE"]
    te_down = (m["delta_log2_TE_relative"] < -1) & (z < -2)
    te_up = (m["delta_log2_TE_relative"] > 1) & (z > 2)
    rec_up = m["delta_D_centered"] > m["delta_D_centered"].quantile(0.8)
    for name, sel in (("TE down", te_down), ("TE up", te_up)):
        tab = [[(sel & rec_up).sum(), (sel & ~rec_up).sum()], [(~sel & rec_up).sum(), (~sel & ~rec_up).sum()]]
        orat, p = stats.fisher_exact(tab)
        rows.append({"pair": pair, "genes": name, "n": int(sel.sum()), "fraction among top-20% recruited": tab[0][0] / max(1, sel.sum()),
                     "odds ratio": orat, "p": p, "median ΔD centered": m.loc[sel, "delta_D_centered"].median()})
quad = pd.DataFrame(rows)
save_tsv(quad, "TE_change_vs_recruitment_change.quadrants.tsv")
display(quad.round(3))
""")

md(r"""
# Basal translation and basal P-body presence

Unstressed conditions: relative log2 TE (public +GLU / untreated) vs the basal relative recruitment of our UT samples (centered $\log_2\alpha_g$ and centered $D_g$). Are poorly translated mRNAs more present in P-bodies already without stress?
""")

code(r"""
sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 4, figsize=(19, 4.2))
k = 0
for ds in sorted(te["dataset"].unique()):
    t = te.loc[(te["dataset"] == ds) & (te["condition"] == "UT")].set_index("gene_id")
    for coh in ("10min;UT", "60min;UT"):
        a = alpha.loc[alpha["cohort"] == coh].set_index("gene_id")
        m = t[["log2_TE_relative"]].join(a[["D_centered", "log2_alpha_mean"]], how="inner").dropna()
        ax = axes[k]; k += 1
        ax.hexbin(m["log2_TE_relative"], m["D_centered"], gridsize=60, bins="log", cmap="Greys", mincnt=1)
        rp = m.loc[m.index.isin(rp_genes)]
        ax.scatter(rp["log2_TE_relative"], rp["D_centered"], s=8, color="orangered")
        rho, p = stats.spearmanr(m["log2_TE_relative"], m["D_centered"])
        ax.set(title=f"{ds.split('_')[0]} untreated TE vs our {coh}\nSpearman ρ = {rho:.2f} (p = {p:.1g}), n = {len(m)}",
               xlabel="relative log2 TE", ylabel="relative P-body recruitment\n(centered log2 PD/Total)")
        ax.title.set_fontsize(9)
fig.tight_layout()
savefig(fig, "TE_basal_vs_recruitment_basal")
plt.show()
""")

md(r"""
# Caveats

- Different labs, strains, media and stress protocols (time points: glucose starvation ~15 min in Zid & O'Shea vs 10 / 60 min here; DTT / tunicamycin durations and doses of van Dalfsen et al. differ from ours; strain "1362" assumed wild type). Matched comparisons are approximate.
- Cycloheximide-treated Ribo-seq (Zid & O'Shea +/-GLU footprints with CHX) shows initiation / termination peaks; the first 15 and last 5 codons are excluded from the counts. Under acute glucose starvation, translation initiation collapses within minutes; ribosome footprints then partly reflect elongating ribosomes that remain.
- TE from fracSanity is defined up to the global factor $F$ (weakly identified); only relative (centered) TE and its changes, or rank correlations, are interpreted. P-body recruitment is likewise used as centered $\log_2(\mathrm{PD}/\mathrm{Total})$ changes; $\Delta\log_2\alpha_g$ is reported alongside (identical ranking within a comparison, shifted by $\Delta\log_2 F$).
- RNA-seq of the public datasets is counted on the same trimmed CDS as the P-sites; UTR reads do not enter TE.
- No public match for CCCP; H2O2 Ribo-seq (e.g. PRJEB25491, PRJNA635625) not processed yet.
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python (sanity_utils)", "language": "python", "name": "sanity_utils"},
                  "language_info": {"name": "python"}}
nbf.write(nb, sys.argv[1])
print("written", sys.argv[1], len(cells), "cells")
