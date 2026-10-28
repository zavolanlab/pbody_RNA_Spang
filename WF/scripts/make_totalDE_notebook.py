"""Generates pbody_RNA_Spang.TotalDE.ipynb: differential expression in the Total fraction (cells below)."""
import sys
import nbformat as nbf

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


md(r"""
# Differential expression in the Total fraction (Sanity_utils)

Total RNA-seq of the yeast stress time course (Spang lab): untreated (**UT**), **DMSO**, glucose starvation (**minusD**), **DTT**, **H2O2**, **CCCP**, tunicamycin (**Tun**), 10 and 60 min, 3 replicates. Companion of `pbody_RNA_Spang.fracSanity.ipynb` (P-body pull-down vs Total) and `pbody_RNA_Spang.current.ipynb` (same environment file, folder layout and plotting style).

- Normalization: C++ Sanity on the Total samples, from the Sanity_utils pipeline `WF/Snakefile-fracSanity` (variant `SU_RUN` below: all expressed genes, length correction, log2 TPM), and ISMARA replicate averages per time point and condition (`method='reml'`).
- Differential expression ([docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/differential_analysis.md)): a gene is differentially expressed if the true change of its log2 relative expression exceeds $\theta$ in absolute value; reported are $\log_2$FC, its SE, the p-value of the threshold test, BH-adjusted `padj` and the local false sign rate `lfsr`. Calls: `padj` < 0.05 with $\theta = 0.5$.
- Comparisons: each stress vs its control (UT for minusD, DTT, H2O2 and DMSO; DMSO for CCCP and Tun), per time point; and 60 vs 10 min within each condition.
- Gene sets ([docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/gene_set_analysis.md)): ISMARA-style activities of GO, Reactome, KEGG and WikiPathways sets fitted on the per-sample log2 expression, with replicate averages and contrasts; CAMERA-like competitive test on the gene-level z-scores as a baseline.

**Kernel.** "Python (sanity_utils)". The gene-set fit takes several minutes; run the notebook on a compute node (sbatch).
""")

code(r"""
import warnings
warnings.simplefilter('ignore')

# general purpose packages
import pandas as pd
import numpy as np
import os
import sys
import csv
import pickle

from dotenv import load_dotenv

import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import display
%matplotlib inline

# Sanity_utils: condition averages, differential expression, gene-set analysis
import sanity_utils as su
from sanity_utils.genesets import load_gene_sets, set_source_table

# general-purpose plots
from zavolab_pyutils import pca_plot
from adjustText import adjust_text
""")

code(r"""
# 1. Load the environment variables
load_dotenv("pbody_RNA_Spang.scicore.env", override=True)

# 2. Reconstruct the subdirs dictionary
subdirs = {
    "main_project_dir": os.getenv("MAIN_PROJECT_DIR"),
    "wf_dir": os.getenv("WF_DIR"),
    "figures_dir": os.getenv("FIGURES_DIR"),
    "tables_dir": os.getenv("TABLES_DIR"),
    "wf_runs_dir": os.getenv("WF_RUNS_DIR"),
}

# outputs of the Sanity_utils pipeline (WF/Snakefile-fracSanity) and of this notebook
WF_version = "v2"
SU_RUN = "Sanity_utils_allTPM"   # all expressed genes, log2 TPM
SU_OUT = os.path.join(subdirs["wf_runs_dir"], WF_version, "output", SU_RUN)
DE_TABLES = os.path.join(subdirs["tables_dir"], "gene_expression_analysis", "bulkRNAseq", SU_RUN, "TotalDE")
DE_FIGS = os.path.join(subdirs["figures_dir"], "gene_expression_analysis", "bulkRNAseq", SU_RUN, "TotalDE")
GENESETS_DIR = os.path.join(subdirs["main_project_dir"], "external_data", "genesets")
for d in (DE_TABLES, DE_FIGS):
    os.makedirs(d, exist_ok=True)

# comparisons (stress vs control, 60 vs 10 min) shared with the pipeline
sys.path.insert(0, os.path.join(subdirs["wf_dir"], "scripts", "sanity_utils_pipeline"))
from common import comparisons

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
        fig.savefig(os.path.join(DE_FIGS, f"{name}.{fmt}"), bbox_inches="tight", dpi=300)
    figure_caption(os.path.join(DE_FIGS, name), formats)

def save_tsv(df, name, index=False):
    df.to_csv(os.path.join(DE_TABLES, name), sep="\t", header=True, index=index, quoting=csv.QUOTE_NONE)

print("Environment loaded and directories verified.")
""")

md("# Define metadata labels and colors uniformly")

code(r"""
def get_metadata_with_labels_and_colors(input_df):
    sel_files = input_df.copy()

    def get_time_point(x):
        if x['merged_sample'].split('_')[1].startswith("26025"):
            return "60min"
        else:
            return "10min"

    present_conditions = {
        "D": "minusD",
        "min": "minusD",
        "DTT": "DTT",
        "U": "UT",
        "Untreat": "UT",
        "CCCP": "CCCP",
        "H2O2": "H2O2",
        "DMSO": "DMSO",
        "Tun": "Tun",
        "Tunica": "Tun",
    }

    def get_condition(x):
        if x['time_point'] == '10min':
            cond_raw = x["merged_sample"].split("-")[5].split("_")[0]
        else:
            cond_raw = x["merged_sample"].split("-")[2].split("_")[0]
        if cond_raw in present_conditions.keys():
            return present_conditions[cond_raw]
        cur_cond = ""
        for elem in list(present_conditions.keys()):
            if cond_raw.startswith(elem) and len(elem) > len(cur_cond):
                cur_cond = elem
        return present_conditions[cur_cond]

    def get_replicate(x):
        unadj = x["merged_sample"].split("_S")[0]
        if x['time_point'] == '10min':
            if (x["condition"] in ["DMSO", "UT"]) and x["fraction"] == "PD":
                return "R" + unadj.split("-")[-2][-1]
            else:
                return "R" + unadj[-1]
        elif x['time_point'] == '60min':
            return 'R' + unadj.split('-')[-1]

    sel_files["time_point"] = sel_files.apply(lambda x: get_time_point(x), 1)
    sel_files["fraction"] = sel_files.apply(lambda x: 'PD' if '-PD-' in x["merged_sample"] else 'Total', 1)
    sel_files["condition"] = sel_files.apply(lambda x: get_condition(x), 1)
    sel_files["replicate"] = sel_files.apply(lambda x: get_replicate(x), 1)

    sel_files["sample_label_within_experiment"] = (
        sel_files["time_point"] + ";" + sel_files["fraction"] + ";" + sel_files["condition"] + ";" + sel_files["replicate"]
    )

    cat = "condition"
    cat_order = ["UT", "DMSO", "minusD", "DTT", "H2O2", "CCCP", "Tun"]
    cat_palette = ["grey", "skyblue", "orangered", "lime", "royalblue", "yellow", "fuchsia"]
    cat_color_dict = {}
    for k, cat_val in enumerate(cat_order):
        cat_color_dict[cat_val] = cat_palette[k]
    sel_files[cat + ";color"] = sel_files[cat].map(cat_color_dict)

    mrk = "time_point"
    mrk_order = ["10min", "60min"]
    mrk_palette = ["o", "s"]
    mrk_marker_dict = {}
    for k, mrk_val in enumerate(mrk_order):
        mrk_marker_dict[mrk_val] = mrk_palette[k]
    sel_files[mrk + ";marker"] = sel_files[mrk].map(mrk_marker_dict)

    sel_files = sel_files.sort_values(
        ["fraction", "condition", "replicate", "merged_sample", "sample"]
    ).reset_index(drop=True)
    return sel_files, cat_color_dict, mrk_marker_dict


start_samples = pd.read_csv(os.path.join(subdirs["wf_runs_dir"], WF_version, "start_samples.tsv"), delimiter="\t")
metadata_df, cat_color_dict, mrk_marker_dict = get_metadata_with_labels_and_colors(start_samples)
metadata_df = metadata_df.drop_duplicates("merged_sample").reset_index(drop=True)
metadata_df["sample"] = metadata_df["merged_sample"]

# don't use replicates 1-3 for CCCP Total 60min (as in pbody_RNA_Spang.current.ipynb)
metadata_df = metadata_df.loc[~((metadata_df['time_point'] == '60min') & (metadata_df['fraction'] == 'Total') &
                                (metadata_df['condition'] == 'CCCP') & (metadata_df['replicate'].isin(['R1', 'R2', 'R3'])))].reset_index(drop=True)
metadata_df["cohort"] = metadata_df["time_point"] + ";" + metadata_df["condition"]
total_metadata = metadata_df.loc[metadata_df["fraction"] == "Total"].reset_index(drop=True)
print(total_metadata.groupby(["time_point", "condition"]).size().unstack())
""")

md(r"""
# Sanity normalization of the Total fraction

Sanity output of the pipeline (rule `sanity_fraction`): per-sample log2 TPM with error bars, and ISMARA replicate averages per cohort (time point;condition).
""")

code(r"""
with open(os.path.join(SU_OUT, "sanity", "Total", "sanity.pkl"), "rb") as fh:
    sanity_total = pickle.load(fh)
res, avg = sanity_total["result"], sanity_total["averages"]
genes_df = pd.read_csv(os.path.join(SU_OUT, "inputs", "genes.tsv"), sep="\t", index_col=0)
gene_names = genes_df["gene_name"].fillna(pd.Series(genes_df.index, index=genes_df.index))

samples_list = [s for s in total_metadata["sample"] if s in res.log2_expression.columns]
assert len(samples_list) == len(total_metadata), "samples of the metadata missing in the Sanity output"
print(res.log2_expression.shape, res.params.get("units"), "| cohorts:", list(avg.means.columns))
""")

md("# PCA")

code(r"""
counts_df = pd.read_csv(os.path.join(SU_OUT, "inputs", "counts.tsv"), sep="\t", index_col=0)
expressed_everywhere_genes = list(counts_df.loc[counts_df[samples_list].min(axis=1) > 0].index)  # at least one read in all samples

for time_point in [None] + list(mrk_marker_dict):
    sub_meta = total_metadata if time_point is None else total_metadata.loc[total_metadata["time_point"] == time_point]
    sub_samples = list(sub_meta["sample"])
    PCA_UMAP_subset = res.log2_expression.loc[expressed_everywhere_genes, sub_samples].reset_index(drop=True)
    label = "all" if time_point is None else time_point
    print(label, len(PCA_UMAP_subset), "genes")
    savefig_path = os.path.join(DE_FIGS, f"PCA_gene_expression.Total.{label}.png")
    pca_plot(
        PCA_UMAP_subset,
        sub_samples,
        sub_meta.reset_index(drop=True),
        "condition",
        savefig_path,
        sns_color_palette=cat_color_dict,   # dict: {condition: color}
        style_feature="time_point",
        markers=mrk_marker_dict,            # dict: {time_point: marker}
        calculate_permanova_R2=True,
        permanova_feature="cohort",
        permanova_R2_ajusted=False,
        add_text_labels_for_samples=True,
        text_label_column_for_samples="replicate",
        s_param=40,
        figsize=(4.2, 4.2),
        legend_title=f"Total fraction ({label})\nconditions:",
        add_2D_KDE_countours=False,
    )
    plt.show()
    figure_caption(savefig_path, formats=("png",))
    PCA_UMAP_subset.to_csv(savefig_path.replace(".png", ".tsv"), sep="\t", header=True, index=None, quoting=csv.QUOTE_NONE)
""")

md(r"""
# Differential expression

Threshold test of $|\Delta| > \theta$ on the ISMARA replicate averages; the fraction of differentially expressed genes per comparison is estimated from the distribution of all $\log_2$FC and their errors (`estimate_de_fraction`), independently of any p-value cutoff.
""")

code(r"""
LFC_THRESHOLD = 0.5   # theta, log2
PADJ = 0.05
comparisons_df = comparisons()

de_results, de_summary = {}, []
for _, row in comparisons_df.iterrows():
    de = su.differential_expression(avg, row["test"], row["reference"], lfc_threshold=LFC_THRESHOLD)
    de["A"] = (de[f"mean_{row['test']}"] + de[f"mean_{row['reference']}"]) / 2
    de["call"] = np.where(de["padj"] < PADJ, np.where(de["log2FC"] > 0, "up", "down"), "n.s.")
    de_results[row["label"]] = de
    frac = su.estimate_de_fraction(de, n_threads=8)
    de_summary.append({"comparison": row["label"], "kind": row["kind"], "test": row["test"], "reference": row["reference"],
                       "n_up": int((de["call"] == "up").sum()), "n_down": int((de["call"] == "down").sum()),
                       "estimated_fraction_DE": frac.fraction_de})
de_summary_df = pd.DataFrame(de_summary)
save_tsv(de_summary_df, f"TotalDE_summary.theta{LFC_THRESHOLD}.tsv")
display(de_summary_df.round(3))

de_long = pd.concat([d.assign(comparison=k) for k, d in de_results.items()]).rename_axis("gene_id").reset_index()
de_long.insert(1, "gene_name", de_long["gene_id"].map(gene_names))
cols = ["gene_id", "gene_name", "comparison", "A", "log2FC", "SE", "z_score", "p_value", "padj", "lfsr", "call"]
save_tsv(de_long[cols], f"TotalDE.theta{LFC_THRESHOLD}.tsv")
""")

code(r"""
# number of differentially expressed genes per comparison
stress_summary = de_summary_df.loc[de_summary_df["kind"] == "stress_vs_control"].copy()
stress_summary["condition"] = stress_summary["test"].str.split(";").str[1]
stress_summary["time_point"] = stress_summary["test"].str.split(";").str[0]
sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
for ax, tp in zip(axes, mrk_marker_dict):
    d = stress_summary.loc[stress_summary["time_point"] == tp].set_index("condition")
    d = d.reindex([c for c in cat_color_dict if c in d.index])
    x = np.arange(len(d))
    ax.bar(x, d["n_up"], color=[cat_color_dict[c] for c in d.index], edgecolor="black")
    ax.bar(x, -d["n_down"], color=[cat_color_dict[c] for c in d.index], edgecolor="black", hatch="//")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, [f"{c}\nvs {d.loc[c, 'reference'].split(';')[1]}" for c in d.index], fontsize=9)
    ax.set_title(tp)
axes[0].set_ylabel(f"# genes up / down\n(padj < {PADJ}, θ = {LFC_THRESHOLD})")
sns.despine()
savefig(fig, f"n_DE_genes.theta{LFC_THRESHOLD}")
""")

md(r"""
## MA plots

$A$: mean log2 TPM of the two conditions; genes with `padj` < 0.05 in color. Labelled: genes discussed by the Spang lab (Wang et al. eLife 2018; thesis slides).
""")

code(r"""
genes_of_interest_list = ["ACT1", "PGK1", "BSC1", "DCP2", "TPI1", "RLM1", "ATP11", "MET10", "RPS6A", "RPS9A", "RPS9B", "RPL19B"]
name_to_id = {n: g for g, n in gene_names.items()}
goi_ids = [name_to_id[n] for n in genes_of_interest_list if n in name_to_id]
call_colors = {"up": "firebrick", "down": "steelblue", "n.s.": "lightgrey"}


def ma_plot(ax, de, title, title_color):
    # returns the gene labels; they are placed with adjust_text once all panels are drawn
    for call in ("n.s.", "down", "up"):
        d = de.loc[de["call"] == call]
        ax.scatter(d["A"], d["log2FC"], s=3 if call == "n.s." else 5, color=call_colors[call], lw=0, rasterized=True,
                   label=f"{call} ({len(d)})")
    for y in (0, LFC_THRESHOLD, -LFC_THRESHOLD):
        ax.axhline(y, color="black", lw=0.8 if y == 0 else 0.5, ls="-" if y == 0 else "--")
    g = de.loc[de.index.intersection(goi_ids)]
    ax.scatter(g["A"], g["log2FC"], s=14, facecolor="none", edgecolor="black", lw=0.8)
    ax.set_title(title, fontsize=10, color="black", bbox=dict(facecolor=title_color, alpha=0.35, lw=0))
    ax.legend(loc="lower right", fontsize=7, markerscale=2, frameon=False, handletextpad=0.1)
    return [ax.text(r["A"], r["log2FC"], gene_names[gid], fontsize=6) for gid, r in g.iterrows()]


sns.set(font_scale=1, style="ticks")
for kind in ("stress_vs_control", "time"):
    rows = comparisons_df.loc[comparisons_df["kind"] == kind]
    ncols = 6 if kind == "stress_vs_control" else 7
    nrows = int(np.ceil(len(rows) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.9 * nrows), sharex=True, sharey=True, squeeze=False)
    labels = []
    for ax, (_, row) in zip(axes.flat, rows.iterrows()):
        cond = row["test"].split(";")[1]
        labels.append((ax, ma_plot(ax, de_results[row["label"]], row["label"].replace("_vs_", " vs "), cat_color_dict[cond])))
    for ax in axes.flat[len(rows):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("A (mean log2 TPM)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\log_2$FC")
    fig.tight_layout()
    for ax, texts in labels:
        if texts:
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5))
    savefig(fig, f"MA_plots.{kind}.theta{LFC_THRESHOLD}")
    plt.show()
""")


md(r"""
# Gene-set analysis

## ISMARA-style gene-set activities of the stress responses

Because the 10 and 60 min samples form separate batches, the activities are fitted not on the log2 expression of the samples but on **per-replicate log2 fold changes of each stress relative to its own control at the same time point** (UT for minusD, DTT, H2O2 and DMSO; DMSO for CCCP and Tun): replicate $i$ of the stress minus replicate $i$ of the control, i.e. pairs from the same experiment (60 min CCCP, replicates R4–R6, is paired with DMSO R1–R3 in replicate order). For each pair column $s$,

$$
\log_2\mathrm{FC}_{gs} = \tilde c_s + \sum_m N_{gm} A_{ms} + \text{noise},
$$

with $N_{gm} = 1$ if gene $g$ belongs to set $m$ (Balwierz et al. 2014; [docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/gene_set_analysis.md)). Only the column constant $\tilde c_s$ is removed (`center="columns"`, which also absorbs global shifts such as those from the rRNA content); gene constants are not removed, so responses shared by all stresses (e.g. repression of ribosomal-protein genes) are kept. $A_{ms}$ is the shift of the log2FC shared by the genes of set $m$ beyond what the other sets of the same fit explain; the ridge penalty $\lambda$ is cross-validated. Replicate pairs are averaged per comparison (`reml`), giving $\langle A_m\rangle$, its error, z-score, BH `padj` and lfsr per comparison; $z_m = \sqrt{\frac1S\sum_s (A_{ms}/\delta A_{ms})^2}$ (step 5 of the docs, ISMARA eq. 14) measures the overall significance of a set across all pairs.

**One fit per database.** Sets of different databases overlap strongly (e.g. GO:BP *cytoplasmic translation*, GO:CC *cytosolic ribosome*, KEGG *Ribosome* and `LIT:RP`). In a joint fit they share their common signal and each gets a smaller, less significant activity. The activities are therefore fitted separately for each database: literature marker sets (`LIT:`, section "Literature validation"), GO:BP, GO:CC, GO:MF, Reactome, KEGG and WikiPathways, each with 5–500 genes in the gene universe (identical sets merged within a database).
- **Biotype sets in every fit.** The biotype sets (`BIOTYPE:`, relative to protein-coding genes) are included in every fit as covariates, because whole RNA classes shift strongly between samples (see the TIN section). They contain almost no genes of the database sets, so their activities hardly depend on the fit; the combined tables report them from the literature fit.
- **Multiple testing.** `padj` is BH over the sets of the same fit (database + biotypes) per comparison.
- **Ranking.** $z_m$ is computed per fit; the ridge penalty $\lambda$ is cross-validated per fit. Overlapping sets of different databases can now all appear among the top sets, so the top-$z_m$ list contains near-duplicates of the same program.
""")

code(r"""
from sanity_utils.genesets import GeneSetCollection

# per-replicate log2 fold changes, stress vs its control at the same time point (replicate i vs replicate i)
sc_comparisons = comparisons_df.loc[comparisons_df["kind"] == "stress_vs_control"]
E = res.log2_expression
lfc_cols, lfc_meta = {}, []
for _, row in sc_comparisons.iterrows():
    t = total_metadata.loc[total_metadata["cohort"] == row["test"]].sort_values("replicate")
    r = total_metadata.loc[total_metadata["cohort"] == row["reference"]].sort_values("replicate")
    for (_, ts), (_, rs) in zip(t.iterrows(), r.iterrows()):
        pid = f"{row['label']}|{ts['replicate']}-{rs['replicate']}"
        lfc_cols[pid] = E[ts["sample"]] - E[rs["sample"]]
        lfc_meta.append({"sample": pid, "comparison": row["label"], "pair": f"{ts['replicate']} vs {rs['replicate']}"})
lfc_values = pd.DataFrame(lfc_cols)
lfc_meta = pd.DataFrame(lfc_meta)
print(lfc_values.shape[1], "replicate pairs:", lfc_meta.groupby("comparison")["pair"].apply(list).to_dict())

# marker sets of stress-response programs from the literature (gene names -> systematic IDs), see "Literature validation"
name_to_id = {n: g for g, n in gene_names.items()}
LIT_MARKERS = {
    "LIT:iESR": ["HSP12", "HSP26", "HSP42", "HSP104", "CTT1", "TPS1", "TPS2", "TSL1", "GSY2", "PGM2", "GPH1", "DDR2",
                 "GRE2", "XBP1", "YGP1"],                                          # Gasch et al. 2000, induced ESR
    "LIT:Yap1_oxidative": ["TRX2", "TRR1", "GSH1", "GLR1", "CCP1", "TSA2", "SRX1", "GPX2", "AHP1"],  # Godon 1998; Lee 1999
    "LIT:UPR_Hac1": ["KAR2", "PDI1", "ERO1", "LHS1", "FKB2", "EUG1", "DER1", "JEM1", "SCJ1"],          # Travers et al. 2000
    "LIT:glucose_repressed": ["HXT2", "HXT4", "HXT6", "HXT7", "SUC2", "ADH2", "JEN1", "FBP1", "PCK1", "ICL1", "MLS1"],  # DeRisi 1997
}
LIT_NAMES = {"LIT:iESR": "induced environmental stress response (iESR) markers",
             "LIT:Yap1_oxidative": "Yap1 oxidative stress response targets",
             "LIT:UPR_Hac1": "unfolded protein response, Hac1 targets",
             "LIT:glucose_repressed": "glucose-repressed genes (derepressed upon glucose depletion)",
             "LIT:RP": "cytoplasmic ribosomal-protein genes (repressed ESR)"}
lit_sets = {k: {name_to_id[n] for n in v if n in name_to_id} for k, v in LIT_MARKERS.items()}
lit_sets["LIT:RP"] = set(gene_names.index[gene_names.str.match(r"^RP[LS]\d+[AB]?$")])
print("marker genes not in the gene universe:", {k: [n for n in v if n not in name_to_id] for k, v in LIT_MARKERS.items()})
lit_coll = GeneSetCollection({k: v & set(lfc_values.index) for k, v in lit_sets.items()}, LIT_NAMES, "literature")

# gene biotypes of the annotation (Ensembl GTF), protein-coding genes as the reference (the biotypes partition the
# genes, so one of them must be left out); plus the mitochondrially encoded genes
universe = set(lfc_values.index)
biotype_sets = {f"BIOTYPE:{b}": set(g.index) & universe for b, g in genes_df.groupby("gene_biotype") if b != "protein_coding"}
biotype_sets["BIOTYPE:mitochondrially_encoded"] = set(genes_df.index[genes_df["chrom"] == "Mito"]) & universe
biotype_sets = {k: v for k, v in biotype_sets.items() if len(v) >= 5}
biotype_coll = GeneSetCollection(biotype_sets, {k: k.split(":")[1].replace("_", " ") + " (biotype)" for k in biotype_sets}, "biotype")
print("biotype sets:", {k: len(v) for k, v in biotype_sets.items()})

# one ISMARA fit per database (sets of different databases overlap and would share their signal in a joint fit),
# each with the biotype sets as covariates (they absorb shifts of whole RNA classes vs protein-coding genes);
# padj: BH over the sets of the same fit
GENESET_DATABASES = ["literature", "GO:BP", "GO:CC", "GO:MF", "Reactome", "KEGG", "WikiPathways"]
db_colls = {db: lit_coll if db == "literature" else
            load_gene_sets(GENESETS_DIR, sources=(db,), genes=lfc_values.index, min_size=5, max_size=500, drop_duplicates=True)
            for db in GENESET_DATABASES}

fits, act_diff, info_all, zm_all, pair_acts = {}, [], [], [], []
for db, c_db in db_colls.items():
    coll_db = c_db + biotype_coll
    a = su.fit_set_activities(lfc_values, coll_db, center="columns", ridge="cv", cv_repeats=5)
    a_avg = a.average_conditions(lfc_meta, sample_col="sample", cond_col="comparison", method="reml")
    fits[db] = {"activities": a, "averages": a_avg, "collection": coll_db}
    info = set_source_table(coll_db).assign(fit=db)
    info["size"] = info["set_id"].map(coll_db.sizes)
    info_all.append(info)
    zm_all.append(info.assign(z_m=info["set_id"].map(a.z_values)))
    pair_acts.append(a.activities.assign(fit=db))
    for _, row in sc_comparisons.iterrows():
        d = su.activity_contrast(a_avg, {row["label"]: 1.0})
        act_diff.append(d[["delta_activity", "SE", "z_score", "p_value", "padj", "lfsr"]]
                        .assign(comparison=row["label"], kind=row["kind"], fit=db).rename_axis("set_id").reset_index())
    print(f"{db}: {len(c_db)} sets (+ {len(biotype_coll)} biotypes); lambda = {a.lambda_:.4g}; "
          f"median explained variance per pair {a.explained_variance.median():.3f}")

# combined tables: database sets from their own fit; biotype sets (fitted in every fit) from the literature fit
keep = lambda df: df.loc[~df["set_id"].str.startswith("BIOTYPE:") | (df["fit"] == "literature")].reset_index(drop=True)
act_diff_all = pd.concat(act_diff)
info_all = pd.concat(info_all)
act_diff_df = keep(act_diff_all).merge(keep(info_all)[["set_id", "database", "name", "size"]], on="set_id", how="left")
set_info = keep(info_all)
zm_df = keep(pd.concat(zm_all)).sort_values("z_m", ascending=False).reset_index(drop=True)
act_pairs = keep(pd.concat(pair_acts).rename_axis("set_id").reset_index()).set_index("set_id").drop(columns="fit")

# check: the biotype activities hardly depend on which database is fitted alongside them
bt_fits = act_diff_all.loc[act_diff_all["set_id"].str.startswith("BIOTYPE:")].pivot_table(
    index=["set_id", "comparison"], columns="fit", values="z_score")
print("biotype z-scores across fits: max range", round(float((bt_fits.max(axis=1) - bt_fits.min(axis=1)).max()), 2),
      "| median range", round(float((bt_fits.max(axis=1) - bt_fits.min(axis=1)).median()), 2))

save_tsv(act_diff_df, "TotalDE.gene_set_activities.log2FC_vs_control.tsv")
save_tsv(zm_df, "TotalDE.gene_set_activities.z_m.tsv")
with open(os.path.join(DE_TABLES, "TotalDE.gene_set_activities.pkl"), "wb") as fh:
    pickle.dump({"fits": fits, "log2FC_pairs": lfc_meta}, fh, protocol=pickle.HIGHEST_PROTOCOL)
display(zm_df.head(15).round(2))
""")


md(r"""
### Sets changed in at least one stress, and in all stresses

Activity z-scores of the stress-vs-control comparisons (comparisons over time are confounded with the 10 / 60 min batch and left out). `database` gives the source of each set (GO:BP, GO:CC, GO:MF, Reactome, KEGG, WikiPathways; `+` for identical sets merged from several sources).
""")

code(r"""
# (1) sets with padj < 0.05 in at least one stress-vs-control comparison, clustered
sc = act_diff_df.loc[act_diff_df["kind"] == "stress_vs_control"].copy()
sc_order = [c for c in comparisons_df["label"] if c in set(sc["comparison"])]
set_label = (set_info.set_index("set_id")["name"].str.slice(0, 60) + " (" + set_info.set_index("set_id")["database"] + ")")
print(sc.loc[sc["padj"] < 0.05].groupby("comparison").size().reindex(sc_order).rename("n sets padj < 0.05").to_string())

sig_any = sc.loc[sc["padj"] < 0.05, "set_id"].unique()
mat_any = sc.loc[sc["set_id"].isin(sig_any)].pivot(index="set_id", columns="comparison", values="z_score")[sc_order]
save_tsv(mat_any.rename_axis("set_id").reset_index().merge(set_info, on="set_id"),
         "TotalDE.gene_set_activities.z.sig_in_any_stress.tsv")
print(len(mat_any), "sets with padj < 0.05 in at least one comparison;",
      mat_any.index.map(set_info.set_index("set_id")["database"]).value_counts().to_dict())

lim = np.nanpercentile(np.abs(mat_any.values), 98)
sns.set(font_scale=0.8, style="white")
g = sns.clustermap(mat_any.rename(index=set_label).fillna(0), cmap="RdBu_r", center=0, vmin=-lim, vmax=lim,
                   col_cluster=False, method="average", metric="correlation",
                   figsize=(0.45 * len(sc_order) + 7, min(0.13 * len(mat_any) + 3, 60)),
                   yticklabels=len(mat_any) <= 400, cbar_kws={"label": "activity z-score"}, dendrogram_ratio=(0.12, 0.02))
g.ax_heatmap.set(xlabel="", ylabel="")
g.savefig(os.path.join(DE_FIGS, "gene_set_activities.sig_in_any_stress.clustermap.pdf"), bbox_inches="tight")
g.savefig(os.path.join(DE_FIGS, "gene_set_activities.sig_in_any_stress.clustermap.png"), bbox_inches="tight", dpi=200)
figure_caption(os.path.join(DE_FIGS, "gene_set_activities.sig_in_any_stress.clustermap"))
plt.show()
""")

code(r"""
# (2) sets with padj < 0.05 in all five stresses vs their controls (DMSO vs UT not required), per time point;
# 'same_direction' marks sets changed in the same direction in all five
stresses = ["minusD", "DTT", "H2O2", "CCCP", "Tun"]
for tp in mrk_marker_dict:
    cols = [c for c in sc_order if c.startswith(tp) and c.split(";")[1].split("_vs_")[0] in stresses]
    d = sc.loc[sc["comparison"].isin(cols)]
    padj = d.pivot(index="set_id", columns="comparison", values="padj")[cols]
    z = d.pivot(index="set_id", columns="comparison", values="z_score")[cols]
    sig_all = padj.index[(padj < 0.05).all(axis=1)]
    if len(sig_all) == 0:
        print(f"{tp}: no set with padj < 0.05 in all of {cols}")
        continue
    mat_all = z.loc[sig_all]
    same_dir = np.sign(mat_all).nunique(axis=1) == 1
    mat_all = mat_all.loc[mat_all.mean(axis=1).sort_values(ascending=False).index]
    save_tsv(mat_all.assign(same_direction=same_dir).rename_axis("set_id").reset_index().merge(set_info, on="set_id"),
             f"TotalDE.gene_set_activities.z.sig_in_all_stresses.{tp}.tsv")
    print(f"{tp}: {len(mat_all)} sets with padj < 0.05 in all five stresses, {int(same_dir.sum())} in the same direction")
    lim = np.nanpercentile(np.abs(mat_all.values), 98)
    sns.set(font_scale=0.9, style="white")
    fig, ax = plt.subplots(figsize=(0.5 * len(cols) + 7, min(0.28 * len(mat_all) + 1.5, 60)))
    sns.heatmap(mat_all.rename(index=lambda s: ("" if same_dir[s] else "≠ ") + set_label[s]), cmap="RdBu_r", center=0,
                vmin=-lim, vmax=lim, ax=ax, cbar_kws={"label": "activity z-score"}, yticklabels=len(mat_all) <= 200)
    ax.set(xlabel="", ylabel="", title=f"{tp}: gene sets changed in all stresses (≠: not all in the same direction)")
    savefig(fig, f"gene_set_activities.sig_in_all_stresses.{tp}.heatmap")
    plt.show()
""")

code(r"""
# (3) sets prioritized by their overall significance z_m across all replicate pairs (ISMARA eq. 14), clustered
N_TOP = 60
top_zm = zm_df.head(N_TOP)["set_id"]
mat_zm = sc.loc[sc["set_id"].isin(top_zm)].pivot(index="set_id", columns="comparison", values="z_score")[sc_order]
zm_label = {s: f"{set_label[s]}  z_m={z:.1f}" for s, z in zm_df.set_index("set_id")["z_m"].items()}
lim = np.nanpercentile(np.abs(mat_zm.values), 98)
sns.set(font_scale=0.8, style="white")
g = sns.clustermap(mat_zm.rename(index=zm_label).fillna(0), cmap="RdBu_r", center=0, vmin=-lim, vmax=lim,
                   col_cluster=False, method="average", metric="correlation",
                   figsize=(0.45 * len(sc_order) + 8, 0.16 * len(mat_zm) + 3),
                   cbar_kws={"label": "activity z-score"}, dendrogram_ratio=(0.12, 0.02))
g.ax_heatmap.set(xlabel="", ylabel="")
g.fig.suptitle(f"top {N_TOP} gene sets by overall significance z_m", y=1.01)
g.savefig(os.path.join(DE_FIGS, "gene_set_activities.top_z_m.clustermap.pdf"), bbox_inches="tight")
g.savefig(os.path.join(DE_FIGS, "gene_set_activities.top_z_m.clustermap.png"), bbox_inches="tight", dpi=200)
figure_caption(os.path.join(DE_FIGS, "gene_set_activities.top_z_m.clustermap"))
plt.show()
""")

code(r"""
# systematic shifts of gene biotypes (relative to protein-coding genes) in each stress vs its control:
# ISMARA activity (left; fitted jointly with all other sets) and median log2FC of the biotype's genes (right, gene level)
bt = sc.loc[sc["set_id"].str.startswith("BIOTYPE:")]
bt_z = bt.pivot(index="set_id", columns="comparison", values="z_score")[sc_order]
bt_d = bt.pivot(index="set_id", columns="comparison", values="delta_activity")[sc_order]
bt_padj = bt.pivot(index="set_id", columns="comparison", values="padj")[sc_order]
annot = bt_d.round(2).astype(str) + np.where(bt_padj < 0.05, "*", "")
biotype_of = genes_df["gene_biotype"].where(genes_df["chrom"] != "Mito", "mitochondrially_encoded")
med = pd.DataFrame({label: de_results[label]["log2FC"].groupby(biotype_of.reindex(de_results[label].index)).median()
                    for label in sc_order})
n_genes = biotype_of.reindex(de_results[sc_order[0]].index).value_counts()
med.index = [f"{b} (n={n_genes.get(b, 0)})" for b in med.index]
sns.set(font_scale=0.9, style="white")
fig, axes = plt.subplots(1, 2, figsize=(2 * (0.7 * len(sc_order) + 3), 0.5 * max(len(bt_z), len(med)) + 2.5))
lim = np.nanmax(np.abs(bt_z.values))
sns.heatmap(bt_z.rename(index=lambda s: s.split(":")[1]), cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=annot.values,
            fmt="", annot_kws={"fontsize": 7}, cbar_kws={"label": "activity z-score"}, ax=axes[0])
axes[0].set(xlabel="", ylabel="", title="ISMARA activity of gene biotypes vs protein-coding\n(numbers: activity = shift of log2FC; * padj < 0.05)")
lim = np.nanmax(np.abs(med.values))
sns.heatmap(med, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".2f", annot_kws={"fontsize": 7},
            cbar_kws={"label": "median log2FC"}, ax=axes[1])
axes[1].set(xlabel="", ylabel="", title="median log2FC of the genes of each biotype")
fig.tight_layout()
savefig(fig, "gene_biotypes.activity_and_median_log2FC")
save_tsv(bt.merge(set_info[["set_id"]], on="set_id"), "TotalDE.gene_biotype_activities.tsv")
plt.show()
""")


md(r"""
### RNA integrity (TIN) of the Total samples and the gene-biotype activities

The biotype sets (tRNA, snoRNA, mitochondrially encoded RNAs, transposable elements) have the highest $z_m$. Besides biology, differences in the integrity of the input RNA between the paired samples could produce them. RIN was not measured for all samples, so the **transcript integrity number** (TIN; [Wang et al. 2016, BMC Bioinformatics](https://doi.org/10.1186/s12859-016-0922-z)) is computed from the alignments instead: for each transcript, the coverage at 100 equally spaced positions of its exons is turned into $\mathrm{TIN} = 100 \cdot \exp(H)/n$, with $H$ the Shannon entropy of the coverage distribution over the $n$ covered positions. TIN is 100 for uniform coverage and decreases as degradation makes coverage uneven along the transcript; transcripts with ≤ 10 distinct read starts get TIN = 0.
- Pipeline: `WF/Snakefile-TIN` (config `WF/config.TIN.yaml`) runs [zavolanlab/tin-score-calculation](https://github.com/zavolanlab/tin-score-calculation) v0.6.3 (a parallel fork of RSeQC `tin.py`) on the UMI-deduplicated genome alignments of each Total sample, with the transcripts of the Ensembl R64-1-1.109 GTF.
- **Sample TIN** = median TIN over nuclear protein-coding genes with TIN > 0 in **all** samples. Genes with TIN = 0 are excluded because that value only means low coverage, so including them would make the median depend on sequencing depth.
- Caveats:
  - TIN is not strand-specific, so overlapping antisense transcripts contribute.
  - pysam's pileup caps the counted depth at 8,000 reads per position, which flattens the coverage of the most highly covered genes and raises their TIN. As a robustness check, the sample TIN is also computed without the top quartile of expression.
  - For short, structured RNAs (tRNAs, snoRNAs), coverage reflects library preparation (fragmentation, reverse-transcription stops at modified nucleotides) more than integrity.

Tests:
1. **Does the integrity differ between conditions and batches?** (sample TIN per cohort)
2. **Do the per-pair activities follow the per-pair ΔTIN** (TIN of the stress replicate minus TIN of its paired control)? The correlation is computed over all 36 replicate pairs and also **within comparisons**, after removing each comparison's mean. The between-comparison part can be biological (a stress can change mRNA decay), but replicate-to-replicate co-variation within a comparison points to sample handling.
3. **Length bias:** the slope of the log2FC of protein-coding genes on log2 transcript length in each pair, vs ΔTIN.
""")

code(r"""
# TIN per gene and sample (WF/Snakefile-TIN); sample TIN = median over nuclear protein-coding genes with TIN > 0 in all samples
TIN_DIR = os.path.join(subdirs["wf_runs_dir"], WF_version, "output", "TIN")
tin = pd.read_csv(os.path.join(TIN_DIR, "TIN.merged.tsv"), sep="\t", index_col=0)
tin_genes = pd.read_csv(os.path.join(TIN_DIR, "annotation", "genes.tsv"), sep="\t", index_col=0)
tin = tin[[s for s in total_metadata["sample"] if s in tin.columns]]
print("samples without TIN:", sorted(set(total_metadata["sample"]) - set(tin.columns)))

biotype_tin = tin_genes["gene_biotype"].where(tin_genes["chrom"] != "Mito", "mitochondrially_encoded")
pc_common = tin.index[(biotype_tin.reindex(tin.index) == "protein_coding") & (tin > 0).all(axis=1)]
# robustness: without the top quartile of expression (pileup depth cap of 8,000 reads per position)
mean_expr = res.log2_expression[tin.columns].mean(axis=1).reindex(pc_common)
pc_low = mean_expr.index[mean_expr <= mean_expr.quantile(0.75)]
tin_sample = (total_metadata.set_index("sample").loc[tin.columns, ["time_point", "condition", "replicate", "cohort"]]
              .assign(TIN=tin.loc[pc_common].median(), TIN_low_expr=tin.loc[pc_low].median(),
                      n_genes_TIN_pos=(tin > 0).sum()))
print(f"{len(pc_common)} protein-coding genes with TIN > 0 in all {tin.shape[1]} samples "
      f"(Spearman of the two sample TINs: {tin_sample['TIN'].corr(tin_sample['TIN_low_expr'], method='spearman'):.2f})")
save_tsv(tin_sample.rename_axis("sample").reset_index(), "TotalDE.TIN_per_sample.tsv")

# median TIN per biotype and sample (genes with TIN > 0 in that sample)
tin_bt = tin.where(tin > 0).groupby(biotype_tin.reindex(tin.index)).median()
tin_bt = tin_bt.loc[(tin > 0).groupby(biotype_tin.reindex(tin.index)).sum().median(axis=1) >= 5]

cond_order = [c for c in cat_color_dict if c in set(tin_sample["condition"])]
sns.set(font_scale=0.9, style="white")
fig, axes = plt.subplots(1, 2, figsize=(16, 4.5), gridspec_kw={"width_ratios": [1, 1.6]})
ax = axes[0]
sns.stripplot(data=tin_sample, x="time_point", y="TIN", hue="condition", hue_order=cond_order, palette=cat_color_dict,
              dodge=True, size=7, edgecolor="black", linewidth=0.5, ax=ax)
ax.set(xlabel="", ylabel="sample TIN (median, protein-coding)", title="RNA integrity of the Total samples")
ax.legend(bbox_to_anchor=(1, 1), loc="upper left", fontsize=8, frameon=False)
s_order = tin_sample.sort_values(["time_point", "condition", "replicate"],
                                 key=lambda c: c.map(cond_order.index) if c.name == "condition" else c).index
col_lab = tin_sample.loc[s_order].apply(lambda r: f"{r.time_point} {r.condition} {r.replicate}", axis=1)
sns.heatmap(tin_bt[s_order].set_axis(col_lab, axis=1), cmap="viridis", annot=False, cbar_kws={"label": "median TIN"},
            xticklabels=True, ax=axes[1])
axes[1].tick_params(axis="x", labelsize=7)
axes[1].set(xlabel="", ylabel="", title="median TIN per biotype (genes with TIN > 0)")
fig.tight_layout()
savefig(fig, "TIN.per_sample")
plt.show()
display(tin_sample.groupby(["time_point", "condition"])["TIN"].agg(["mean", "min", "max"]).round(1).unstack(0))
""")

code(r"""
# per replicate pair (same pairing as the activity fit): delta TIN = TIN(stress replicate) - TIN(control replicate)
pairs = []
for _, row in sc_comparisons.iterrows():
    t = total_metadata.loc[total_metadata["cohort"] == row["test"]].sort_values("replicate")
    r = total_metadata.loc[total_metadata["cohort"] == row["reference"]].sort_values("replicate")
    for (_, ts), (_, rs) in zip(t.iterrows(), r.iterrows()):
        pairs.append({"pair_id": f"{row['label']}|{ts['replicate']}-{rs['replicate']}", "comparison": row["label"],
                      "condition": ts["condition"], "time_point": ts["time_point"],
                      "dTIN": tin_sample.loc[ts["sample"], "TIN"] - tin_sample.loc[rs["sample"], "TIN"]})
pairs = pd.DataFrame(pairs).set_index("pair_id")
A = act_pairs[pairs.index]   # per-pair activities, each set from its own database fit

# length bias per pair: slope of log2FC on log2 exonic length, nuclear protein-coding genes
pc = lfc_values.index[lfc_values.index.isin(tin_genes.index[biotype_tin == "protein_coding"])]
log_len = np.log2(tin_genes.loc[pc, "exonic_length"])
pairs["length_slope"] = [np.polyfit(log_len, lfc_values.loc[pc, p], 1)[0] for p in pairs.index]

def corr_within(x, y, groups):
    # correlation after removing the mean of each comparison (replicate-to-replicate co-variation)
    xc, yc = x - x.groupby(groups).transform("mean"), y - y.groupby(groups).transform("mean")
    return xc.corr(yc)

top_sets = list(dict.fromkeys(zm_df.head(40)["set_id"].tolist() + [s for s in A.index if s.startswith("BIOTYPE:")]))
tin_cor = pd.DataFrame({"r_all": [A.loc[s].corr(pairs["dTIN"]) for s in top_sets],
                        "r_within": [corr_within(A.loc[s], pairs["dTIN"], pairs["comparison"]) for s in top_sets]},
                       index=top_sets)
# share of the across-pair variance of an activity that is within comparisons, and the part of it explained by delta TIN
# (the rest, between comparisons, is the stress effect that the averaged activities report)
within_share = pd.Series({s: (A.loc[s] - A.loc[s].groupby(pairs["comparison"]).transform("mean")).var() / A.loc[s].var()
                          for s in top_sets})
tin_cor["within_share"] = within_share
tin_cor["var_explained_by_dTIN"] = tin_cor["r_within"] ** 2 * within_share
tin_cor = tin_cor.join(zm_df.set_index("set_id")[["database", "name", "z_m"]])
save_tsv(tin_cor.rename_axis("set_id").reset_index(), "TotalDE.gene_set_activities.vs_dTIN.tsv")

# permutation null of r_within: delta TIN shuffled within comparisons
rng = np.random.default_rng(0)
null = []
for _ in range(2000):
    perm = pairs["dTIN"].groupby(pairs["comparison"]).transform(lambda v: rng.permutation(v.values))
    null.append(corr_within(A.loc[top_sets[0]], perm, pairs["comparison"]))
print(f"|r_within| expected under no association (95% quantile, 3 pairs per comparison): {np.quantile(np.abs(null), 0.95):.2f}")
r_len = pairs["length_slope"].corr(pairs["dTIN"]); r_len_w = corr_within(pairs["length_slope"], pairs["dTIN"], pairs["comparison"])
print(f"length slope vs delta TIN: r = {r_len:.2f} (within comparisons {r_len_w:.2f})")
display(tin_cor.loc[[s for s in top_sets if s.startswith("BIOTYPE:")] + ["LIT:iESR"],
                    ["z_m", "r_all", "r_within", "within_share", "var_explained_by_dTIN"]].round(2))
""")

code(r"""
# (a) biotype activities vs delta TIN per replicate pair; (b) correlation of the top-z_m sets with delta TIN; (c) length bias
bt_sets = [s for s in top_sets if s.startswith("BIOTYPE:")]
markers = {"10min": "o", "60min": "s"}
sns.set(font_scale=0.85, style="white")
ncol = 4
nrow = int(np.ceil((len(bt_sets) + 1) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.4 * nrow), squeeze=False)
for ax, s in zip(axes.flat, bt_sets + ["length_slope"]):
    y = pairs["length_slope"] if s == "length_slope" else A.loc[s]
    for (cond, tp), d in pairs.groupby(["condition", "time_point"]):
        ax.scatter(d["dTIN"], y[d.index], color=cat_color_dict.get(cond, "grey"), marker=markers.get(tp, "o"),
                   edgecolor="black", linewidth=0.4, s=40, label=f"{cond} {tp}")
    r = y.corr(pairs["dTIN"]); rw = corr_within(y, pairs["dTIN"], pairs["comparison"])
    title = "log2FC vs log2 length (slope)" if s == "length_slope" else s.split(":")[1]
    ax.set(title=f"{title}\nr = {r:.2f}, within comparisons {rw:.2f}", xlabel="ΔTIN (stress − control)",
           ylabel="slope" if s == "length_slope" else "activity")
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
for ax in axes.flat[len(bt_sets) + 1:]:
    ax.axis("off")
h, l = axes.flat[0].get_legend_handles_labels()
fig.legend(h, l, bbox_to_anchor=(1.0, 0.5), loc="center left", fontsize=7, frameon=False)
fig.tight_layout()
savefig(fig, "TIN.biotype_activities_vs_dTIN")
plt.show()

d = tin_cor.sort_values("z_m", ascending=True)
lab = d["name"].str.slice(0, 55) + " (" + d["database"] + ")"
fig, ax = plt.subplots(figsize=(8, 0.22 * len(d) + 1.5))
y = np.arange(len(d))
ax.barh(y + 0.2, d["r_all"], height=0.4, color="lightgrey", edgecolor="black", lw=0.3, label="all pairs")
ax.barh(y - 0.2, d["r_within"], height=0.4, color=np.where(d.index.str.startswith("BIOTYPE:"), "tab:red", "tab:blue"),
        label="within comparisons")
ax.set_yticks(y, lab, fontsize=7)
ax.axvline(0, color="black", lw=0.5)
ax.set(xlabel="Pearson r of the per-pair activity with ΔTIN", xlim=(-1, 1),
       title="top sets by $z_m$ (top) and biotype sets: association with RNA integrity")
ax.legend(loc="lower right", fontsize=8)
fig.tight_layout()
savefig(fig, "TIN.top_zm_sets_vs_dTIN")
plt.show()
""")

md(r"""
**Result (run of 2026-10-27).**
- **RNA integrity is high and uniform.**
  - Sample TIN ranges from 79.7 to 84.3 across all 42 Total samples.
  - The only systematic difference is a batch offset: about 82.7 at 10 min vs about 80.6 at 60 min. It cancels in the comparisons, which are all within a time point.
  - Within a time point, no condition differs from its control by more than the spread of its replicates.
  - Paired ΔTIN values are small (−2.6 to +2.7).
  - Removing the most highly expressed quarter of the genes (pileup depth cap) gives the same ranking (Spearman 1.00).
- **ΔTIN does not explain the biotype activities.**
  - For tRNA, snoRNA, mitochondrially encoded RNAs and Ty elements, 70–87% of the across-pair variance of the activity lies *between* comparisons. That is, the shifts are reproducible across the three replicate pairs of a stress.
  - ΔTIN explains ≤ 3% of the total variance for each of these sets, 0% for Ty elements and 13% for rRNA.
  - Example: the strong tRNA decrease in 60 min DTT, H2O2 and DMSO vs UT (activity −1.5 to −2) occurs at ΔTIN ≈ 0.
- **ΔTIN does track part of the replicate-to-replicate noise.**
  - Within comparisons, ΔTIN correlates with the length bias of the log2FC (r = 0.58).
  - Several broad sets also correlate with ΔTIN within comparisons (|r_within| 0.5–0.67; the 95% permutation bound is ≈ 0.38): rRNA, vacuole, protein maturation, cell wall, hydrolases.
  - For those sets this is a small part of the total variance (≤ 13%), and the averaged activities are dominated by the stress effect.
- **What TIN cannot rule out.** TIN measures how uniform the coverage is along mRNAs. It does not measure other technical factors that act on whole biotypes, such as the efficiency of rRNA depletion, retention of small RNAs (tRNAs, snoRNAs) during size selection, or the mature/precursor ratio of tRNAs. The biotype activities therefore remain to be interpreted with care. The likely explanation is no longer RNA degradation but either biology (Pol III/Maf1 repression, RiBi repression, Ty1 induction) or library composition. Of note, the Ty elements have a lower TIN at 60 than at 10 min (heatmap), but this is a batch-wide difference that does not differ between stresses and controls.
""")





code(r"""
# focus sets: ribosome, ribosome biogenesis, TOR and stress-related programs
focus_patterns = {
    "cytoplasmic translation (GO)": r"^cytoplasmic translation$",
    "ribosome biogenesis (GO)": r"^ribosome biogenesis$",
    "rRNA processing (GO)": r"^rRNA processing$",
    "Cytoplasmic ribosomal proteins (WP)": r"^Cytoplasmic ribosomal proteins$",
    "TOR signaling (GO)": r"^TOR signaling$",
    "proteasome": r"^proteasome complex$|^Proteasome$",
    "response to oxidative stress (GO)": r"^response to oxidative stress$",
    "unfolded protein response": r"endoplasmic reticulum unfolded protein response|^Unfolded Protein Response",
    "glycolytic process (GO)": r"^glycolytic process$",
    "aerobic respiration (GO)": r"^aerobic respiration$",
}
focus_ids = {}
for label, pat in focus_patterns.items():
    hit = set_info.loc[set_info["name"].str.contains(pat, case=False, regex=True, na=False)]
    if len(hit):
        focus_ids[label] = hit.sort_values("size", ascending=False)["set_id"].iloc[0]
print({k: v for k, v in focus_ids.items()})

sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(11, 0.45 * len(focus_ids) + 1.5), sharey=True)
for ax, tp in zip(axes, mrk_marker_dict):
    for i, (label, sid) in enumerate(focus_ids.items()):
        d = sc.loc[(sc["set_id"] == sid) & sc["comparison"].str.startswith(tp)]
        for j, (_, r) in enumerate(d.iterrows()):
            cond = r["comparison"].split(";")[1].split("_vs_")[0]
            y = i + (j - (len(d) - 1) / 2) * 0.12
            ax.errorbar(r["delta_activity"], y, xerr=1.96 * r["SE"], fmt=mrk_marker_dict[tp], ms=5,
                        color=cat_color_dict[cond], mec="black", mew=0.4, ecolor=cat_color_dict[cond], lw=0.8,
                        label=r["comparison"].split(";")[1].replace("_vs_", " vs ") if i == 0 else None)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(range(len(focus_ids)), list(focus_ids))
    ax.set(title=tp, xlabel=r"$\Delta$ activity (stress − control) ± 95% CI")
axes[0].invert_yaxis()
axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
sns.despine()
savefig(fig, "gene_set_activities.focus_sets")
plt.show()
""")

md(r"""
## Literature validation: do the steady-state transcriptome changes reproduce known stress responses?

Expected direction of the activity of stress-response programs in each stress vs its control, from transcriptome studies of *S. cerevisiae*:

| program | gene sets used | expected | references |
|---|---|---|---|
| repressed environmental stress response (rESR): ribosomal proteins, ribosome biogenesis | `LIT:RP`, GO:BP *cytoplasmic translation*, *ribosome biogenesis*, *rRNA processing* | down in glucose starvation, H2O2, DTT; by analogy in Tun and CCCP | Gasch et al. 2000 ([MBoC](https://doi.org/10.1091/mbc.11.12.4241)); Causton et al. 2001 ([MBoC](https://doi.org/10.1091/mbc.12.2.323)); Zid & O'Shea 2014 ([Nature](https://doi.org/10.1038/nature13578)) |
| induced ESR (chaperones, trehalose / glycogen metabolism, antioxidants) | `LIT:iESR`, GO:BP *trehalose biosynthetic process*, *protein refolding* | up in glucose starvation, H2O2, DTT; by analogy in Tun and CCCP | Gasch et al. 2000; Causton et al. 2001; Zid & O'Shea 2014 |
| oxidative stress response (Yap1/Skn7 regulon) | `LIT:Yap1_oxidative`, GO:BP *cellular response to oxidative stress* | up in H2O2 | Godon et al. 1998 ([JBC](https://doi.org/10.1074/jbc.273.35.22480)); Lee et al. 1999 ([JBC](https://doi.org/10.1074/jbc.274.23.16040)); Gasch et al. 2000 |
| unfolded protein response (Ire1–Hac1) | `LIT:UPR_Hac1`, GO:BP *endoplasmic reticulum unfolded protein response*, *ERAD pathway* | up in DTT and Tun (stronger at 60 min) | Travers et al. 2000 ([Cell](https://doi.org/10.1016/S0092-8674(00)80835-1)) |
| glucose derepression (Mig1 targets, gluconeogenesis, glyoxylate cycle, high-affinity hexose transporters) | `LIT:glucose_repressed`, GO:BP *gluconeogenesis*, *glyoxylate cycle* | up in glucose starvation | DeRisi et al. 1997 ([Science](https://doi.org/10.1126/science.278.5338.680)); Zid & O'Shea 2014 |
| none (vehicle) | all of the above | no change in DMSO vs UT | negative control |

"By analogy": the general stress response (ESR) is expected for any acute stress (Gasch et al. 2000), but no matched transcriptome reference was used for that stress. The P-body imaging, translation and TORC1 data of the same stresses are in Mookherjee, Rommel et al. (bioRxiv 2026; `materials/2026.03.18.712576v1.full.pdf`): strong translational shutoff in glucose starvation, H2O2 and CCCP, weak in DTT and Tun.

Calls per gene set and comparison (activity z-score of the comparison, BH `padj` over all sets): **✓** expected direction and padj < 0.05; **·** expected direction, not significant; **✗** opposite direction and padj < 0.05; **○** opposite direction, not significant; for DMSO vs UT, ✓ = no significant change, ✗ = significant change.
""")

code(r"""
EXPECTED = {   # program -> stress -> (expected sign, basis)
    "rESR": {"minusD": (-1, "direct"), "H2O2": (-1, "direct"), "DTT": (-1, "direct"), "Tun": (-1, "analogy"), "CCCP": (-1, "analogy")},
    "iESR": {"minusD": (1, "direct"), "H2O2": (1, "direct"), "DTT": (1, "direct"), "Tun": (1, "analogy"), "CCCP": (1, "analogy")},
    "oxidative": {"H2O2": (1, "direct")},
    "UPR": {"DTT": (1, "direct"), "Tun": (1, "direct")},
    "glucose derepression": {"minusD": (1, "direct")},
}
PROGRAM_SETS = {   # program -> literature marker set and GO proxies (name regex, matched within merged names)
    "rESR": ["LIT:RP", r"cytoplasmic translation", r"ribosome biogenesis", r"rRNA processing"],
    "iESR": ["LIT:iESR", r"trehalose biosynthetic process", r"protein refolding"],
    "oxidative": ["LIT:Yap1_oxidative", r"cellular response to oxidative stress"],
    "UPR": ["LIT:UPR_Hac1", r"endoplasmic reticulum unfolded protein response", r"ERAD pathway"],
    "glucose derepression": ["LIT:glucose_repressed", r"gluconeogenesis", r"glyoxylate cycle"],
}
rows = []
for prog, items in PROGRAM_SETS.items():
    for it in items:
        if it.startswith("LIT:"):
            ids = [it] if it in set(set_info["set_id"]) else []
        else:
            hit = set_info.loc[set_info["database"].str.startswith("GO:BP") &
                               set_info["name"].str.contains(rf"(?:^|\| ){it}(?: \||$)", regex=True)]
            ids = hit.sort_values("size", ascending=False)["set_id"].head(1).tolist()
        for sid in ids:
            rows.append({"program": prog, "set_id": sid})
lit_rows = pd.DataFrame(rows).drop_duplicates().merge(set_info[["set_id", "database", "name", "size"]], on="set_id")
print(lit_rows.groupby("program")["name"].apply(lambda x: [n[:50] for n in x]).to_dict())


def call(expected, z, padj):
    if expected == 0:
        return "✓" if padj >= 0.05 else "✗"
    same = np.sign(z) == expected
    return ("✓" if same else "✗") if padj < 0.05 else ("·" if same else "○")


obs = lit_rows.merge(sc[["set_id", "comparison", "z_score", "padj", "delta_activity"]], on="set_id")
obs["stress"] = obs["comparison"].str.split(";").str[1].str.split("_vs_").str[0]
obs["time_point"] = obs["comparison"].str.split(";").str[0]
exp = []
for _, r in obs.iterrows():
    if r["stress"] == "DMSO":
        exp.append((0, "negative control"))
    else:
        exp.append(EXPECTED[r["program"]].get(r["stress"], (np.nan, "")))
obs["expected"] = [e[0] for e in exp]
obs["basis"] = [e[1] for e in exp]
obs["call"] = [call(e, z, p) if not np.isnan(e) else "" for e, z, p in zip(obs["expected"], obs["z_score"], obs["padj"])]
save_tsv(obs, "TotalDE.literature_validation.tsv")

# heatmap of the activity z-scores with the calls
obs["row"] = obs["program"] + " | " + obs["name"].str.slice(0, 45) + " (" + obs["database"] + ")"
row_order = obs.drop_duplicates("row").set_index("row").loc[:, "program"].sort_values(
    key=lambda s: s.map(list(PROGRAM_SETS).index), kind="stable").index
zmat = obs.pivot(index="row", columns="comparison", values="z_score").loc[row_order, sc_order]
cmat = obs.pivot(index="row", columns="comparison", values="call").loc[row_order, sc_order]
lim = np.nanpercentile(np.abs(zmat.values), 98)
sns.set(font_scale=0.9, style="white")
fig, ax = plt.subplots(figsize=(0.75 * len(sc_order) + 6, 0.42 * len(zmat) + 2))
sns.heatmap(zmat, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=cmat, fmt="", annot_kws={"fontsize": 11},
            cbar_kws={"label": "activity z-score"}, ax=ax, linewidths=0.3, linecolor="white")
ax.set(xlabel="", ylabel="", title="observed activity (color) vs literature expectation (✓ · ✗ ○)")
savefig(fig, "gene_set_activities.literature_validation")
plt.show()

# agreement summary over the programs' sets with an expectation
s = obs.loc[obs["call"] != ""]
summary = (s.assign(confirmed=s["call"] == "✓", same_direction=s["call"].isin(["✓", "·"]), contradicted=s["call"] == "✗")
            .groupby(["time_point", "stress", "basis"])[["confirmed", "same_direction", "contradicted"]].mean().mul(100).round(0))
summary["n sets"] = s.groupby(["time_point", "stress", "basis"]).size()
display(summary)
d = s.loc[s["basis"] == "direct"]
print(f"direct expectations: {100 * (d['call'] == '✓').mean():.0f}% confirmed, "
      f"{100 * d['call'].isin(['✓', '·']).mean():.0f}% in the expected direction, {100 * (d['call'] == '✗').mean():.0f}% contradicted "
      f"(n = {len(d)} set x comparison tests)")
""")

code(r"""
from scipy import stats

# gene-level check of the marker sets: median log2FC of their genes (replicate averages, stress vs control)
# relative to all other genes (Mann-Whitney), independent of the ISMARA fit
rows = []
for sid, genes in lit_sets.items():
    for label in sc_order:
        de = de_results[label]
        inside = de.index.isin(genes)
        u = stats.mannwhitneyu(de.loc[inside, "log2FC"], de.loc[~inside, "log2FC"])
        rows.append({"set": sid, "comparison": label, "n genes": int(inside.sum()),
                     "median log2FC": de.loc[inside, "log2FC"].median(), "median log2FC, other genes": de.loc[~inside, "log2FC"].median(),
                     "p (Mann-Whitney)": u.pvalue})
marker_check = pd.DataFrame(rows)
save_tsv(marker_check, "TotalDE.literature_markers.gene_level.tsv")
mc = marker_check.pivot(index="set", columns="comparison", values="median log2FC")[sc_order]
sns.set(font_scale=0.9, style="white")
fig, ax = plt.subplots(figsize=(0.75 * len(sc_order) + 4, 3.2))
lim = np.nanmax(np.abs(mc.values))
sns.heatmap(mc, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".1f", cbar_kws={"label": "median log2FC"}, ax=ax)
ax.set(xlabel="", ylabel="", title="marker genes of the literature programs: median log2FC vs control")
savefig(fig, "literature_markers.median_log2FC")
plt.show()
""")


md(r"""
## CAMERA-like competitive test (baseline)

Per comparison, the gene-level z-scores of the differential expression are compared between the genes of a set and the other genes, with the variance inflation of an inter-gene correlation of 0.01 (Wu & Smyth 2012).
""")

code(r"""
camera_res = []
for _, row in comparisons_df.iterrows():
    z = de_results[row["label"]]["z_score"]
    c = su.camera(z, coll, inter_gene_cor=0.01)
    camera_res.append(c.assign(comparison=row["label"], kind=row["kind"]).rename_axis("set_id").reset_index())
camera_df = pd.concat(camera_res)
save_tsv(camera_df, "TotalDE.camera.tsv")
for label in comparisons_df.loc[comparisons_df["kind"] == "stress_vs_control", "label"]:
    d = camera_df.loc[(camera_df["comparison"] == label) & (camera_df["padj"] < 0.05)].sort_values("p_value")
    print(f"{label}: {len(d)} sets with padj < 0.05")
    display(d.head(8)[["set_id", "name", "size", "direction", "delta", "padj"]].round(4))
""")

code(r"""
# agreement of the two gene-set methods: activity z vs CAMERA t, stress vs control
m = act_diff_df.merge(camera_df[["set_id", "comparison", "t"]], on=["set_id", "comparison"])
m = m.loc[m["kind"] == "stress_vs_control"]
print(m.groupby("comparison").apply(lambda d: d["z_score"].corr(d["t"], method="spearman")).round(2).rename("Spearman(activity z, CAMERA t)"))
""")

md(r"""
# Expression of selected genes relative to the controls

The 10 and 60 min samples were processed as separate batches: even UT 60 vs 10 min has many differentially expressed genes. Each stress is therefore shown relative to its own control at the same time point (UT for minusD, DTT, H2O2 and DMSO; DMSO for CCCP and Tun). Points: log2 TPM of each replicate of the stress minus the replicate average of the control; markers: $\log_2$FC $= A(\text{stress}) - A(\text{control})$ with 95% credible interval $\pm 1.96\sqrt{SE_{stress}^2 + SE_{control}^2}$ (`sanity_utils.plotting.plot_gene_log2fc`).
""")

code(r"""
from sanity_utils import plotting as su_plotting

genes_of_interest_list = ["ACT1", "PGK1", "BSC1", "DCP2", "TPI1", "RLM1", "ATP11", "MET10", "RPS6A", "RPS9A", "RPS9B", "RPL19B"]

# to define the order
name_to_id = {n: g for g, n in gene_names.items()}
selected_genes = [name_to_id[n] for n in genes_of_interest_list if n in name_to_id]
print("not found:", [n for n in genes_of_interest_list if n not in name_to_id])

# stress vs control comparisons, ordered by condition, then time point; colored by condition
sc_rows = comparisons_df.loc[comparisons_df["kind"] == "stress_vs_control"].copy()
sc_rows["condition"] = sc_rows["test"].str.split(";").str[1]
sc_rows["time_point"] = sc_rows["test"].str.split(";").str[0]
sc_rows = sc_rows.assign(o=sc_rows["condition"].map(list(cat_color_dict).index)).sort_values(["o", "time_point"])
lfc_comparisons = [(f"{r.condition} vs {r.reference.split(';')[1]}, {r.time_point}", r.test, r.reference)
                   for r in sc_rows.itertuples()]
lfc_palette = {label: cat_color_dict[test.split(";")[1]] for label, test, _ in lfc_comparisons}

sns.set(font_scale=1, style="white")
fig = su_plotting.plot_gene_log2fc(
    res, avg, total_metadata, selected_genes, lfc_comparisons,
    sample_col="sample", cond_col="cohort",
    ci=0.95,
    gene_labels=gene_names.to_dict(),
    palette=lfc_palette,
    ncols=4,
    panel_size=(2.8, 4.0),
    savefig_path=os.path.join(DE_FIGS, "over_genes", "_".join(genes_of_interest_list) + ".log2FC_vs_control.png"),
)
figure_caption(os.path.join(DE_FIGS, "over_genes", "_".join(genes_of_interest_list) + ".log2FC_vs_control.png"), formats=("png",))
plt.show()
""")

nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python (sanity_utils)", "language": "python", "name": "sanity_utils"},
    "language_info": {"name": "python"},
}
nbf.write(nb, sys.argv[1])
print("written", sys.argv[1], len(cells), "cells")
