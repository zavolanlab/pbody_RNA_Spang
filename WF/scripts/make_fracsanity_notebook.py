"""Generates pbody_RNA_Spang.fracSanity.ipynb (cells below)."""
import sys
import nbformat as nbf

cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


md(r"""
# P-body recruitment of mRNAs under stress: fracSanity with Sanity_utils

Analysis of the yeast P-body pull-down (**PD**, Dcp2-HBH cCLAP) and **Total** RNA-seq (Spang lab; stresses as in Mookherjee, Rommel et al., bioRxiv 2026): untreated (**UT**), **DMSO** (vehicle), glucose starvation (**minusD**), **DTT**, **H2O2**, **CCCP**, tunicamycin (**Tun**), at **10 min** and **60 min**, 3 replicates per fraction, condition and time point.

This notebook complements `pbody_RNA_Spang.current.ipynb` (same environment file, folder layout and plotting style). It replaces the pySanity-based fracSanity section with the C++ Sanity run through [Sanity_utils](https://github.com/nimwegenLab/Sanity_utils) and leaves out the DESeq-style analysis. All heavy steps run in `WF/Snakefile-fracSanity` on SLURM; this notebook reads their outputs.

**Quantities.**
- $\alpha_g$: fraction of the mRNA of gene $g$ that is in the pulled-down (P-body) fraction, per condition and time point (fracSanity, [docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/fracSanity.md)).
- $F$: fraction of all cytoplasmic mRNA molecules in the pulled-down fraction (abundance-weighted mean of the $\alpha_g$). Because the counts are length-corrected (log2 TPM), $F$ refers to molecules, not to sequenced nucleotides.
- Differential recruitment: $\Delta\log_2\alpha_g$ of a stress relative to its control. Controls (as in the manuscript): **UT** for minusD, DTT, H2O2 (and DMSO); **DMSO** for CCCP and Tun.
- A gene counts as **more recruited** if $\Delta\log_2\alpha_g > \tau = 0.5$ and its local false sign rate (lfsr) is < 0.05; **less recruited** if $\Delta\log_2\alpha_g < -\tau$ and lfsr < 0.05.

**Gene universe.** All genes with at least one read (including rRNA and other non-coding RNAs, as in the pySanity analysis of `pbody_RNA_Spang.current.ipynb`), with length correction (log2 TPM). The per-cohort offset of $D_g$, and hence $F$, depends on this choice; the variant is set by `GENE_UNIVERSE`, `LENGTH_CORRECTION` and `SU_RUN` below.

**Kernel.** Run with the `sanity_utils` conda environment ("Python (sanity_utils)" kernel, which sets `PYTHONNOUSERSITE=1`).
""")

code(r"""
import warnings
warnings.simplefilter('ignore')

# general purpose packages
import pandas as pd
import numpy as np
import os
import sys
import json
import re
import csv
import pickle
import subprocess
import itertools
from pathlib import Path
from collections import Counter

import scipy.stats as stats
from scipy.stats import mannwhitneyu, spearmanr, pearsonr, bootstrap, fisher_exact
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

from dotenv import load_dotenv

import seaborn as sns; sns.set()
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from IPython.display import display, Image
%matplotlib inline

# for working with yaml files
import ruamel.yaml

# Sanity_utils: C++ Sanity interface, condition averages, fracSanity, gene-set analysis
import sanity_utils as su
from sanity_utils import plotting as su_plotting
from sanity_utils.genesets import load_gene_sets, set_source_table

# general-purpose plots
from zavolab_pyutils import pca_plot, plot_mean_vs_cv

# overlap plots
from matplotlib_venn import venn2, venn3


def upset_plot(contents, colors, max_subsets=30):
    # matrix layout (UpSet, Lex et al. 2014) in plain matplotlib: exclusive intersections sorted by size
    # (upsetplot 0.9 fails with pandas 3 copy-on-write)
    names = list(contents)
    counts = {}
    for g in set().union(*contents.values()):
        key = tuple(g in contents[k] for k in names)
        counts[key] = counts.get(key, 0) + 1
    subsets = sorted(counts.items(), key=lambda kv: -kv[1])[:max_subsets]
    x, y = np.arange(len(subsets)), np.arange(len(names))
    fig = plt.figure(figsize=(max(6, 0.35 * len(subsets) + 3), 3 + 0.4 * len(names)))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[3, 0.4 * len(names)], wspace=0.3, hspace=0.05)
    ax_bar, ax_mat, ax_tot = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 0])
    ax_bar.bar(x, [c for _, c in subsets], color="black")
    for xi, (_, c) in zip(x, subsets):
        ax_bar.text(xi, c, str(c), ha="center", va="bottom", fontsize=8)
    ax_bar.set(xticks=[], xlim=(-0.6, len(subsets) - 0.4), ylabel="# genes")
    for xi, (key, _) in zip(x, subsets):
        on = y[list(key)]
        ax_mat.scatter(np.full(len(y), xi), y, color="lightgrey", s=40)
        ax_mat.scatter(np.full(len(on), xi), on, color="black", s=40, zorder=3)
        ax_mat.plot([xi, xi], [on.min(), on.max()], color="black", lw=2)
    ax_mat.set(xticks=[], yticks=[], xlim=ax_bar.get_xlim(), ylim=(len(names) - 0.5, -0.5))
    ax_tot.barh(y, [len(contents[k]) for k in names], color=[colors[k] for k in names], edgecolor="black")
    ax_tot.set(ylim=ax_mat.get_ylim(), xlabel="set size")
    ax_tot.set_yticks(y, names)
    ax_tot.invert_xaxis()
    ax_tot.yaxis.tick_right()
    for ax in (ax_bar, ax_mat, ax_tot):
        sns.despine(ax=ax, left=ax is not ax_bar, bottom=ax is ax_mat)
    return fig
""")

code(r"""
print('sanity_utils', su.__version__, su.__file__)
print('Sanity binary:', su.find_sanity_binary(), '| version:', su.get_sanity_version(), '| commit:', su.get_sanity_commit())
print('numpy', np.__version__, '| pandas', pd.__version__, '| python', sys.version.split()[0])
""")

code(r"""
# 1. Load the environment variables
load_dotenv("pbody_RNA_Spang.scicore.env", override=True)

# 2. Reconstruct the subdirs dictionary
subdirs = {
    "lab_group_dir": os.getenv("LAB_GROUP_DIR"),
    "raw_sequencing_data_dir": os.getenv("RAW_SEQUENCING_DATA_DIR"),
    "main_project_dir": os.getenv("MAIN_PROJECT_DIR"),
    "wf_dir": os.getenv("WF_DIR"),
    "shared_project_dir": os.getenv("SHARED_PROJECT_DIR"),
    "temp_dir": os.getenv("TEMP_DIR"),
    "slurm_dir": os.getenv("SLURM_DIR"),
    "slurm_scripts_dir": os.getenv("SLURM_SCRIPTS_DIR"),
    "figures_dir": os.getenv("FIGURES_DIR"),
    "tables_dir": os.getenv("TABLES_DIR"),
    "fastq_dir": os.getenv("FASTQ_DIR"),
    "metadata_dir": os.getenv("METADATA_DIR"),
    "wf_runs_dir": os.getenv("WF_RUNS_DIR"),
    "configs_dir": os.getenv("CONFIGS_DIR"),
}

# 3. Reconstruct the file_paths dictionary
file_paths = {
    "genome_file": os.getenv("YEAST_GENOME_FILE"),
    "annotation_file": os.getenv("YEAST_ANNOTATION_FILE"),
}

for path in subdirs.values():
    if path:
        os.makedirs(path, exist_ok=True)

# outputs of the Sanity_utils pipeline (WF/Snakefile-fracSanity) and of this notebook
WF_version = "v2"
# pipeline variant: gene universe ('all' expressed genes or nuclear 'protein_coding') and length correction (TPM or CPM)
GENE_UNIVERSE, LENGTH_CORRECTION = "all", True
SU_RUN = "Sanity_utils_allTPM"
SU_OUT = os.path.join(subdirs["wf_runs_dir"], WF_version, "output", SU_RUN)
SU_TABLES = os.path.join(subdirs["tables_dir"], "gene_expression_analysis", "bulkRNAseq", SU_RUN)
SU_FIGS = os.path.join(subdirs["figures_dir"], "gene_expression_analysis", "bulkRNAseq", SU_RUN)
GENESETS_DIR = os.path.join(subdirs["main_project_dir"], "external_data", "genesets")
for d in (SU_TABLES, SU_FIGS):
    os.makedirs(d, exist_ok=True)

# shared helpers of the pipeline (metadata parsing, comparisons, colors)
sys.path.insert(0, os.path.join(subdirs["wf_dir"], "scripts", "sanity_utils_pipeline"))
from common import CONTROL_OF, STRESSES, CONDITION_ORDER, CONDITION_COLORS, TIME_POINTS, TIME_MARKERS, comparisons

from IPython.display import HTML

def figure_caption(path, formats=("png", "pdf")):
    # file name of a saved figure, shown with the plot; the links point to the full-size files relative to the
    # HTML report (WF/scripts/make_html_reports.py copies them to <report>/figures/, same subfolders as FIGURES_DIR)
    path = os.path.splitext(path)[0] if path.endswith((".png", ".pdf")) else path
    rel = os.path.relpath(path, subdirs["figures_dir"])
    links = " | ".join(f'<a href="figures/{rel}.{f}" target="_blank">{f}</a>' for f in formats)
    display(HTML(f'<p><b>Figure:</b> <code>{os.path.basename(path)}</code> &nbsp;({links})</p>'))

def savefig(fig, subdir, name, formats=("png", "pdf")):
    out = os.path.join(SU_FIGS, subdir)
    os.makedirs(out, exist_ok=True)
    for fmt in formats:
        fig.savefig(os.path.join(out, f"{name}.{fmt}"), bbox_inches="tight", dpi=300)
    figure_caption(os.path.join(out, name), formats)
    return os.path.join(out, name)

def save_tsv(df, subdir, name, index=False):
    out = os.path.join(SU_TABLES, subdir)
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, name), sep="\t", header=True, index=index, quoting=csv.QUOTE_NONE)
    return os.path.join(out, name)

print("Environment loaded and directories verified.")
""")

md(r"""
# Prepare .yaml config file and run WF

The read processing (UMI extraction, mapping, deduplication, featureCounts) was done with `WF/Snakefile-basic` (`v2`; see `pbody_RNA_Spang.current.ipynb`). The Sanity_utils steps run in `WF/Snakefile-fracSanity` on the featureCounts tables:

1. `prepare_inputs`: counts of all expressed genes (`gene_universe: all`; `protein_coding` keeps nuclear protein-coding genes only), gene lengths, sample metadata;
2. `sanity_fraction`: C++ Sanity (`v_method='MAP'`) on the Total and on the PD samples separately, with length correction (log2 TPM; `length_correction: false` gives log2 CPM), and ISMARA replicate averages per time point and condition (`method='reml'`);
3. `fit_cohort`: fracSanity per time point and condition ($F$, $\alpha_g$);
4. `differential`: differential recruitment, stress vs control and 60 vs 10 min;
5. `set_activities`: ISMARA-style gene-set activities on PD and Total together (GO, Reactome, KEGG, WikiPathways downloaded by `fetch_genesets`);
6. `gene_features`: sequence features per gene (UTRs, tAI, 5' end, RNAfold energies, RBP motifs);
7. `collect_tables`: final tables in `tables/gene_expression_analysis/bulkRNAseq/<SU_RUN>/`.

Run snakemake under the conda environment `snakemake` (with the SLURM executor plugin), from the vscode or a login node (the gene-set download needs internet).
""")

code(r"""
# load the config template, fill in the paths and save
yaml = ruamel.yaml.YAML()
yaml.preserve_quotes = True
with open(subdirs["wf_dir"] + "config.fracSanity.template.yaml") as f_read:
    data = yaml.load(f_read)

env_prefix = os.path.dirname(os.path.dirname(sys.executable))  # the sanity_utils conda env of this kernel
data["samples_file"] = os.path.join(subdirs["wf_runs_dir"], WF_version, "start_samples.tsv")
data["featurecounts_dir"] = os.path.join(subdirs["wf_runs_dir"], WF_version, "output",
                                         "gene_expression_quantification", "FeatureCounts_genomic_DefaultGTF") + "/"
data["gtf_file"] = file_paths["annotation_file"]
data["genome_file"] = file_paths["genome_file"]
data["gtf_CDS_and_UTRs"] = os.path.join(subdirs["wf_runs_dir"], "v1", "modified_annotation", "CDS_and_UTRs.gtf")
data["output_dir"] = SU_OUT + "/"
data["local_log"] = os.path.join(SU_OUT, "local_log") + "/"
data["cluster_log"] = os.path.join(SU_OUT, "cluster_log") + "/"
data["tables_dir"] = SU_TABLES + "/"
data["genesets_dir"] = GENESETS_DIR + "/"
data["sanity_utils_python"] = os.path.join(env_prefix, "bin", "python")
data["sanity_utils_bin"] = os.path.join(env_prefix, "bin")
data["averaging_method"] = "reml"
data["gene_universe"] = GENE_UNIVERSE
data["length_correction"] = LENGTH_CORRECTION

for dir_path in [data["output_dir"], data["local_log"], data["cluster_log"]]:
    os.makedirs(dir_path, exist_ok=True)

config_path = os.path.join(subdirs["wf_runs_dir"], WF_version, f"modified_config_fracSanity.{SU_RUN.replace('Sanity_utils_', '')}.yaml")
with open(config_path, "w") as f_write:
    yaml.dump(data, f_write)

command = (
    "snakemake --snakefile " + subdirs["wf_dir"] + "Snakefile-fracSanity"
    + " --configfile " + config_path
    + " --executor slurm --profile " + subdirs["wf_dir"] + "profile"
    + " --printshellcmds --rerun-incomplete --keep-going --nolock -np"
)
print(command)
""")

md("# Define metadata labels and colors uniformly")

code(r"""
metadata_df = pd.read_csv(os.path.join(SU_OUT, "inputs", "metadata.tsv"), sep="\t")
genes_df = pd.read_csv(os.path.join(SU_OUT, "inputs", "genes.tsv"), sep="\t", index_col=0)
gene_names = genes_df["gene_name"].fillna(pd.Series(genes_df.index, index=genes_df.index))

cat_color_dict = dict(CONDITION_COLORS)   # condition -> color (as in pbody_RNA_Spang.current.ipynb)
mrk_marker_dict = dict(TIME_MARKERS)      # time point -> marker
comparisons_df = comparisons()
stress_comparisons = comparisons_df.loc[comparisons_df["kind"] == "stress_vs_control"].reset_index(drop=True)

metadata_df["t"] = 1
display(metadata_df.groupby(["fraction", "time_point", "condition"]).agg({"t": "sum"}).unstack())
display(comparisons_df)
""")

md(r"""
# mRNA Expression
## Analysis of read mapping stats
""")

code(r"""
dir_path = os.path.join(subdirs["wf_runs_dir"], WF_version, "output", "mapping_stats")

def parse_mapping_stats(path):
    # custom format: '>' header lines followed by one integer each
    res, key = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:]
            elif line and key is not None:
                res[key] = int(float(line))
                key = None
    return res

rows = []
for s in metadata_df["sample"]:
    f = os.path.join(dir_path, s + ".mapping_stats.txt")
    if os.path.exists(f) and os.path.exists(os.path.join(dir_path, s + ".success")):
        rows.append({"sample": s, **parse_mapping_stats(f)})
mapping_stats_df = pd.merge(metadata_df, pd.DataFrame(rows), how="inner", on="sample")
mapping_stats_df["number of MM reads; genomic mapping"] = (
    mapping_stats_df["total number of mapped reads; genomic mapping"] - mapping_stats_df["number of UM reads; genomic mapping"])
mapping_stats_df["number of MM reads; UMI dedup"] = (
    mapping_stats_df["total number of mapped reads; UMI dedup"] - mapping_stats_df["number of UM reads; UMI dedup"])
for elem in [c for c in mapping_stats_df.columns if c.startswith(("total number", "number of"))]:
    mapping_stats_df[elem + "; mln"] = np.round(mapping_stats_df[elem] / 10**6, 2)
mapping_stats_df["perc_of_UM; genomic mapping"] = (mapping_stats_df["number of UM reads; genomic mapping"]
    / mapping_stats_df["total number of mapped reads; genomic mapping"] * 100)
mapping_stats_df["perc_of_UM; UMI dedup"] = (mapping_stats_df["number of UM reads; UMI dedup"]
    / mapping_stats_df["total number of mapped reads; UMI dedup"] * 100)
mapping_stats_df = mapping_stats_df.sort_values(["time_point", "fraction", "condition", "replicate"])
print(len(mapping_stats_df), "samples with mapping stats")
""")

code(r"""
features = [
    "total number of mapped reads; genomic mapping; mln",
    "number of UM reads; genomic mapping; mln",
    "total number of mapped reads; UMI dedup; mln",
    "number of UM reads; UMI dedup; mln",
]
palette = ["black", "green", "grey", "royalblue"]

sns.set(font_scale=1)
sns.set_style("white")
fig, axes = plt.subplots(1, 1, sharey=False, sharex=True, figsize=(4, 15))
for i, feature in enumerate(features):
    ax = sns.pointplot(data=mapping_stats_df, x=feature, y="sample_label_within_experiment",
                       color=palette[i], label=feature, linestyles="none")
ax.legend(bbox_to_anchor=(1.05, 1.0), loc=2, borderaxespad=0.0, title="", ncols=1)
ax.tick_params(left=True, bottom=True)
ax.set(ylabel="", xlabel="# reads, mln")
savefig(fig, "QC", "read_mapping.read_numbers")
""")

code(r"""
# composition of the counted reads (deduplicated, uniquely mapped, assigned by featureCounts) by gene class,
# from the counts of all genes; independent of the gene universe used for Sanity (GENE_UNIVERSE)
counts_all = pd.read_csv(os.path.join(SU_OUT, "inputs", "counts.tsv"), sep="\t", index_col=0)
genes_all = pd.read_csv(os.path.join(SU_OUT, "inputs", "genes.tsv"), sep="\t", index_col=0).reindex(counts_all.index)
gene_class = np.select(
    [genes_all["chrom"] == "Mito", genes_all["gene_biotype"] == "protein_coding", genes_all["gene_biotype"] == "rRNA",
     genes_all["gene_biotype"] == "tRNA", genes_all["gene_biotype"].isin(["snoRNA", "snRNA", "ncRNA"])],
    ["mitochondrial", "nuclear protein-coding", "rRNA", "tRNA", "snoRNA, snRNA, ncRNA"], default="other")
composition = 100 * counts_all.groupby(gene_class).sum().T.div(counts_all.sum(), axis=0)
composition = pd.merge(metadata_df, composition.rename_axis("sample").reset_index(), on="sample")
classes = ["nuclear protein-coding", "rRNA", "tRNA", "snoRNA, snRNA, ncRNA", "mitochondrial", "other"]
composition = composition.sort_values(["time_point", "fraction", "condition", "replicate"])
if GENE_UNIVERSE != "all":
    print(f"note: counts.tsv holds only the '{GENE_UNIVERSE}' genes; the composition is relative to them")

sns.set(font_scale=1)
sns.set_style("white")
fig, axes = plt.subplots(1, 2, figsize=(8, 15), sharey=True)
for ax, feature in zip(axes, ["nuclear protein-coding", "rRNA"]):
    sns.pointplot(data=composition, x=feature, y="sample_label_within_experiment",
                  hue="condition", palette=cat_color_dict, linestyles="none", ax=ax)
    ax.set(xlabel=f"% of counted reads: {feature}", ylabel="", xlim=(0, 100))
    ax.tick_params(left=True, bottom=True)
axes[0].legend_.remove()
axes[1].legend(bbox_to_anchor=(1.05, 1.0), loc=2, borderaxespad=0.0, title="")
savefig(fig, "QC", "fraction_of_counts_per_gene_class")
save_tsv(composition[["sample", "sample_label_within_experiment"] + [c for c in classes if c in composition]],
         "QC", "fraction_of_counts_per_gene_class.tsv")
display(composition.groupby(["fraction", "time_point"])[[c for c in classes if c in composition]].median().round(1))
""")
code(r"""
# read composition by gene biotype of the annotation (Ensembl GTF 'gene_biotype'; mitochondrial genes by their biotype)
biotype_composition = 100 * counts_all.groupby(genes_all["gene_biotype"].to_numpy()).sum().T.div(counts_all.sum(), axis=0)
biotype_order = ["protein_coding", "transposable_element", "rRNA", "snoRNA", "snRNA", "ncRNA", "tRNA", "pseudogene"]
biotype_order = [b for b in biotype_order if b in biotype_composition] + [b for b in biotype_composition if b not in biotype_order]
biotype_colors = dict(zip(biotype_order, ["lightgrey", "crimson", "royalblue", "mediumseagreen", "gold", "darkorange",
                                          "plum", "sienna"] + ["black"] * len(biotype_order)))
lab = metadata_df.set_index("sample")["sample_label_within_experiment"]
plot_df = biotype_composition.loc[composition["sample"]][biotype_order]   # same sample order as above
plot_df.index = plot_df.index.map(lab)

sns.set(font_scale=1, style="white")
fig, ax = plt.subplots(1, 1, figsize=(6, 15))
left = np.zeros(len(plot_df))
for b in biotype_order:
    ax.barh(np.arange(len(plot_df)), plot_df[b], left=left, color=biotype_colors[b], edgecolor="black", lw=0.3, label=b)
    left += plot_df[b].to_numpy()
ax.set_yticks(np.arange(len(plot_df)), plot_df.index, fontsize=8)
ax.set_ylim(len(plot_df) - 0.5, -0.5)
ax.set(xlabel="% of counted reads", xlim=(0, 100))
ax.legend(bbox_to_anchor=(1.02, 1.0), loc=2, borderaxespad=0.0, title="gene biotype")
savefig(fig, "QC", "fraction_of_counts_per_biotype")
save_tsv(biotype_composition.rename_axis("sample").reset_index(), "QC", "fraction_of_counts_per_biotype.tsv")
display(pd.merge(metadata_df[["sample", "fraction", "time_point"]], biotype_composition.rename_axis("sample").reset_index(), on="sample")
        .groupby(["fraction", "time_point"])[biotype_order].median().round(1))
""")



md(r"""
## Sanity normalization (C++ Sanity via Sanity_utils)

Sanity was run separately on the Total and on the PD samples (the pull-down adds its own variability), with `v_method='MAP'` and gene-length correction: `log2_expression` is log2 TPM (transcripts per million transcripts of the universe). Replicates were averaged per time point and condition with ISMARA's replicate averaging, `method='reml'` (recommended for 3 replicates, [condition_averages.md](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/condition_averages.md)).
""")

code(r"""
def load_pickle(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)

sanity = {f: load_pickle(os.path.join(SU_OUT, "sanity", f, "sanity.pkl")) for f in ("Total", "PD")}
for f, d in sanity.items():
    res = d["result"]
    print(f, res.log2_expression.shape, "| units:", res.params.get("units"), "| v_method:", res.params.get("v_method"),
          "| Sanity commit:", res.params.get("sanity_commit"))
""")

code(r"""
for FRACTION in ("Total", "PD"):
    res = sanity[FRACTION]["result"]
    meta_f = sanity[FRACTION]["metadata"].copy()
    fig = su_plotting.plot_variance_vs_expression(res)
    savefig(fig, "QC", f"variance_vs_expr.{FRACTION}")
    plt.show()

    # CV vs mean of the relative expression (natural scale, fractions of the transcriptome)
    natScale = 2 ** (res.log2_expression - np.log2(1e6))
    out = os.path.join(SU_FIGS, "QC")
    plot_mean_vs_cv(natScale, meta_f, os.path.join(out, f"cv_vs_mean_plot.per_cohort.{FRACTION}.png"),
                    sample_col="sample", cond_col="cohort")
    figure_caption(os.path.join(out, f"cv_vs_mean_plot.per_cohort.{FRACTION}.png"), formats=("png",))
    meta_f["all"] = "all"
    plot_mean_vs_cv(natScale, meta_f, os.path.join(out, f"cv_vs_mean_plot.all_cohorts.{FRACTION}.png"),
                    sample_col="sample", cond_col="all")
    figure_caption(os.path.join(out, f"cv_vs_mean_plot.all_cohorts.{FRACTION}.png"), formats=("png",))
""")

code(r"""
# PCA on Sanity log2 TPM, genes with at least one read in all samples
counts_df = pd.read_csv(os.path.join(SU_OUT, "inputs", "counts.tsv"), sep="\t", index_col=0)
expressed_everywhere_genes = list(counts_df.index[(counts_df > 0).all(axis=1)])
print(len(expressed_everywhere_genes), "genes with reads in all samples")

for FRACTION in ("Total", "PD"):
    res = sanity[FRACTION]["result"]
    sub_meta = sanity[FRACTION]["metadata"].copy().reset_index(drop=True)
    genes = [g for g in expressed_everywhere_genes if g in res.log2_expression.index]
    PCA_UMAP_subset = res.log2_expression.loc[genes, sub_meta["sample"]].reset_index(drop=True)
    savefig_path = os.path.join(SU_FIGS, "QC", f"PCA_gene_expression.{FRACTION}.png")
    pca_plot(
        PCA_UMAP_subset, list(sub_meta["sample"]), sub_meta, "condition", savefig_path,
        sns_color_palette=cat_color_dict, hue_order=[c for c in CONDITION_ORDER],
        style_feature="time_point", markers=mrk_marker_dict,
        calculate_permanova_R2=True, permanova_feature="cohort", permanova_R2_ajusted=False,
        add_text_labels_for_samples=True, text_label_column_for_samples="replicate",
        s_param=40, figsize=(4.2, 4.2), legend_title=FRACTION + " fraction\nconditions:",
        add_2D_KDE_countours=False,
    )
    figure_caption(savefig_path, formats=("png",))
    PCA_UMAP_subset.to_csv(str(Path(savefig_path).with_suffix(".tsv")), sep="\t", header=True, index=None)
    plt.show()
""")

code(r"""
# normalized expression tables (log2 TPM and TPM) with gene annotation, both fractions
annot = genes_df[["gene_name", "gene_biotype", "chrom", "length"]].copy()
log2_tpm = pd.concat([sanity["Total"]["result"].log2_expression, sanity["PD"]["result"].log2_expression], axis=1)
save_tsv(annot.join(log2_tpm, how="right").reset_index().rename(columns={"index": "gene_id"}),
         "normalized_expression", "Sanity_log2TPM.tsv")
save_tsv(annot.join(2 ** log2_tpm, how="right").reset_index().rename(columns={"index": "gene_id"}),
         "normalized_expression", "Sanity_TPM.tsv")
save_tsv(metadata_df.drop(columns=["t"]), "normalized_expression", "Sanity.metadata.tsv")
for FRACTION in ("Total", "PD"):
    avg = sanity[FRACTION]["averages"]
    save_tsv(avg.means.add_prefix("mean;").join(avg.errors.add_prefix("error;")).reset_index()
             .rename(columns={"index": "gene_id"}), "normalized_expression", f"Sanity_condition_averages.{FRACTION}.tsv")
""")

md(r"""
# fracSanity: global fraction $F$ of mRNA in P-bodies

For each time point and condition, fracSanity fits $D_g = \langle A^{PD}_g\rangle - \langle A^{Total}_g\rangle = \log_2\alpha_g - \log_2 F + \varepsilon_g$ with $\alpha_g \sim$ Beta$(a, b)$ and a contaminated-Gaussian error model, by profile likelihood over $\log_2 F$ (uniform prior on $F$). $F$ is reported as the posterior median with a 95% credible interval.
""")

code(r"""
globalF_df = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_global_F.tsv"), sep="\t")
deltaF_df = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_delta_global_F.tsv"), sep="\t")
display(globalF_df[["cohort", "F", "F_CI_lower", "F_CI_upper", "F_identifiable", "a", "b", "pi_outlier", "tau", "ks_stat", "n_genes"]].round(4))
display(deltaF_df.round(3))
""")

code(r"""
# fit diagnostics per cohort (profile likelihood, PIT histogram, alpha_g vs fitted Beta)
for cohort in sorted(globalF_df["cohort"]):
    fit = load_pickle(os.path.join(SU_OUT, "fracSanity", "fits", cohort.replace(";", "__") + ".pkl"))
    fig = su_plotting.plot_frac_sanity_fit(fit, condition_name=cohort)
    savefig(fig, os.path.join("fracSanity", "fit_diagnostics"), "fit." + cohort.replace(";", "_"), formats=("png",))
    plt.close(fig)
""")

code(r"""
# F per condition and time point, with 95% credible intervals
sns.set(font_scale=1, style="ticks")
fig, ax = plt.subplots(1, 1, figsize=(6, 3.5))
data = globalF_df.copy()
data["condition"] = pd.Categorical(data["condition"], CONDITION_ORDER, ordered=True)
data = data.sort_values(["condition", "time_point"])
x_labels = []
for i, (_, row) in enumerate(data.iterrows()):
    x = list(CONDITION_ORDER).index(row["condition"]) + (-0.18 if row["time_point"] == "10min" else 0.18)
    ax.errorbar(x, row["F"], yerr=[[row["F"] - row["F_CI_lower"]], [row["F_CI_upper"] - row["F"]]],
                fmt=mrk_marker_dict[row["time_point"]], color=cat_color_dict[row["condition"]],
                mec="black", ms=8, ecolor="black", capsize=2)
ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels(CONDITION_ORDER)
ax.set(ylabel=r"Global fraction $\hat{F}$ in PD", xlabel="")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0], [0], marker=m, color="w", markerfacecolor="grey", markeredgecolor="black", ms=8, label=t)
                   for t, m in mrk_marker_dict.items()], loc="upper left", bbox_to_anchor=(1.02, 1.0), title="time point")
savefig(fig, "fracSanity", "globalF_per_condition")
""")

md(r"""
## Ground truth: P-body imaging (Dcp2-GFP)

Imaging data of Mookherjee, Rommel et al. (bioRxiv 2026): P-body brightness = mean gray value of the P-body normalized to the cytoplasm after subtracting autofluorescence, one value per P-body (`materials/pbody_intensity.tsv`, stresses only), and the number of P-bodies per cell (`materials/pbody_number_*.tsv`, all conditions). Means with 95% bootstrap intervals; the means reproduce the values printed in Fig. 1D/F of the preprint.

If the pull-down measures the mRNA content of P-bodies, $F$ should increase with P-body brightness across conditions (right panel), and, within each stress, the change of $F$ from 10 to 60 min should have the same sign as the change of the mean brightness.
""")

code(r"""
agg_func = np.mean

def boot_ci(x, n_resamples=10000, seed=0):
    x = x.dropna().to_numpy()
    ci = bootstrap((x,), agg_func, n_resamples=n_resamples, confidence_level=0.95,
                   method="percentile", rng=np.random.default_rng(seed)).confidence_interval
    return pd.Series({"y_CI_lower": ci.low, "y_CI_upper": ci.high})

pbody_intensity_df = pd.read_csv(os.path.join(subdirs["main_project_dir"], "materials", "pbody_intensity.tsv"), sep="\t")
grouped = pbody_intensity_df.groupby(["time_point", "condition"])["pb/cyto"]
pbody_intensity_gr_df = pd.merge(grouped.agg(agg_func).reset_index(), grouped.apply(boot_ci).unstack().reset_index(),
                                 on=["time_point", "condition"])
pbody_intensity_gr_df["n_pbodies"] = grouped.size().to_numpy()

cols10 = ["sub_image", "UT", "DTT", "minusD", "H2O2", "DMSO", "Tun", "CCCP"]
cols60 = ["sub_image", "UT", "DMSO", "DTT", "minusD", "H2O2", "Tun", "CCCP"]
n10 = pd.read_csv(os.path.join(subdirs["main_project_dir"], "materials", "pbody_number_10min.tsv"), sep="\t")
n10.columns = cols10
n10["time_point"] = "10min"
n60 = pd.read_csv(os.path.join(subdirs["main_project_dir"], "materials", "pbody_number_60min.tsv"), sep="\t")
n60.columns = cols60
n60["time_point"] = "60min"
pbody_number_melted_df = pd.melt(pd.concat([n10[cols10 + ["time_point"]], n60[cols10 + ["time_point"]]]),
                                 id_vars=["sub_image", "time_point"], value_vars=CONDITION_ORDER,
                                 var_name="condition", value_name="pbody_abundance").dropna()
grouped = pbody_number_melted_df.groupby(["time_point", "condition"])["pbody_abundance"]
pbody_number_gr_df = pd.merge(grouped.agg(agg_func).reset_index(), grouped.apply(boot_ci).unstack().reset_index(),
                              on=["time_point", "condition"])
display(pbody_intensity_gr_df.round(3))
display(pbody_number_gr_df.round(3))
""")

code(r"""
def p_to_stars(p):
    for thr, s in [(1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")]:
        if p <= thr:
            return s
    return "ns"


def violin_with_stats(ax0, df_plot, y_feature, ylim, reference, ylabel, round_medians=2):
    # violins per time point and condition, median labels, Mann-Whitney vs reference within time point
    x_feature = "time_point"
    x_order = list(mrk_marker_dict)
    hue_order = [c for c in cat_color_dict if c in set(df_plot["condition"])]
    sns.violinplot(ax=ax0, data=df_plot, x=x_feature, y=y_feature, order=x_order, hue="condition",
                   hue_order=hue_order, palette=cat_color_dict, bw_adjust=0.2, cut=0, inner="quart")
    ax0.legend_.remove()
    n_hue, dodge_width = len(hue_order), 0.8

    def x_pos(x_cat, hue_cat):
        i, j = x_order.index(x_cat), hue_order.index(hue_cat)
        return i + dodge_width * ((j + 0.5) / n_hue - 0.5)

    grouped = df_plot.dropna(subset=[y_feature]).groupby([x_feature, "condition"])[y_feature]
    medians, maxima = np.round(grouped.median(), round_medians), grouped.max()
    label_offset = 0.02 * (ylim[1] - ylim[0])
    for (tp, cond), med in medians.items():
        ax0.text(x_pos(tp, cond), maxima[(tp, cond)] + label_offset, f"{med:g}", ha="center", va="bottom",
                 fontsize=8, color="dimgray")
    if reference in hue_order:
        bar_step, bar_gap, bar_tick = (0.07 * (ylim[1] - ylim[0]), 0.08 * (ylim[1] - ylim[0]),
                                       0.015 * (ylim[1] - ylim[0]))
        y_top = ylim[1]
        for tp in x_order:
            ref_vals = df_plot.loc[(df_plot[x_feature] == tp) & (df_plot["condition"] == reference), y_feature]
            j_ref, prev_y = hue_order.index(reference), -np.inf
            for j, cond in enumerate(hue_order):
                if cond == reference:
                    continue
                vals = df_plot.loc[(df_plot[x_feature] == tp) & (df_plot["condition"] == cond), y_feature]
                if len(ref_vals) == 0 or len(vals) == 0:
                    continue
                stars = p_to_stars(mannwhitneyu(ref_vals, vals, alternative="two-sided").pvalue)
                if stars == "ns":
                    continue
                spanned = [maxima.get((tp, c), -np.inf) for c in hue_order[min(j, j_ref):max(j, j_ref) + 1]]
                y = max(max(spanned) + bar_gap, prev_y + bar_step)
                x1, x2 = x_pos(tp, reference), x_pos(tp, cond)
                ax0.plot([x1, x1, x2, x2], [y - bar_tick, y, y, y - bar_tick], lw=1, c="black")
                ax0.text((x1 + x2) / 2, y, stars, ha="center", va="bottom", fontsize=9)
                prev_y = y
                y_top = max(y_top, y + bar_step)
        ax0.set_ylim(ylim[0], max(ylim[1], y_top))
    else:
        ax0.set_ylim(*ylim)
    ax0.set(ylabel=ylabel, xlabel="Time point")


def F_scatter(ax1, data, y_feature, ylim, ylabel=""):
    sns.scatterplot(ax=ax1, data=data, x="F", y=y_feature, hue="condition",
                    hue_order=[c for c in cat_color_dict if c in set(data["condition"])], palette=cat_color_dict,
                    style="time_point", style_order=list(mrk_marker_dict), markers=mrk_marker_dict, s=80,
                    edgecolor="black", zorder=3)
    ax1.errorbar(data["F"], data[y_feature],
                 xerr=[data["F"] - data["F_CI_lower"], data["F_CI_upper"] - data["F"]],
                 yerr=[data[y_feature] - data["y_CI_lower"], data["y_CI_upper"] - data[y_feature]],
                 fmt="none", ecolor="black", capsize=0.1, alpha=0.7, zorder=2)
    ax1.set(ylim=ylim, ylabel=ylabel, xlabel=r"Estimated Global Fraction ($\hat{F}$)")
    corr_df = data[["F", y_feature]].dropna()
    rho, p_rho = stats.spearmanr(corr_df["F"], corr_df[y_feature])
    r, p_r = stats.pearsonr(corr_df["F"], corr_df[y_feature])
    ax1.set_title(f"Spearman ρ = {rho:.2f} (p = {p_rho:.2g})\nPearson r = {r:.2f} (p = {p_r:.2g}), n = {len(corr_df)}",
                  fontsize=10)
    ax1.legend(loc="upper left", bbox_to_anchor=(1.05, 1.05), frameon=True, facecolor="white", title="")
    return rho, p_rho, r, p_r


data_int = pd.merge(globalF_df, pbody_intensity_gr_df, on=["time_point", "condition"], how="inner")
sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1.5, 1]})
violin_with_stats(axes[0], pbody_intensity_df, "pb/cyto", (0, 30), reference="UT",
                  ylabel="normalized mean gray value\n(PB / cytoplasm)")
corr_intensity = F_scatter(axes[1], data_int, "pb/cyto", (0, 30))
savefig(fig, "fracSanity", "globalF_vs_pBodyIntensity.Sanity_utils")
save_tsv(data_int, "fracSanity", "globalF_vs_pBodyIntensity.Sanity_utils.tsv")
""")

code(r"""
data_num = pd.merge(globalF_df, pbody_number_gr_df, on=["time_point", "condition"], how="inner")
sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1.5, 1]})
violin_with_stats(axes[0], pbody_number_melted_df, "pbody_abundance", (0, 10), reference="UT",
                  ylabel="# PBs / cell", round_medians=1)
corr_number = F_scatter(axes[1], data_num, "pbody_abundance", (0, 5), ylabel="mean # PBs / cell")
savefig(fig, "fracSanity", "globalF_vs_pBodyAbundance.Sanity_utils")
save_tsv(data_num, "fracSanity", "globalF_vs_pBodyAbundance.Sanity_utils.tsv")
""")

code(r"""
# Within each stress: change from 10 to 60 min of F (fracSanity posterior) and of the mean P-body brightness,
# both as log2 fold changes (brightness: log2 of the ratio of the means, bootstrap over P-bodies, the two time points
# resampled independently), so that the two axes are directly comparable
rng = np.random.default_rng(1)
rows = []
for cond in STRESSES:
    x10 = pbody_intensity_df.loc[(pbody_intensity_df["condition"] == cond) & (pbody_intensity_df["time_point"] == "10min"), "pb/cyto"].dropna().to_numpy()
    x60 = pbody_intensity_df.loc[(pbody_intensity_df["condition"] == cond) & (pbody_intensity_df["time_point"] == "60min"), "pb/cyto"].dropna().to_numpy()
    boots = np.array([np.log2(rng.choice(x60, len(x60)).mean() / rng.choice(x10, len(x10)).mean()) for _ in range(10000)])
    d = deltaF_df.loc[deltaF_df["comparison"] == f"{cond};60min_vs_10min"].iloc[0]
    rows.append({
        "condition": cond,
        "delta_mean_intensity": x60.mean() - x10.mean(),
        "log2FC_intensity": np.log2(x60.mean() / x10.mean()),
        "log2FC_intensity_CI_lower": np.quantile(boots, 0.025), "log2FC_intensity_CI_upper": np.quantile(boots, 0.975),
        "prob_intensity_increase": float(np.mean(boots > 0)),
        "mannwhitney_p": mannwhitneyu(x60, x10).pvalue,
        "F_10min": d["F_reference"], "F_60min": d["F_test"], "delta_F": d["F_test"] - d["F_reference"],
        "delta_log2_F": d["delta_log2_F"], "delta_log2_F_CI_lower": d["delta_log2_F_CI_lower"],
        "delta_log2_F_CI_upper": d["delta_log2_F_CI_upper"], "prob_F_increase": d["prob_delta_log2_F_positive"],
    })
consistency_df = pd.DataFrame(rows)
consistency_df["sign_intensity"] = np.sign(consistency_df["log2FC_intensity"])
consistency_df["sign_F"] = np.sign(consistency_df["delta_log2_F"])
consistency_df["signs_agree"] = consistency_df["sign_intensity"] == consistency_df["sign_F"]
# a sign is 'credible' if the 95% interval excludes 0
consistency_df["intensity_change_credible"] = (consistency_df["log2FC_intensity_CI_lower"] > 0) | (consistency_df["log2FC_intensity_CI_upper"] < 0)
consistency_df["F_change_credible"] = (consistency_df["delta_log2_F_CI_lower"] > 0) | (consistency_df["delta_log2_F_CI_upper"] < 0)
save_tsv(consistency_df, "fracSanity", "deltaF_vs_deltaIntensity_60_vs_10min.tsv")
display(consistency_df.round(3))
print(f"signs agree in {consistency_df['signs_agree'].sum()} of {len(consistency_df)} stresses")

sns.set(font_scale=1, style="ticks")
fig, ax = plt.subplots(1, 1, figsize=(5.2, 5.0))
for _, r in consistency_df.iterrows():
    ax.errorbar(r["log2FC_intensity"], r["delta_log2_F"],
                xerr=[[r["log2FC_intensity"] - r["log2FC_intensity_CI_lower"]], [r["log2FC_intensity_CI_upper"] - r["log2FC_intensity"]]],
                yerr=[[r["delta_log2_F"] - r["delta_log2_F_CI_lower"]], [r["delta_log2_F_CI_upper"] - r["delta_log2_F"]]],
                fmt="o", ms=9, color=cat_color_dict[r["condition"]], mec="black", ecolor="black", capsize=2, label=r["condition"])
lim = np.nanmax(np.abs(consistency_df[["log2FC_intensity_CI_lower", "log2FC_intensity_CI_upper",
                                       "delta_log2_F_CI_lower", "delta_log2_F_CI_upper"]].to_numpy())) * 1.05
ax.fill_between([0, lim], 0, lim, color="green", alpha=0.06, lw=0)
ax.fill_between([-lim, 0], -lim, 0, color="green", alpha=0.06, lw=0)
ax.plot([-lim, lim], [-lim, lim], color="grey", ls="--", lw=1, label="y = x")
ax.axhline(0, color="grey", lw=0.8)
ax.axvline(0, color="grey", lw=0.8)
ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal",
       xlabel="log2FC of mean PB brightness (60 vs 10 min)", ylabel=r"$\Delta\log_2\hat{F}$ (60 vs 10 min)")
ax.set_title(f"same sign in {consistency_df['signs_agree'].sum()} / {len(consistency_df)} stresses", fontsize=10)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), title="")
savefig(fig, "fracSanity", "deltaF_vs_deltaIntensity_60_vs_10min")
""")


md(r"""
# Recruitment of individual mRNAs relative to the controls

Each stress is compared with its own control at the same time point (UT for minusD, DTT, H2O2 and DMSO; DMSO for CCCP and Tun), as in `pbody_RNA_Spang.TotalDE.ipynb`. Per gene and cohort, $D_g = \langle A^{PD}_g\rangle - \langle A^{Total}_g\rangle = \log_2\alpha_g - \log_2 F + \varepsilon_g$ (Sanity replicate averages). Because $F$ is weakly identified (see the caveats), the recruitment of a gene is measured relative to the other genes of the same cohort: $D_g$ minus its median over genes. Its change between stress and control, $\Delta D_g$, is then independent of $F$ and equals $\Delta\log_2\alpha_g$ up to a constant per comparison ($\Delta\log_2 F$). It is tested with the threshold test of the differential expression ([docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/differential_analysis.md)): **more / less recruited** if $|\Delta D_g| > \tau = 0.5$ at `padj` < 0.05. The fracSanity estimate $\Delta\log_2\alpha_g$ (with the uncertainty of $F$) is reported alongside.
""")

code(r"""
from adjustText import adjust_text
from types import SimpleNamespace

TAU = 0.5
tot_avg, pd_avg = sanity["Total"]["averages"], sanity["PD"]["averages"]
genes_rec = tot_avg.means.index.intersection(pd_avg.means.index)
D_raw = pd_avg.means.loc[genes_rec] - tot_avg.means.loc[genes_rec]
D_avg = D_raw - D_raw.median()                       # relative recruitment per cohort, independent of F
D_se = np.sqrt(pd_avg.errors.loc[genes_rec] ** 2 + tot_avg.errors.loc[genes_rec] ** 2)
dalpha = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_differential_recruitment.stress_vs_control.tsv"), sep="\t")

sc_comparisons = comparisons_df.loc[comparisons_df["kind"] == "stress_vs_control"]
rec_results, rec_summary = {}, []
for _, row in sc_comparisons.iterrows():
    t, r = row["test"], row["reference"]
    d = su.differential.z_test(D_avg[t], D_se[t], D_avg[r], D_se[r], index=genes_rec, lfc_threshold=TAU)
    d["A_Total"] = (tot_avg.means.loc[genes_rec, t] + tot_avg.means.loc[genes_rec, r]) / 2
    d["call"] = np.where(d["padj"] < 0.05, np.where(d["log2FC"] > 0, "more recruited", "less recruited"), "n.s.")
    a = dalpha.loc[dalpha["comparison"] == row["label"]].set_index("gene_id")
    d["delta_log2_alpha"], d["lfsr_alpha"] = a["delta_log2_alpha"].reindex(genes_rec), a["lfsr"].reindex(genes_rec)
    rec_results[row["label"]] = d
    rec_summary.append({"comparison": row["label"], "test": t, "reference": r,
                        "n more recruited": int((d["call"] == "more recruited").sum()),
                        "n less recruited": int((d["call"] == "less recruited").sum()),
                        "median delta_log2_alpha (= delta log2 F shift)": d["delta_log2_alpha"].median()})
rec_summary_df = pd.DataFrame(rec_summary)
display(rec_summary_df.round(2))
rec_long = pd.concat([d.assign(comparison=k) for k, d in rec_results.items()]).rename_axis("gene_id").reset_index()
rec_long.insert(1, "gene_name", rec_long["gene_id"].map(gene_names))
rec_long = rec_long.rename(columns={"log2FC": "delta_D_centered"})
save_tsv(rec_long, "fracSanity", f"recruitment_vs_control.delta_D_centered.tau{TAU}.tsv")
for d in rec_results.values():
    d.rename(columns={"log2FC": "delta_D_centered"}, inplace=True)

# number of more / less recruited genes per comparison
s = rec_summary_df.assign(condition=rec_summary_df["test"].str.split(";").str[1], time_point=rec_summary_df["test"].str.split(";").str[0])
sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
for ax, tp in zip(axes, mrk_marker_dict):
    d = s.loc[s["time_point"] == tp].set_index("condition")
    d = d.reindex([c for c in cat_color_dict if c in d.index])
    x = np.arange(len(d))
    ax.bar(x, d["n more recruited"], color=[cat_color_dict[c] for c in d.index], edgecolor="black")
    ax.bar(x, -d["n less recruited"], color=[cat_color_dict[c] for c in d.index], edgecolor="black", hatch="//")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, [f"{c}\nvs {d.loc[c, 'reference'].split(';')[1]}" for c in d.index], fontsize=9)
    ax.set_title(tp)
axes[0].set_ylabel(f"# genes more / less recruited\n(padj < 0.05, τ = {TAU})")
sns.despine()
savefig(fig, "recruitment", f"n_recruited_genes.tau{TAU}")
plt.show()
""")

md(r"""
## Recruitment change vs expression

$A$: mean log2 TPM in Total of the two conditions; genes with `padj` < 0.05 in color. Labelled: genes discussed by the Spang lab (Wang et al. eLife 2018; thesis slides).
""")

code(r"""
genes_of_interest_list = ["ACT1", "PGK1", "BSC1", "DCP2", "TPI1", "RLM1", "ATP11", "MET10", "RPS6A", "RPS9A", "RPS9B", "RPL19B"]
name_to_id = {n: g for g, n in gene_names.items()}
goi_ids = [name_to_id[n] for n in genes_of_interest_list if n in name_to_id]
call_colors = {"more recruited": "firebrick", "less recruited": "steelblue", "n.s.": "lightgrey"}


def ma_plot(ax, d, title, title_color):
    for call in ("n.s.", "less recruited", "more recruited"):
        x = d.loc[d["call"] == call]
        ax.scatter(x["A_Total"], x["delta_D_centered"], s=3 if call == "n.s." else 5, color=call_colors[call], lw=0,
                   rasterized=True, label=f"{call} ({len(x)})")
    for y in (0, TAU, -TAU):
        ax.axhline(y, color="black", lw=0.8 if y == 0 else 0.5, ls="-" if y == 0 else "--")
    g = d.loc[d.index.intersection(goi_ids)]
    ax.scatter(g["A_Total"], g["delta_D_centered"], s=14, facecolor="none", edgecolor="black", lw=0.8)
    ax.set_title(title, fontsize=10, color="black", bbox=dict(facecolor=title_color, alpha=0.35, lw=0))
    ax.legend(loc="lower right", fontsize=7, markerscale=2, frameon=False, handletextpad=0.1)
    return [ax.text(r["A_Total"], r["delta_D_centered"], gene_names[gid], fontsize=6) for gid, r in g.iterrows()]


sns.set(font_scale=1, style="ticks")
nrows, ncols = 2, 6
fig, axes = plt.subplots(nrows, ncols, figsize=(3 * ncols, 2.9 * nrows), sharex=True, sharey=True, squeeze=False)
labels = []
for ax, (_, row) in zip(axes.flat, sc_comparisons.iterrows()):
    cond = row["test"].split(";")[1]
    labels.append((ax, ma_plot(ax, rec_results[row["label"]], row["label"].replace("_vs_", " vs "), cat_color_dict[cond])))
for ax in axes[-1]:
    ax.set_xlabel("A (mean log2 TPM, Total)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\Delta D$ (relative recruitment)")
fig.tight_layout()
for ax, texts in labels:
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5))
savefig(fig, "recruitment", f"MA_plots.recruitment.tau{TAU}")
plt.show()
""")

md(r"""
## Recruitment of selected genes relative to the controls

Points: per replicate, relative recruitment $D$ of the stress sample (PD replicate $i$ minus Total replicate $i$, centred over genes) minus the control average; markers: $\Delta D$ with 95% credible interval (`sanity_utils.plotting.plot_gene_log2fc`).
""")

code(r"""
# per-sample relative recruitment: PD replicate i minus Total replicate i of the same cohort (replicate order)
E_tot, E_pd = sanity["Total"]["result"].log2_expression, sanity["PD"]["result"].log2_expression
meta_t, meta_p = sanity["Total"]["metadata"], sanity["PD"]["metadata"]
D_sample_raw, D_meta = {}, []
for coh in sorted(meta_t["cohort"].unique()):
    t = meta_t.loc[meta_t["cohort"] == coh].sort_values("replicate")
    p = meta_p.loc[meta_p["cohort"] == coh].sort_values("replicate")
    for i, ((_, ts), (_, ps)) in enumerate(zip(t.iterrows(), p.iterrows()), start=1):
        sid = f"{coh}|{i}"
        D_sample_raw[sid] = E_pd.loc[genes_rec, ps["sample"]] - E_tot.loc[genes_rec, ts["sample"]]
        D_meta.append({"sample": sid, "cohort": coh, "pair": f"PD {ps['replicate']} / Total {ts['replicate']}"})
D_sample_raw = pd.DataFrame(D_sample_raw)
D_meta = pd.DataFrame(D_meta)

res_D = SimpleNamespace(log2_expression=D_sample_raw - D_sample_raw.median(), params={})
avg_D = SimpleNamespace(means=D_avg, errors=D_se, errors_absolute=D_se)
sc_rows = sc_comparisons.copy()
sc_rows["condition"] = sc_rows["test"].str.split(";").str[1]
sc_rows["time_point"] = sc_rows["test"].str.split(";").str[0]
sc_rows = sc_rows.assign(o=sc_rows["condition"].map(list(cat_color_dict).index)).sort_values(["o", "time_point"])
lfc_comparisons = [(f"{r.condition} vs {r.reference.split(';')[1]}, {r.time_point}", r.test, r.reference) for r in sc_rows.itertuples()]
lfc_palette = {label: cat_color_dict[test.split(";")[1]] for label, test, _ in lfc_comparisons}
sns.set(font_scale=1, style="white")
fig = su_plotting.plot_gene_log2fc(res_D, avg_D, D_meta, goi_ids, lfc_comparisons, sample_col="sample", cond_col="cohort",
                                   ci=0.95, gene_labels=gene_names.to_dict(), palette=lfc_palette, ncols=4, panel_size=(2.8, 4.0),
                                   savefig_path=os.path.join(SU_FIGS, "recruitment", "over_genes",
                                                             "_".join(genes_of_interest_list) + ".delta_D_vs_control.png"))
for ax in fig.axes:
    if ax.get_xlabel():
        ax.set_xlabel(r"$\Delta D$ vs control")
figure_caption(os.path.join(SU_FIGS, "recruitment", "over_genes", "_".join(genes_of_interest_list) + ".delta_D_vs_control.png"),
               formats=("png",))
plt.show()
""")

md(r"""
# Gene-set analysis of the recruitment changes

ISMARA-style activities (as in `pbody_RNA_Spang.TotalDE.ipynb`) fitted on the **per-replicate changes of the relative recruitment**: for each stress-vs-control comparison and replicate pair $i$, $\Delta D_{gi} = (A^{PD} - A^{Total})_{\text{stress},i} - (A^{PD} - A^{Total})_{\text{control},i}$ (PD and Total paired within a cohort, stress and control paired in replicate order). Only the column constant is removed (`center="columns"`), which also removes $\Delta\log_2 F$: the activity of a set is the shift of the recruitment of its genes relative to the other genes.

**One fit per database**, as in `pbody_RNA_Spang.TotalDE.ipynb`: literature (cytoplasmic ribosomal-protein genes, `LIT:RP`), GO:BP, GO:CC, GO:MF, Reactome, KEGG, WikiPathways (5–500 genes each). Sets of different databases overlap strongly and would share their signal in a joint fit. The gene biotypes of the annotation (vs protein-coding genes) and the mitochondrially encoded genes (`BIOTYPE:`) are included in **every** fit as covariates; the combined tables report them from the literature fit. `padj` is BH over the sets of the same fit; $z_m$ and $\lambda$ are per fit.
""")

code(r"""
from sanity_utils.genesets import GeneSetCollection

GS_SUB = "recruitment_gene_sets"
DE_FIGS = os.path.join(SU_FIGS, GS_SUB)
os.makedirs(DE_FIGS, exist_ok=True)
savefig_gs = lambda fig, name: savefig(fig, GS_SUB, name)
save_tsv_gs = lambda df, name: save_tsv(df, GS_SUB, name)

dD_cols, dD_meta = {}, []
for _, row in sc_comparisons.iterrows():
    t_ids = D_meta.loc[D_meta["cohort"] == row["test"], "sample"].tolist()
    r_ids = D_meta.loc[D_meta["cohort"] == row["reference"], "sample"].tolist()
    for i, (ts, rs) in enumerate(zip(t_ids, r_ids), start=1):
        pid = f"{row['label']}|{i}"
        dD_cols[pid] = D_sample_raw[ts] - D_sample_raw[rs]
        dD_meta.append({"sample": pid, "comparison": row["label"]})
dD_values = pd.DataFrame(dD_cols)
dD_meta = pd.DataFrame(dD_meta)
print(dD_values.shape[1], "replicate pairs on", len(dD_values), "genes")

universe = set(dD_values.index)
biotype_sets = {f"BIOTYPE:{b}": set(g.index) & universe for b, g in genes_df.groupby("gene_biotype") if b != "protein_coding"}
biotype_sets["BIOTYPE:mitochondrially_encoded"] = set(genes_df.index[genes_df["chrom"] == "Mito"]) & universe
biotype_sets = {k: v for k, v in biotype_sets.items() if len(v) >= 5}
biotype_coll = GeneSetCollection(biotype_sets, {k: k.split(":")[1].replace("_", " ") + " (biotype)" for k in biotype_sets},
                                 "biotype")
lit_coll = GeneSetCollection({"LIT:RP": set(gene_names.index[gene_names.str.match(r"^RP[LS]\d+[AB]?$")]) & universe},
                             {"LIT:RP": "cytoplasmic ribosomal-protein genes, in-house"}, "literature")

# one ISMARA fit per database (sets of different databases overlap and would share their signal in a joint fit),
# each with the biotype sets as covariates; padj: BH over the sets of the same fit
GENESET_DATABASES = ["literature", "GO:BP", "GO:CC", "GO:MF", "Reactome", "KEGG", "WikiPathways"]
db_colls = {db: lit_coll if db == "literature" else
            load_gene_sets(GENESETS_DIR, sources=(db,), genes=dD_values.index, min_size=5, max_size=500, drop_duplicates=True)
            for db in GENESET_DATABASES}

fits, act_diff, info_all, zm_all, pair_acts = {}, [], [], [], []
for db, c_db in db_colls.items():
    coll_db = c_db + biotype_coll
    a = su.fit_set_activities(dD_values, coll_db, center="columns", ridge="cv", cv_repeats=5)
    a_avg = a.average_conditions(dD_meta, sample_col="sample", cond_col="comparison", method="reml")
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
bt_fits = act_diff_all.loc[act_diff_all["set_id"].str.startswith("BIOTYPE:")].pivot_table(
    index=["set_id", "comparison"], columns="fit", values="z_score")
print("biotype z-scores across fits: max range", round(float((bt_fits.max(axis=1) - bt_fits.min(axis=1)).max()), 2),
      "| median range", round(float((bt_fits.max(axis=1) - bt_fits.min(axis=1)).median()), 2))
save_tsv_gs(act_diff_df, "recruitment.gene_set_activities.delta_D_vs_control.tsv")
save_tsv_gs(zm_df, "recruitment.gene_set_activities.z_m.tsv")
de_results = {k: d.rename(columns={"delta_D_centered": "log2FC"}) for k, d in rec_results.items()}   # for the biotype view
display(zm_df.head(15).round(2))
""")

md(r"""
### Sets changed in at least one stress, and in all stresses

Activity z-scores of the stress-vs-control comparisons (stress vs control only). `database` gives the source (= the fit) of each set: literature, GO:BP, GO:CC, GO:MF, Reactome, KEGG, WikiPathways, and biotype.
""")

code(r"""
# (1) sets with padj < 0.05 in at least one stress-vs-control comparison, clustered
sc = act_diff_df.loc[act_diff_df["kind"] == "stress_vs_control"].copy()
sc_order = [c for c in comparisons_df["label"] if c in set(sc["comparison"])]
set_label = (set_info.set_index("set_id")["name"].str.slice(0, 60) + " (" + set_info.set_index("set_id")["database"] + ")")
print(sc.loc[sc["padj"] < 0.05].groupby("comparison").size().reindex(sc_order).rename("n sets padj < 0.05").to_string())

sig_any = sc.loc[sc["padj"] < 0.05, "set_id"].unique()
mat_any = sc.loc[sc["set_id"].isin(sig_any)].pivot(index="set_id", columns="comparison", values="z_score")[sc_order]
save_tsv_gs(mat_any.rename_axis("set_id").reset_index().merge(set_info, on="set_id"),
         "recruitment.gene_set_activities.z.sig_in_any_stress.tsv")
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
    save_tsv_gs(mat_all.assign(same_direction=same_dir).rename_axis("set_id").reset_index().merge(set_info, on="set_id"),
             f"recruitment.gene_set_activities.z.sig_in_all_stresses.{tp}.tsv")
    print(f"{tp}: {len(mat_all)} sets with padj < 0.05 in all five stresses, {int(same_dir.sum())} in the same direction")
    lim = np.nanpercentile(np.abs(mat_all.values), 98)
    sns.set(font_scale=0.9, style="white")
    fig, ax = plt.subplots(figsize=(0.5 * len(cols) + 7, min(0.28 * len(mat_all) + 1.5, 60)))
    sns.heatmap(mat_all.rename(index=lambda s: ("" if same_dir[s] else "≠ ") + set_label[s]), cmap="RdBu_r", center=0,
                vmin=-lim, vmax=lim, ax=ax, cbar_kws={"label": "activity z-score"}, yticklabels=len(mat_all) <= 200)
    ax.set(xlabel="", ylabel="", title=f"{tp}: gene sets changed in all stresses (≠: not all in the same direction)")
    savefig_gs(fig, f"gene_set_activities.sig_in_all_stresses.{tp}.heatmap")
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
# ISMARA activity (left; from the literature fit, see the fit cell) and median ΔD of the biotype's genes (right, gene level)
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
axes[0].set(xlabel="", ylabel="", title="ISMARA activity of gene biotypes vs protein-coding\n(numbers: activity = shift of ΔD; * padj < 0.05)")
lim = np.nanmax(np.abs(med.values))
sns.heatmap(med, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".2f", annot_kws={"fontsize": 7},
            cbar_kws={"label": "median ΔD"}, ax=axes[1])
axes[1].set(xlabel="", ylabel="", title="median ΔD of the genes of each biotype")
fig.tight_layout()
savefig_gs(fig, "gene_biotypes.activity_and_median_log2FC")
save_tsv_gs(bt.merge(set_info[["set_id"]], on="set_id"), "recruitment.gene_biotype_activities.tsv")
plt.show()
""")

md(r"""
# Per-gene fractions $\alpha_g$

`fracSanity_alpha_per_cohort.tsv`: per gene and cohort the posterior median of $\alpha_g$ with 95% credible interval, the posterior mean and variance of $\log_2\alpha_g$ (averaged over the posterior of $F$), and the input log2 PD/Total ratio $D_g$ with its variance. `fracSanity_alpha_matrix.tsv`: genes × cohorts.
""")

code(r"""
alpha_df = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_alpha_per_cohort.tsv"), sep="\t")
alpha_df["condition"] = pd.Categorical(alpha_df["condition"], CONDITION_ORDER, ordered=True)
alpha_df["log2_alpha_sd"] = np.sqrt(alpha_df["log2_alpha_var"])
alpha_df["log2_alpha_median"] = np.log2(alpha_df["alpha"])
alpha_df["CI_width_log2"] = np.log2(alpha_df["alpha_CI_upper"]) - np.log2(alpha_df["alpha_CI_lower"])

sns.set(font_scale=1, style="ticks")
fig, ax = plt.subplots(1, 1, figsize=(8, 3.5))
sns.boxenplot(data=alpha_df, x="condition", y="log2_alpha_mean", hue="time_point", order=CONDITION_ORDER,
              hue_order=list(mrk_marker_dict), palette={"10min": "lightgrey", "60min": "dimgrey"}, ax=ax, showfliers=False)
ax.set(ylabel=r"posterior mean $\log_2 \alpha_g$", xlabel="")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), title="time point")
savefig(fig, "fracSanity", "log2_alpha_per_cohort")
display(alpha_df.groupby(["time_point", "condition"], observed=True)[["alpha", "CI_width_log2"]].median().round(4).unstack(0))
""")

code(r"""
# genes discussed by the Spang lab (Wang et al. eLife 2018; thesis slides): validated P-body mRNAs under
# glucose starvation (BSC1, TPI1, RLM1 decayed; ATP11, ILM1, MRPL38, AIM2 stabilised), negative control MET10,
# housekeeping ACT1, PGK1, the P-body factor DCP2, and ribosomal-protein mRNAs highlighted in the slides
genes_of_interest = ["BSC1", "TPI1", "RLM1", "ATP11", "ILM1", "MRPL38", "AIM2", "MET10", "ACT1", "PGK1", "DCP2",
                     "RPS9A", "RPS9B", "RPS6A", "RPL19B"]
name_to_id = {v: k for k, v in gene_names.items()}
goi = alpha_df.loc[alpha_df["gene_name"].isin(genes_of_interest)].copy()
goi["gene_name"] = pd.Categorical(goi["gene_name"], [g for g in genes_of_interest if g in set(goi["gene_name"])], ordered=True)

sns.set(font_scale=0.9, style="ticks")
g = sns.FacetGrid(goi.sort_values(["gene_name", "condition"]), col="gene_name", col_wrap=5, height=2.3, aspect=1.25, sharey=False)
def _pp(data, **kw):
    ax = plt.gca()
    for tp, dx in (("10min", -0.15), ("60min", 0.15)):
        d = data.loc[data["time_point"] == tp]
        x = [CONDITION_ORDER.index(c) + dx for c in d["condition"]]
        ax.errorbar(x, d["alpha"], yerr=[d["alpha"] - d["alpha_CI_lower"], d["alpha_CI_upper"] - d["alpha"]],
                    fmt=mrk_marker_dict[tp], ms=4, color="black", ecolor="grey", lw=0.8, capsize=1)
        ax.scatter(x, d["alpha"], c=[cat_color_dict[c] for c in d["condition"]], marker=mrk_marker_dict[tp], s=25, zorder=3, edgecolor="black", lw=0.5)
    ax.set_xticks(range(len(CONDITION_ORDER)))
    ax.set_xticklabels(CONDITION_ORDER, rotation=90)
g.map_dataframe(_pp)
g.set_titles("{col_name}")
g.set_axis_labels("", r"$\alpha_g$ (95% CI)")
savefig(g.fig, "fracSanity", "alpha_genes_of_interest")
""")

md(r"""
# Differential recruitment (stress vs control)

Per gene and comparison: $\Delta\log_2\alpha_g$ (posterior mean, including the uncertainty of $F$ of both conditions), 95% credible interval, `prob_increase` $= P(\Delta > 0)$ and the local false sign rate `lfsr`. Calls: $|\Delta\log_2\alpha_g| > \tau = 0.5$ and lfsr < 0.05.
""")

code(r"""
TAU = 0.5
LFSR = 0.05
diff_df = pd.read_csv(os.path.join(SU_TABLES, "fracSanity", "fracSanity_differential_recruitment.tsv"), sep="\t")
diff_df["call"] = np.select(
    [(diff_df["delta_log2_alpha"] > TAU) & (diff_df["lfsr"] < LFSR), (diff_df["delta_log2_alpha"] < -TAU) & (diff_df["lfsr"] < LFSR)],
    ["increased", "decreased"], default="n.s.")
diff_df["time_point"] = diff_df["comparison"].str.extract(r"^(\d+min);")[0]
diff_df["condition"] = diff_df["comparison"].str.extract(r";(\w+)_vs_")[0]
sc = diff_df.loc[diff_df["kind"] == "stress_vs_control"].copy()
counts_calls = sc.groupby(["comparison", "call"]).size().unstack(fill_value=0)
counts_calls = counts_calls.reindex(stress_comparisons["label"])
save_tsv(counts_calls.reset_index(), "fracSanity", f"n_genes_differentially_recruited.tau{TAU}.tsv")
display(counts_calls)
display(diff_df.loc[diff_df["kind"] == "time"].groupby(["comparison", "call"]).size().unstack(fill_value=0))
""")

code(r"""
sns.set(font_scale=0.9, style="ticks")
labels = list(stress_comparisons["label"])
fig, axes = plt.subplots(2, 6, figsize=(20, 7), sharex=True, sharey=True)
pal = {"increased": "red", "decreased": "royalblue", "n.s.": "lightgrey"}
rp_mask = diff_df["gene_name"].fillna("").str.match(r"^RP[LS]\d+[AB]?$")
for ax, lab in zip(axes.flat, labels):
    d = diff_df.loc[diff_df["comparison"] == lab].copy()
    d["-log10 lfsr"] = -np.log10(d["lfsr"].clip(lower=1e-12))
    sns.scatterplot(data=d, x="delta_log2_alpha", y="-log10 lfsr", hue="call", palette=pal, s=5, alpha=0.6,
                    linewidth=0, ax=ax, legend=False, hue_order=["n.s.", "decreased", "increased"])
    rp = d.loc[rp_mask.loc[d.index]]
    ax.scatter(rp["delta_log2_alpha"], rp["-log10 lfsr"], s=6, facecolor="none", edgecolor="black", lw=0.5, label="RP genes")
    ax.axvline(TAU, ls="--", c="grey", lw=0.7)
    ax.axvline(-TAU, ls="--", c="grey", lw=0.7)
    ax.axhline(-np.log10(LFSR), ls="--", c="grey", lw=0.7)
    n_up, n_dn = (d["call"] == "increased").sum(), (d["call"] == "decreased").sum()
    ax.set_title(f"{lab}\n↑ {n_up}  ↓ {n_dn}", fontsize=9)
    ax.set(xlabel=r"$\Delta\log_2\alpha$", ylabel=r"$-\log_{10}$ lfsr", xlim=(-4, 4))
axes.flat[0].legend(loc="upper left", fontsize=7)
savefig(fig, "fracSanity", f"volcano_differential_recruitment.tau{TAU}")
""")

md(r"""
# Overlaps between stresses

Genes **more recruited** into P-bodies (Δlog2α > τ, lfsr < 0.05) in each stress relative to its own control, per time point. Matrix layout ("UpSet", Lex & Gehlenborg, Nat Methods 2014) for the five stresses, and three-set Venn diagrams for the stresses with the untreated control (minusD, DTT, H2O2) and for those with the DMSO control (CCCP, Tun) together with the vehicle effect itself (DMSO vs UT).
""")

code(r"""
def gene_sets_from_calls(time_point, call="increased", conditions=STRESSES + ["DMSO"]):
    out = {}
    for cond in conditions:
        lab = f"{time_point};{cond}_vs_{CONTROL_OF[cond]}"
        d = diff_df.loc[(diff_df["comparison"] == lab) & (diff_df["call"] == call)]
        out[cond] = set(d["gene_id"])
    return out

overlap_summary, membership_tables = [], []
for tp in TIME_POINTS:
    for call in ("increased", "decreased"):
        sets_tp = gene_sets_from_calls(tp, call)
        five = {k: sets_tp[k] for k in STRESSES}
        union = set().union(*five.values())
        in_all = set.intersection(*five.values())
        n_sets = Counter(g for s in five.values() for g in s)
        specific = {k: {g for g in v if n_sets[g] == 1} for k, v in five.items()}
        overlap_summary.append({"time_point": tp, "call": call, "n_union": len(union), "n_in_all_5_stresses": len(in_all),
                                "n_in_>=3_stresses": sum(v >= 3 for v in n_sets.values()),
                                **{f"n_{k}": len(v) for k, v in five.items()},
                                **{f"n_{k}_specific": len(v) for k, v in specific.items()}})
        m = pd.DataFrame({k: [g in v for g in sorted(union)] for k, v in five.items()}, index=sorted(union))
        m["n_stresses"] = m.sum(axis=1)
        m["also_DMSO_vs_UT"] = [g in sets_tp["DMSO"] for g in m.index]
        m.insert(0, "gene_name", m.index.map(gene_names))
        m.insert(0, "call", call)
        m.insert(0, "time_point", tp)
        membership_tables.append(m.reset_index().rename(columns={"index": "gene_id"}))
overlap_summary_df = pd.DataFrame(overlap_summary)
membership_df = pd.concat(membership_tables)
save_tsv(overlap_summary_df, "overlaps", f"overlap_summary.tau{TAU}.tsv")
save_tsv(membership_df, "overlaps", f"gene_membership_across_stresses.tau{TAU}.tsv")
display(overlap_summary_df.T)
""")

code(r"""
for tp in TIME_POINTS:
    for call in ("increased", "decreased"):
        sets_tp = gene_sets_from_calls(tp, call)
        contents = {k: sets_tp[k] for k in STRESSES if len(sets_tp[k])}
        if len(contents) < 2:
            continue
        sns.set(font_scale=1, style="white")
        fig = upset_plot(contents, {k: cat_color_dict[k] for k in contents})
        fig.suptitle(f"{tp}: genes {call} in P-bodies (Δlog2α {'>' if call == 'increased' else '<'} {'' if call == 'increased' else '-'}{TAU}, lfsr < {LFSR})", fontsize=11)
        savefig(fig, "overlaps", f"upset.{call}.{tp}.tau{TAU}")
        plt.show()
""")

code(r"""
for call in ("increased", "decreased"):
    sns.set(font_scale=1, style="white")
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for i, tp in enumerate(TIME_POINTS):
        sets_tp = gene_sets_from_calls(tp, call)
        ut = ["minusD", "DTT", "H2O2"]
        v = venn3([sets_tp[k] for k in ut], set_labels=ut, ax=axes[i, 0],
                  set_colors=[cat_color_dict[k] for k in ut], alpha=0.6)
        axes[i, 0].set_title(f"{tp}, control UT")
        dm = ["CCCP", "Tun", "DMSO"]
        venn3([sets_tp[k] for k in dm], set_labels=["CCCP vs DMSO", "Tun vs DMSO", "DMSO vs UT"], ax=axes[i, 1],
              set_colors=[cat_color_dict[k] for k in dm], alpha=0.6)
        axes[i, 1].set_title(f"{tp}, control DMSO")
    fig.suptitle(f"genes {call} in P-bodies (|Δlog2α| > {TAU}, lfsr < {LFSR})")
    savefig(fig, "overlaps", f"venn3.{call}.tau{TAU}")
    plt.show()
""")

md(r"""
# Functional enrichment

Three complementary analyses:

1. **ISMARA-style gene-set activities** (Sanity_utils `fit_set_activities`; [docs](https://github.com/nimwegenLab/Sanity_utils/blob/main/docs/gene_set_analysis.md)), fitted on the log2 TPM of all PD and Total samples together: per set and cohort the recruitment activity (PD − Total), and its change between stress and control, with errors from replicate averaging. Overlapping sets compete for the same genes, and replicate-to-replicate co-variation of the genes of a set enters the error bars.
2. **Over-representation** (hypergeometric) of the gene lists of the overlap analysis, e.g. genes more recruited in all five stresses, or in one stress only. Universe: all genes with a fracSanity estimate. Genes are treated as independent, so the p-values are optimistic; they describe the lists.
3. **CAMERA-like** competitive test on the gene-level z-scores of differential recruitment (inter-gene correlation 0.01).

Gene sets: GO (biological process, cellular component, molecular function; SGD annotations propagated over is_a/part_of), Reactome, KEGG and WikiPathways, 5–500 genes in the universe; sets with identical members merged.
""")

code(r"""
universe = sorted(set(diff_df["gene_id"]))
coll = load_gene_sets(GENESETS_DIR, sources=("GO:BP", "GO:CC", "GO:MF", "Reactome", "KEGG", "WikiPathways"),
                      genes=universe, min_size=5, max_size=500, drop_duplicates=True)
set_table = set_source_table(coll)
set_table["size"] = set_table["set_id"].map(coll.sizes)
print(len(coll), "gene sets;", set_table["database"].value_counts().to_dict())

def set_label(set_id, maxlen=55):
    name = str(coll.names.get(set_id, set_id)).split(" | ")[0]
    return (name[:maxlen] + "…" if len(name) > maxlen else name) + f" ({set_id})"
""")

code(r"""
# 1. ISMARA-style activities from the pipeline
sa_dir = os.path.join(SU_TABLES, "gene_set_activities")
sa_diff = pd.read_csv(os.path.join(sa_dir, "differential_recruitment.tsv"), sep="\t").rename(columns={"set": "set_id", "index": "set_id"})
sa_rec = pd.read_csv(os.path.join(sa_dir, "recruitment_per_cohort.tsv"), sep="\t").rename(columns={"set": "set_id", "index": "set_id"})
sa_expr = pd.read_csv(os.path.join(sa_dir, "differential_expression_total.tsv"), sep="\t").rename(columns={"set": "set_id", "index": "set_id"})
cv = pd.read_csv(os.path.join(sa_dir, "cv_curve.tsv"), sep="\t")
fit_ps = pd.read_csv(os.path.join(sa_dir, "fit_per_sample.tsv"), sep="\t")
lam = cv.loc[cv["test_mse"].idxmin(), "lambda"]

sns.set(font_scale=1, style="ticks")
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
axes[0].errorbar(cv["lambda"], cv["test_mse"], yerr=cv["test_mse_sd"], fmt="o-", ms=3, color="black", ecolor="grey")
axes[0].axvline(lam, color="red", ls="--")
axes[0].set(xscale="log", xlabel=r"$\lambda$", ylabel="test MSE (80/20 CV)")
sns.histplot(fit_ps["explained_variance"], ax=axes[1], bins=20, color="grey")
axes[1].set(xlabel="fraction of variance explained per sample")
savefig(fig, "enrichment", "set_activities.cv_and_fit")

sa_sc = sa_diff.loc[sa_diff["kind"] == "stress_vs_control"].copy()
z = sa_sc.pivot(index="set_id", columns="comparison", values="z_score")[list(stress_comparisons["label"])]
padj = sa_sc.pivot(index="set_id", columns="comparison", values="padj")[list(stress_comparisons["label"])]
top = z.abs().max(axis=1).sort_values(ascending=False).index[:40]
zt = z.loc[top].copy()
zt.index = [set_label(s, 45) for s in zt.index]
sns.set(font_scale=0.8, style="white")
cg = sns.clustermap(zt, cmap="RdBu_r", center=0, vmin=-8, vmax=8, col_cluster=False, figsize=(10, 12),
                    cbar_kws={"label": "z (differential recruitment activity)"}, yticklabels=True)
savefig(cg.fig, "enrichment", "set_activities.differential_recruitment.top40_heatmap")
save_tsv(sa_sc, "enrichment", "set_activities.differential_recruitment.stress_vs_control.tsv")
""")

code(r"""
# ribosome, ribosome biogenesis, TOR and other focus sets
focus_patterns = {
    "cytosolic ribosome": r"^cytosolic ribosome$|cytosolic large ribosomal subunit|cytosolic small ribosomal subunit",
    "structural constituent of ribosome": r"structural constituent of ribosome",
    "KEGG Ribosome": r"^Ribosome$",
    "ribosome biogenesis": r"^ribosome biogenesis$|Ribosome biogenesis in eukaryotes|^rRNA processing$",
    "TORC1 signaling": r"TORC1 signaling|TOR signaling",
    "P-body": r"^P-body$",
    "mitochondrial respiration": r"oxidative phosphorylation|respiratory electron transport chain|^aerobic respiration$",
    "glycolysis": r"glycolytic process$|Glycolysis",
}
focus_sets = {}
for label, pat in focus_patterns.items():
    hits = [k for k in coll.sets if re.search(pat, str(coll.names[k]), flags=re.IGNORECASE)]
    focus_sets[label] = hits
print({k: [(h, coll.names[h][:40], len(coll[h])) for h in v][:4] for k, v in focus_sets.items()})
focus_ids = sorted({h for v in focus_sets.values() for h in v})

d = sa_sc.loc[sa_sc["set_id"].isin(focus_ids)].copy()
d["set"] = d["set_id"].map(lambda s: set_label(s, 40))
d["time_point"] = d["comparison"].str.extract(r"^(\d+min);")[0]
d["condition"] = d["comparison"].str.extract(r";(\w+)_vs_")[0]
sns.set(font_scale=0.85, style="ticks")
sets_order = list(dict.fromkeys(d.sort_values("set")["set"]))
fig, axes = plt.subplots(1, 2, figsize=(14, 0.35 * len(sets_order) + 1.5), sharey=True)
for ax, tp in zip(axes, TIME_POINTS):
    dd = d.loc[d["time_point"] == tp]
    conds = [c for c in CONTROL_OF if c in set(dd["condition"])]
    for j, cond in enumerate(conds):
        x = dd.loc[dd["condition"] == cond]
        y = [sets_order.index(s) + (j - len(conds) / 2) * 0.12 for s in x["set"]]
        ax.errorbar(x["delta_activity"], y, xerr=1.96 * x["SE"], fmt="o", ms=4, color=cat_color_dict[cond],
                    mec="black", mew=0.4, ecolor=cat_color_dict[cond], lw=0.8, label=f"{cond} vs {CONTROL_OF[cond]}")
    ax.axvline(0, color="grey", lw=0.8)
    ax.set_yticks(range(len(sets_order)))
    ax.set_yticklabels(sets_order)
    ax.set(title=tp, xlabel="Δ recruitment activity (stress − control, ±95%)")
axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), title="")
savefig(fig, "enrichment", "set_activities.focus_sets.differential_recruitment")
save_tsv(d, "enrichment", "set_activities.focus_sets.differential_recruitment.tsv")

r = sa_rec.loc[sa_rec["set_id"].isin(focus_ids)].copy()
r["set"] = r["set_id"].map(lambda s: set_label(s, 40))
save_tsv(r, "enrichment", "set_activities.focus_sets.recruitment_per_cohort.tsv")
display(r.pivot(index="set", columns="cohort", values="z_score").round(1))
""")

code(r"""
# ribosomal-protein (RP) and ribosome-biogenesis (RiBi) mRNAs at the gene level
rp_genes = set(gene_names.index[gene_names.fillna("").str.match(r"^RP[LS]\d+[AB]?$")]) & set(universe)
ribi_ids = [k for k in coll.sets if re.search(r"^ribosome biogenesis$", str(coll.names[k]).split(" | ")[0], flags=re.I)]
ribi_genes = (set().union(*[coll[k] for k in ribi_ids]) - rp_genes) if ribi_ids else set()
print(len(rp_genes), "RP genes;", len(ribi_genes), "RiBi genes (GO:0042254, without RP genes)")
diff_df["gene_class"] = np.where(diff_df["gene_id"].isin(rp_genes), "RP",
                         np.where(diff_df["gene_id"].isin(ribi_genes), "RiBi", "other"))

rows = []
for lab in diff_df["comparison"].unique():
    d = diff_df.loc[diff_df["comparison"] == lab]
    for cls in ("RP", "RiBi"):
        a, b = d.loc[d["gene_class"] == cls, "delta_log2_alpha"], d.loc[d["gene_class"] == "other", "delta_log2_alpha"]
        rows.append({"comparison": lab, "class": cls, "n": len(a), "median_class": a.median(), "median_other": b.median(),
                     "median_difference": a.median() - b.median(), "mannwhitney_p": mannwhitneyu(a, b).pvalue,
                     "frac_increased": (d.loc[d["gene_class"] == cls, "call"] == "increased").mean(),
                     "frac_increased_other": (d.loc[d["gene_class"] == "other", "call"] == "increased").mean()})
rp_tests = pd.DataFrame(rows)
rp_tests["padj"] = multipletests(rp_tests["mannwhitney_p"], method="fdr_bh")[1]
save_tsv(rp_tests, "enrichment", "RP_RiBi_vs_other.delta_log2_alpha.tsv")
display(rp_tests.round(4))

sns.set(font_scale=0.9, style="ticks")
d = diff_df.loc[diff_df["kind"] == "stress_vs_control"].copy()
fig, ax = plt.subplots(1, 1, figsize=(14, 4))
sns.boxplot(data=d, x="comparison", y="delta_log2_alpha", hue="gene_class", order=list(stress_comparisons["label"]),
            hue_order=["other", "RiBi", "RP"], palette={"other": "lightgrey", "RiBi": "orange", "RP": "firebrick"},
            showfliers=False, ax=ax)
ax.axhline(0, color="black", lw=0.8)
ax.set(ylabel=r"$\Delta\log_2\alpha$ (stress − control)", xlabel="")
ax.tick_params(axis="x", rotation=45)
for i, lab in enumerate(stress_comparisons["label"]):
    t = rp_tests.loc[(rp_tests["comparison"] == lab) & (rp_tests["class"] == "RP")]
    if len(t):
        ax.text(i + 0.27, ax.get_ylim()[1] * 0.92, p_to_stars(t["padj"].iloc[0]), ha="center", fontsize=9, color="firebrick")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="")
savefig(fig, "enrichment", "RP_RiBi_delta_log2_alpha")

# basal recruitment (log2 alpha) of RP genes in each cohort
alpha_df["gene_class"] = np.where(alpha_df["gene_id"].isin(rp_genes), "RP",
                          np.where(alpha_df["gene_id"].isin(ribi_genes), "RiBi", "other"))
fig, ax = plt.subplots(1, 1, figsize=(14, 3.5))
alpha_df["cohort_label"] = alpha_df["time_point"] + ";" + alpha_df["condition"].astype(str)
order = [f"{t};{c}" for t in TIME_POINTS for c in CONDITION_ORDER]
sns.boxplot(data=alpha_df, x="cohort_label", y="log2_alpha_mean", hue="gene_class", order=order,
            hue_order=["other", "RiBi", "RP"], palette={"other": "lightgrey", "RiBi": "orange", "RP": "firebrick"},
            showfliers=False, ax=ax)
ax.tick_params(axis="x", rotation=45)
ax.set(ylabel=r"$\log_2\alpha$", xlabel="")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), title="")
savefig(fig, "enrichment", "RP_RiBi_log2_alpha_per_cohort")
""")

code(r"""
# 2. over-representation of the overlap gene lists
ora_results = []
for tp in TIME_POINTS:
    m = membership_df.loc[(membership_df["time_point"] == tp) & (membership_df["call"] == "increased")]
    lists = {"in_all_5_stresses": m.loc[m["n_stresses"] == 5, "gene_id"],
             "in_>=3_stresses": m.loc[m["n_stresses"] >= 3, "gene_id"]}
    for k in STRESSES:
        lists[f"{k}_specific"] = m.loc[(m["n_stresses"] == 1) & m[k], "gene_id"]
        lists[f"{k}_all"] = m.loc[m[k], "gene_id"]
    for name, genes in lists.items():
        if len(genes) < 5:
            continue
        r = su.overrepresentation(list(genes), universe, coll, min_size=5, max_size=500)
        r.insert(0, "list", name)
        r.insert(0, "time_point", tp)
        r.insert(2, "n_list", len(genes))
        ora_results.append(r.reset_index())
ora_df = pd.concat(ora_results) if ora_results else pd.DataFrame()
save_tsv(ora_df.loc[ora_df["padj"] < 0.1] if len(ora_df) else ora_df, "enrichment", "ORA.increased_gene_lists.padj0.1.tsv")

if len(ora_df):
    sig = ora_df.loc[(ora_df["padj"] < 0.05) & (ora_df["overlap"] >= 3)].copy()
    top_terms = (sig.sort_values("p_value").groupby(["time_point", "list"]).head(6)["set_id"].unique())
    dd = ora_df.loc[ora_df["set_id"].isin(top_terms) & (ora_df["padj"] < 0.05)].copy()
    if len(dd):
        dd["term"] = dd["set_id"].map(lambda s: set_label(s, 50))
        dd["column"] = dd["time_point"] + " | " + dd["list"]
        dd["-log10 padj"] = -np.log10(dd["padj"].clip(lower=1e-30))
        sns.set(font_scale=0.8, style="whitegrid")
        fig, ax = plt.subplots(1, 1, figsize=(0.45 * dd["column"].nunique() + 6, 0.3 * dd["term"].nunique() + 2))
        sns.scatterplot(data=dd, x="column", y="term", size="-log10 padj", hue="fold_enrichment", palette="viridis",
                        sizes=(15, 200), ax=ax, edgecolor="black")
        ax.tick_params(axis="x", rotation=90)
        ax.set(xlabel="", ylabel="")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
        savefig(fig, "enrichment", "ORA.increased_gene_lists.dotplot")
    display(sig.sort_values("p_value").groupby(["time_point", "list"]).head(5)[["time_point", "list", "n_list", "name", "size", "overlap", "fold_enrichment", "padj"]])
""")

code(r"""
# 3. CAMERA-like test on gene-level z-scores, and agreement with the set activities
cam_rows = []
for lab in stress_comparisons["label"]:
    d = diff_df.loc[diff_df["comparison"] == lab].set_index("gene_id")["z_score"]
    r = su.camera(d, coll, inter_gene_cor=0.01)
    r.insert(0, "comparison", lab)
    cam_rows.append(r.reset_index())
camera_df = pd.concat(cam_rows)
save_tsv(camera_df.loc[camera_df["padj"] < 0.1], "enrichment", "CAMERA.stress_vs_control.padj0.1.tsv")

m = pd.merge(camera_df[["comparison", "set_id", "t", "padj"]],
             sa_sc[["comparison", "set_id", "z_score", "padj"]], on=["comparison", "set_id"], suffixes=("_camera", "_ismara"))
rho = m.groupby("comparison").apply(lambda x: spearmanr(x["t"], x["z_score"])[0])
display(rho.rename("Spearman(CAMERA t, ISMARA z)").round(2))
display(camera_df.loc[camera_df["set_id"].isin(focus_ids)].pivot(index="set_id", columns="comparison", values="t").round(1)
        .rename(index=lambda s: set_label(s, 40)))
""")

md(r"""
# Stress specificity vs P-body properties and translation

Per stress and time point: global $F$, number of more/less recruited genes, and the differential recruitment activity of the cytosolic-ribosome set, against P-body brightness and number (imaging), and against the translation shutoff (L-HPG incorporation relative to control) and TORC1 activity (phospho-Rps6 / Rps6) reported by the Spang lab. The latter two are approximate values read from figures (manuscript Fig. 4B/C; thesis slides 18–19) and not available for every condition; with at most 10 points the correlations are descriptive.
""")

code(r"""
lit = pd.DataFrame([
    # condition, time, translation (L-HPG rel. to control), TORC1 (pRps6/Rps6 rel. to control); approximate
    ("minusD", "10min", 0.04, 0.03), ("H2O2", "10min", 0.04, 0.27), ("CCCP", "10min", 0.35, 0.29),
    ("DTT", "10min", 0.73, 1.10), ("Tun", "10min", 0.93, 0.88),
    ("minusD", "60min", np.nan, 0.0), ("H2O2", "60min", np.nan, 0.0), ("CCCP", "60min", 0.37, 1.10),
    ("DTT", "60min", 0.44, 0.34), ("Tun", "60min", 0.95, 0.65),
], columns=["condition", "time_point", "translation_rel_control", "TORC1_pRps6_rel_control"])
lit["source"] = "Mookherjee, Rommel et al. bioRxiv 2026 (Fig. 4) and Lesslauer MSc defense 2026 (slides 18-19); read from plots"

ribo_set = [k for k in focus_sets.get("cytosolic ribosome", [])][:1]
rows = []
for tp in TIME_POINTS:
    for cond in STRESSES:
        lab = f"{tp};{cond}_vs_{CONTROL_OF[cond]}"
        r = {"time_point": tp, "condition": cond,
             "F": globalF_df.loc[globalF_df["cohort"] == f"{tp};{cond}", "F"].iloc[0],
             "F_control": globalF_df.loc[globalF_df["cohort"] == f"{tp};{CONTROL_OF[cond]}", "F"].iloc[0],
             "n_increased": int((diff_df.loc[diff_df["comparison"] == lab, "call"] == "increased").sum()),
             "n_decreased": int((diff_df.loc[diff_df["comparison"] == lab, "call"] == "decreased").sum())}
        if ribo_set:
            x = sa_sc.loc[(sa_sc["comparison"] == lab) & (sa_sc["set_id"] == ribo_set[0])]
            r["ribosome_set_delta_activity_z"] = x["z_score"].iloc[0] if len(x) else np.nan
        rows.append(r)
spec = pd.DataFrame(rows)
spec["log2_F_over_control"] = np.log2(spec["F"] / spec["F_control"])
spec = spec.merge(pbody_intensity_gr_df[["time_point", "condition", "pb/cyto"]], how="left")
spec = spec.merge(pbody_number_gr_df[["time_point", "condition", "pbody_abundance"]], how="left")
spec = spec.merge(lit, how="left")
save_tsv(spec, "stress_specificity", "per_stress_summary.tsv")
display(spec.round(3))

xs = ["pb/cyto", "pbody_abundance", "translation_rel_control", "TORC1_pRps6_rel_control"]
ys = ["F", "log2_F_over_control", "n_increased", "ribosome_set_delta_activity_z"]
ys = [y for y in ys if y in spec.columns]
sns.set(font_scale=0.85, style="ticks")
fig, axes = plt.subplots(len(ys), len(xs), figsize=(3.1 * len(xs), 2.7 * len(ys)))
for i, y in enumerate(ys):
    for j, x in enumerate(xs):
        ax = axes[i, j]
        dd = spec[[x, y, "condition", "time_point"]].dropna()
        for _, rr in dd.iterrows():
            ax.scatter(rr[x], rr[y], color=cat_color_dict[rr["condition"]], marker=mrk_marker_dict[rr["time_point"]],
                       s=50, edgecolor="black")
        if len(dd) >= 4:
            rho_, p_ = spearmanr(dd[x], dd[y])
            ax.set_title(f"ρ = {rho_:.2f} (p = {p_:.2g}), n = {len(dd)}", fontsize=8)
        ax.set(xlabel=x if i == len(ys) - 1 else "", ylabel=y if j == 0 else "")
savefig(fig, "stress_specificity", "per_stress_summary_scatter")
""")

md(r"""
# Sequence and structure determinants of recruitment

Per-gene features (`gene_features.tsv`, pipeline rule `gene_features`): lengths and GC of 5'UTR, CDS and 3'UTR (UTRs from TIF-seq longest full-ORF transcripts, so the 5' ends are approximate), GC3 and tRNA adaptation index (tAI, codon optimality), the 5' end (first nucleotide, pyrimidine content, mammalian-style TOP motif: C followed by ≥ 4 pyrimidines), RNAfold minimum free energies (5'UTR, start-codon window, 3'UTR), and counts of RBP motifs (Puf3, Puf4, Whi3, Ssd1). In vivo DMS-seq structure (Rouskin et al. Nature 2014, lifted to R64; rule `structure_features`): window Gini of A/C reactivities in vivo minus that of the denatured control, per region (UTRs capped to 120 nt upstream of the start and 150 nt downstream of the stop codon, where the DMS signal is), for genes with ≥ 15 reads per A/C. Expression in Total (log2 TPM) is included as a covariate because many features correlate with it.

Relation to recruitment: univariate Spearman correlations with $\Delta\log_2\alpha$ (stress vs control) and with the basal $\log_2\alpha$ (UT), and a multivariate weighted regression ($1/\mathrm{Var}$ weights) on standardized features.
""")

code(r"""
feat = pd.read_csv(os.path.join(SU_TABLES, "gene_features", "gene_features.tsv"), sep="\t", index_col=0)
# DMS-seq structure (pipeline rule structure_features): in vivo window Gini minus that of the denatured control
# (corrects for the uneven DMS reactivity of unstructured RNA); higher = more structured in vivo
struct_path = os.path.join(SU_TABLES, "structure", "structure_features.wide.tsv")
if os.path.exists(struct_path):
    st = pd.read_csv(struct_path, sep="\t", index_col=0)
    for region in ("5UTR", "start_region", "CDS", "3UTR"):
        col = f"in_vivo_minus_denatured_{region}_delta_mean_window_gini"
        if col in st.columns:
            feat[f"gini_struct_{region}"] = st[col].reindex(feat.index)
    print("structure columns:", {c: int(feat[c].notna().sum()) for c in feat.columns if c.startswith("gini_")})
else:
    print("structure features not available:", struct_path)

feat["log_len_5utr"] = np.log10(feat["len_5utr"].clip(lower=1))
feat["log_len_cds"] = np.log10(feat["len_cds"])
feat["log_len_3utr"] = np.log10(feat["len_3utr"].clip(lower=1))
for k in ("Puf3", "Puf4", "Whi3", "Ssd1"):
    feat[f"has_{k}_3utr"] = (feat[f"n_{k}_3utr"] > 0).astype(float)
feat["top_like"] = feat["top_like"].astype(float)
feat["RP_gene"] = feat.index.isin(rp_genes).astype(float)
tot_avg = sanity["Total"]["averages"].means
feat["log2_TPM_Total_UT10"] = tot_avg["10min;UT"].reindex(feat.index)

feature_cols = ["log_len_5utr", "log_len_cds", "log_len_3utr", "gc_5utr", "gc_cds", "gc_3utr", "gc3", "tai",
                "pyr_frac_first10", "top_like", "mfe_5utr_per_nt", "mfe_start_window_per_nt", "mfe_3utr_per_nt",
                "has_Puf3_3utr", "has_Puf4_3utr", "has_Whi3_3utr", "has_Ssd1_3utr", "RP_gene", "log2_TPM_Total_UT10"]
feature_cols += [c for c in feat.columns if c.startswith("gini_")]
save_tsv(feat.reset_index().rename(columns={"index": "gene_id"}), "features", "gene_features_with_structure.tsv")
display(feat[feature_cols].describe().T.round(3))

# RP vs other genes
cmp = feat.groupby(feat["RP_gene"].map({1.0: "RP", 0.0: "other"}))[feature_cols].median().T
display(cmp.round(3))
""")

code(r"""
targets = {}
for lab in stress_comparisons["label"]:
    targets[lab] = diff_df.loc[diff_df["comparison"] == lab].set_index("gene_id")[["delta_log2_alpha", "posterior_variance"]]
for cohort in ("10min;UT", "60min;UT", "10min;DMSO", "60min;DMSO"):
    a = alpha_df.loc[alpha_df["cohort"] == cohort].set_index("gene_id")
    targets[f"basal log2α {cohort}"] = a[["log2_alpha_mean", "log2_alpha_var"]].rename(
        columns={"log2_alpha_mean": "delta_log2_alpha", "log2_alpha_var": "posterior_variance"})

uni = pd.DataFrame({lab: {c: spearmanr(*feat[[c]].join(t["delta_log2_alpha"]).dropna().to_numpy().T)[0]
                          for c in feature_cols} for lab, t in targets.items()})
sns.set(font_scale=0.8, style="white")
fig, ax = plt.subplots(1, 1, figsize=(0.55 * uni.shape[1] + 3, 0.35 * uni.shape[0] + 1.5))
sns.heatmap(uni, cmap="RdBu_r", center=0, vmin=-0.4, vmax=0.4, annot=True, fmt=".2f", annot_kws={"fontsize": 6},
            cbar_kws={"label": "Spearman ρ"}, ax=ax)
savefig(fig, "features", "features_vs_recruitment.spearman")
save_tsv(uni.reset_index().rename(columns={"index": "feature"}), "features", "features_vs_recruitment.spearman.tsv")
""")

code(r"""
# multivariate weighted least squares on standardized features (genes with all features)
coef_rows = []
base_cols = [c for c in feature_cols if not c.startswith("gini_")]
for lab, t in targets.items():
    d = feat[base_cols].join(t).dropna()
    X = (d[base_cols] - d[base_cols].mean()) / d[base_cols].std()
    X = sm.add_constant(X)
    w = 1.0 / d["posterior_variance"].clip(lower=d["posterior_variance"].quantile(0.05))
    res = sm.WLS(d["delta_log2_alpha"], X, weights=w).fit(cov_type="HC3")
    for c in base_cols:
        coef_rows.append({"target": lab, "feature": c, "coef": res.params[c], "se": res.bse[c], "p_value": res.pvalues[c],
                          "n_genes": len(d), "R2": res.rsquared})
coef_df = pd.DataFrame(coef_rows)
coef_df["padj"] = multipletests(coef_df["p_value"], method="fdr_bh")[1]
save_tsv(coef_df, "features", "features_vs_recruitment.WLS_coefficients.tsv")

cm = coef_df.pivot(index="feature", columns="target", values="coef").loc[base_cols, list(targets)]
stars = coef_df.assign(s=coef_df["padj"].map(p_to_stars).replace("ns", "")).pivot(index="feature", columns="target", values="s").loc[base_cols, list(targets)]
sns.set(font_scale=0.8, style="white")
fig, ax = plt.subplots(1, 1, figsize=(0.55 * cm.shape[1] + 3, 0.35 * cm.shape[0] + 1.5))
sns.heatmap(cm, cmap="RdBu_r", center=0, annot=stars, fmt="", cbar_kws={"label": "coefficient (log2 units per s.d.)"}, ax=ax)
r2 = coef_df.drop_duplicates("target").set_index("target")["R2"]
ax.set_title("WLS, standardized features; R² per target: " + ", ".join(f"{v:.2f}" for v in r2[list(targets)]), fontsize=8)
savefig(fig, "features", "features_vs_recruitment.WLS_coefficients")
""")

code(r"""
# TOP-like 5' ends and RBP motifs among genes more recruited in each stress (Fisher's exact test)
rows = []
for lab in stress_comparisons["label"]:
    d = diff_df.loc[diff_df["comparison"] == lab].set_index("gene_id").join(feat.drop(columns="gene_name", errors="ignore"), how="inner")
    up = d["call"] == "increased"
    for f in ["top_like", "has_Puf3_3utr", "has_Puf4_3utr", "has_Whi3_3utr", "has_Ssd1_3utr", "RP_gene"]:
        x = d[f].astype(bool)
        tab = [[(up & x).sum(), (up & ~x).sum()], [(~up & x).sum(), (~up & ~x).sum()]]
        orat, p = fisher_exact(tab)
        rows.append({"comparison": lab, "feature": f, "n_up_with": tab[0][0], "n_up": int(up.sum()),
                     "frac_up_with": tab[0][0] / max(1, up.sum()), "frac_other_with": tab[1][0] / max(1, (~up).sum()),
                     "odds_ratio": orat, "p_value": p})
fisher_df = pd.DataFrame(rows)
fisher_df["padj"] = multipletests(fisher_df["p_value"], method="fdr_bh")[1]
save_tsv(fisher_df, "features", "motifs_TOP_RP_in_increased_genes.fisher.tsv")
display(fisher_df.pivot(index="feature", columns="comparison", values="odds_ratio").round(2))
display(fisher_df.pivot(index="feature", columns="comparison", values="padj").round(4))
""")

md(r"""
# Translation (public Ribo-seq)

Public yeast ribosome-profiling data with matched RNA-seq, processed with `WF/Snakefile-Riboseq` / `WF/Snakefile-RNAseq` (adapted from the Arf1_TORC1 project): glucose starvation 10 min (Kershaw et al. RNA Biol 2021) and 5/15/30 min (Zedan et al. bioRxiv 2024), H2O2 30 min (Blevins et al. Sci Rep 2019), DTT and tunicamycin 1 h (van Dalfsen et al. Dev Cell 2018). No CCCP data exist. Ribosome footprints (P-sites in CDS) and RNA counts are normalised with Sanity per assay, and the log2 translation efficiency of a condition is $\langle A^{Ribo}\rangle - \langle A^{RNA}\rangle$ (the fracSanity-style difference, up to a constant per condition); its change between stress and control is compared with $\Delta\log_2\alpha$ of the closest stress and time point of our data.
""")

code(r"""
ribo_root = os.path.join(subdirs["tables_dir"], "riboseq_public")
te_results = {}
if os.path.isdir(ribo_root):
    for ds in sorted(os.listdir(ribo_root)):
        p = os.path.join(ribo_root, ds)
        if not all(os.path.exists(os.path.join(p, f)) for f in ("ribo_counts.tsv", "rna_counts.tsv", "metadata.tsv")):
            continue
        meta = pd.read_csv(os.path.join(p, "metadata.tsv"), sep="\t")
        out = {}
        for assay, fname in (("Ribo", "ribo_counts.tsv"), ("RNA", "rna_counts.tsv")):
            cnt = pd.read_csv(os.path.join(p, fname), sep="\t", index_col=0)
            cnt = cnt.loc[cnt.index.isin(universe)]
            cnt = cnt.loc[cnt.sum(axis=1) > 0]
            m = meta.loc[meta["sample"].isin(cnt.columns)]
            res = su.sanity_normalization(cnt[m["sample"].tolist()], n_threads=4)
            out[assay] = su.average_conditions(res, m, sample_col="sample", cond_col="condition", method="reml")
        genes_c = out["Ribo"].means.index.intersection(out["RNA"].means.index)
        te = out["Ribo"].means.loc[genes_c] - out["RNA"].means.loc[genes_c]
        te_err = np.sqrt(out["Ribo"].errors.loc[genes_c] ** 2 + out["RNA"].errors.loc[genes_c] ** 2)
        te_results[ds] = {"te": te, "te_err": te_err, "meta": meta}
        print(ds, te.shape, list(te.columns))
else:
    print("no Ribo-seq tables yet in", ribo_root)
""")

code(r"""
# pairs: (dataset, stress condition, control condition) -> our comparison(s)
te_pairs = [
    ("Kershaw2021", "minusD_10min", "untreated", ["10min;minusD_vs_UT"]),
    ("Zedan2024", "minusD_5min", "untreated", ["10min;minusD_vs_UT"]),
    ("Zedan2024", "minusD_15min", "untreated", ["10min;minusD_vs_UT"]),
    ("Zedan2024", "minusD_30min", "untreated", ["10min;minusD_vs_UT", "60min;minusD_vs_UT"]),
    ("Blevins2019", "H2O2_30min", "untreated", ["10min;H2O2_vs_UT", "60min;H2O2_vs_UT"]),
    ("vanDalfsen2018", "DTT_1h", "untreated", ["60min;DTT_vs_UT"]),
    ("vanDalfsen2018", "Tm_1h", "untreated", ["60min;Tun_vs_DMSO"]),
]
rows, merged_all = [], []
for ds, stress, ctrl, ours in te_pairs:
    if ds not in te_results:
        continue
    te = te_results[ds]["te"]
    if stress not in te.columns or ctrl not in te.columns:
        print(ds, "missing", stress, ctrl, list(te.columns))
        continue
    dte = (te[stress] - te[ctrl]).rename("delta_log2_TE")
    for lab in ours:
        d = diff_df.loc[diff_df["comparison"] == lab].set_index("gene_id")[["gene_name", "delta_log2_alpha", "call", "gene_class"]]
        j = d.join(dte, how="inner").join(te[ctrl].rename("log2_TE_control"), how="inner").dropna()
        rho1, p1 = spearmanr(j["delta_log2_TE"], j["delta_log2_alpha"])
        rho2, p2 = spearmanr(j["log2_TE_control"], j["delta_log2_alpha"])
        up = j["call"] == "increased"
        rows.append({"dataset": ds, "riboseq_stress": stress, "our_comparison": lab, "n_genes": len(j),
                     "spearman_dTE_vs_dlog2alpha": rho1, "p_dTE": p1,
                     "spearman_TEcontrol_vs_dlog2alpha": rho2, "p_TEcontrol": p2,
                     "median_dTE_increased": j.loc[up, "delta_log2_TE"].median(),
                     "median_dTE_other": j.loc[~up, "delta_log2_TE"].median(),
                     "mannwhitney_p_dTE": mannwhitneyu(j.loc[up, "delta_log2_TE"], j.loc[~up, "delta_log2_TE"]).pvalue if up.sum() > 2 else np.nan})
        merged_all.append(j.assign(dataset=ds, riboseq_stress=stress, our_comparison=lab).reset_index())
if rows:
    te_summary = pd.DataFrame(rows)
    save_tsv(te_summary, "translation", "deltaTE_vs_delta_log2_alpha.summary.tsv")
    te_merged = pd.concat(merged_all)
    save_tsv(te_merged, "translation", "deltaTE_vs_delta_log2_alpha.per_gene.tsv")
    display(te_summary.round(4))
    n = len(te_summary)
    sns.set(font_scale=0.85, style="ticks")
    fig, axes = plt.subplots(1, n, figsize=(3.3 * n, 3.3), squeeze=False)
    for ax, (_, r) in zip(axes[0], te_summary.iterrows()):
        j = te_merged.loc[(te_merged["dataset"] == r["dataset"]) & (te_merged["riboseq_stress"] == r["riboseq_stress"])
                          & (te_merged["our_comparison"] == r["our_comparison"])]
        ax.hexbin(j["delta_log2_TE"], j["delta_log2_alpha"], gridsize=40, cmap="Greys", mincnt=1, bins="log")
        rp = j.loc[j["gene_class"] == "RP"]
        ax.scatter(rp["delta_log2_TE"], rp["delta_log2_alpha"], s=6, color="firebrick", label="RP")
        ax.axhline(0, color="grey", lw=0.6)
        ax.axvline(0, color="grey", lw=0.6)
        ax.set(xlabel=f"Δlog2 TE ({r['dataset']}, {r['riboseq_stress']})", ylabel=f"Δlog2α ({r['our_comparison']})")
        ax.set_title(f"ρ = {r['spearman_dTE_vs_dlog2alpha']:.2f} (p = {r['p_dTE']:.1g})", fontsize=9)
    axes[0][0].legend(loc="lower left", fontsize=7)
    savefig(fig, "translation", "deltaTE_vs_delta_log2_alpha")
""")

md(r"""
# Summary

See the final section of the analysis report (filled in after the run).
""")

caveat_start = len(cells)
md(r"""
# Caveats and considerations

Status of 2026-10-26. Pipeline variant `Sanity_utils_allTPM`: all expressed genes, length correction (log2 TPM; gene lengths = summed exon lengths of `YEAST_ANNOTATION_FILE`, Ensembl R64-1-1.109, from featureCounts), ISMARA replicate averages (`reml`).

## 1. The global fraction $F$ is weakly identified

- PD and Total are normalized separately (each to $10^6$ TPM), so the input $D_g = \langle A^{PD}_g\rangle - \langle A^{Total}_g\rangle$ is defined only up to a constant per cohort. In $D_g = \log_2\alpha_g - \log_2 F + \varepsilon_g$ this constant is absorbed by $F$; $F$ is pinned only by the shape of the Beta distribution of $\alpha_g$ with the bound $\alpha_g \le 1$.
- Consequently `F_identifiable` is False in every cohort (profile likelihood flat within the 1.92 threshold), and the 95% credible intervals of $F$ span one to two orders of magnitude. This was the case already in the pySanity analysis of `pbody_RNA_Spang.current.ipynb`.
- The contaminated error model assigns most genes to the outlier component in several cohorts ($\pi_{outlier} \approx 0.95$–$1$, e.g. 10 min H2O2, 10 min minusD, 10 min CCCP), and Beta parameters reach the bounds of the optimizer ($a = 500$ or $b = 0.02$) in some cohorts: the Beta model of $\alpha_g$ does not describe the distribution of $D_g$ well there, and $F$ of these cohorts is the least reliable (10 min CCCP: $\hat F = 0.26$, interval ~0–0.65).
- External information on $F$ would make it identifiable: spike-ins added before the pull-down, a measured yield of the pull-down, or an informative prior on $F$ from imaging (`prior_F` of `fit_frac_sanity`).

## 2. $F$ depends on the normalization (gene universe)

Reproduction of `globalF_vs_pBodyIntensity.v2.pdf` (figures in `figures/.../Sanity_utils/validation/`):

| run | Spearman ρ (p), $F$ vs mean PB brightness, n = 10 | same sign ΔF vs Δbrightness, 60 vs 10 min |
|---|---|---|
| pySanity + `frac_sanity_v2` (reference) | 0.20 (0.58) | 4 / 4 |
| pySanity inputs + Sanity_utils fit | 0.20 (0.58), identical $F$ | 4 / 4 |
| C++ Sanity, all genes, CPM | 0.22 (0.53) | 4 / 4 |
| C++ Sanity, all genes, TPM (this notebook) | 0.33 (0.35) | 4 / 4 |
| C++ Sanity, protein-coding genes, TPM | −0.08 (0.83) | 3 / 4 |

(minusD is not counted: its brightness does not change, Δ = −0.07.)

- The fracSanity fit of Sanity_utils is identical to `frac_sanity_v2`; the differences between runs come from the normalization only. At the gene level the inputs agree (correlation of $D_g$ between runs 0.97–0.997); what changes is the per-cohort offset of $D_g$, i.e. $F$.
- Restricting the genes to nuclear protein-coding genes changes $F$ strongly (rRNA and other non-coding RNAs are a large and variable part of the reads, and their share differs between PD and Total); length correction (TPM vs CPM) changes $F$ little. Which universe is the right reference for $F$ is a modelling choice, not a technicality; all expressed genes reproduce the earlier analysis.
- Composition of the counted reads (section "Analysis of read mapping stats", medians): Total 10 / 60 min: nuclear protein-coding 51 / 36%, rRNA 22 / 41%, sno/sn/ncRNA 10%, transposable elements and pseudogenes 12%; PD 10 / 60 min: nuclear protein-coding 33 / 30%, rRNA 7 / 15%, **transposable elements (Ty retrotransposon mRNAs) and pseudogenes 52 / 44%**. Ty mRNAs are strongly enriched in the pull-down (Ty1 RNA associates with P-body components, Checkley et al. MCB 2010), and, being among the most abundant PD RNAs, they weigh heavily on the PD normalization and hence on $F$ with the all-genes universe. rRNA content varies widely between Total libraries (8–76%).

## 3. Comparison with the imaging ground truth

- Cross-condition correlation of $F$ with P-body brightness is positive but not significant (ρ = 0.33, p = 0.35, n = 10), with large uncertainty of $F$ (section 1). P-body brightness and the fraction of mRNA in P-bodies need not be proportional across stresses (P-body number and size also differ: weak stresses give many dim, persistent P-bodies, strong stresses few bright ones).
- Within each stress, the change of $F$ from 10 to 60 min has the same sign as the change of brightness in all 4 stresses with a brightness change (DTT, Tun up; H2O2, CCCP down). Under random signs the probability is 1/16. Only the CCCP change of $F$ has a 95% interval excluding 0, and it rests on the unstable 10 min CCCP fit.
- There are no brightness data for UT and DMSO (the stresses' controls), so the controls are not part of the comparison. The reference analysis excluded 60 min DMSO from the scatter; here it is absent anyway because it has no brightness data.
- Brightness is summarized by the mean over P-bodies (heavy-tailed distributions; medians are printed on the violins).
- **Brightness vs $F$ (60 vs 10 min).** On a common log2 scale the changes of mean P-body brightness are small (|log2FC| ≤ 0.32) while $\Delta\log_2\hat F$ ranges from −6.6 to +1.8 with intervals of ±5 log2 units: the signs agree, the magnitudes do not, and the points lie far from the diagonal. $\hat F$ is not a quantitative proxy of P-body brightness.
- **Recruitment of individual mRNAs** is reported as the change of the centred $D_g$ (F-independent; "Recruitment of individual mRNAs relative to the controls"); per-replicate pairing of PD with Total and of stress with control uses the replicate order within each cohort (60 min CCCP Total: R4–R6).

## 4. Per-gene fractions and differential recruitment

- $\alpha_g$ and $\Delta\log_2\alpha_g$ inherit the uncertainty and the normalization dependence of $F$: a shift of $\log_2 F$ moves all genes of a cohort together. Calls of "more / less recruited" genes based on a fixed threshold on $\Delta\log_2\alpha_g$ can therefore be global (all or no genes), as in the earlier protein-coding run. The ranking of genes within a comparison is robust; absolute values and thresholds are not.
- Quantities independent of $F$: $D_g$ centered per cohort (relative recruitment), its differences between cohorts after centering, and the ISMARA-style gene-set and motif activities of the pipeline (samples are centered before the fit).
- The pull-down (Dcp2-HBH cCLAP, formaldehyde crosslinking) enriches RNA crosslinked to Dcp2-containing complexes, not only P-bodies; the PD libraries also contain background RNA. Wang et al. (eLife 2018) defined P-body mRNAs as PD(stress) vs PD(unstressed), without Total: their lists are not directly comparable with $\alpha_g$.

## 5. Batch and replicate structure

- The 10 and 60 min samples were sequenced as separate batches: in the Total fraction even UT 60 vs 10 min has many differentially expressed genes. Comparisons within a time point (stress vs its control) are unaffected; comparisons over time (60 vs 10 min, including ΔF) are confounded with the batch.
- CCCP Total 60 min replicates R1–R3 are excluded (as in the earlier analysis). 10 min Tun and 10 min DMSO have large replicate variability (large errors of their averages).
- In the Total fraction, the inclusion of rRNA in the normalization shifts the log2FC of all genes in some comparisons (e.g. 10 min Tun vs DMSO, median −0.62), see `pbody_RNA_Spang.TotalDE.ipynb`.

## 6. Annotation of UTRs

- UTRs come from the TIF-seq "longest full-ORF" transcripts (Pelechano et al. 2013), YPD and galactose merged (most upstream start, most downstream end): 5' ends are the most upstream TSS, and UTRs are longer than those of the dominant isoforms (median 5'UTR 117 nt, 3'UTR 196 nt). The TOP-like 5'-end feature and 5'UTR lengths are approximate; the conditions (YPD/galactose, unstressed) differ from the stresses here.

## 7. Structure data (DMS-seq)

- Rouskin et al. (Nature 2014): unstressed BY4741 in log phase (and in vitro refolded, denatured, ATP-depleted), not the stresses of this study; lifted from sacCer2 to R64.
- Coverage limits the gene set (≥ 15 reads per A/C): ~2,160 genes with a CDS value, ~1,170–1,400 for the start region, ~950 for 5'UTRs and ~1,140 for 3'UTRs, biased towards highly expressed genes. The DMS signal is concentrated near the start and stop codons, so UTRs are capped (120 nt upstream of the start, 150 nt downstream of the stop codon). Positions of non-unique 25-mers are masked (affects paralog pairs such as many ribosomal-protein genes).
- Gini indices increase with low coverage; use the within-gene difference to the denatured control (`in_vivo_minus_denatured_*_delta_mean_window_gini`) and include expression as a covariate.

## 8. Motifs

- Motif counts grow with region length; the region length is a covariate of the motif-activity fits, and expression should be one in follow-up tests.
- Short consensus motifs (Nab3 UCUU, Nrd1 GUAR) are frequent and mainly reflect base composition; UAUAUA (Hrp1) is also the polyadenylation efficiency element, so 3'UTR signals near the poly(A) site may reflect 3'-end processing or annotation rather than P-body targeting.
- The MEME motif databases contain only four yeast RBP motifs (Hrp1, Nab2, Rna15, Vts1, two of them in both files); most yeast RBP specificities are represented by consensus sequences from the literature. 4,096 6-mers per region require multiple-testing correction and replication across stresses.
- Motif and structure associations are correlational; testing a mechanism needs perturbations (e.g. deletions of the RBP, mutated reporters), as for Slf1 and the P-body brightness in the thesis.

## 9. Not yet done

- Translation (public Ribo-seq) is not processed: matching data exist for glucose starvation and H2O2, for DTT and tunicamycin only at ≥ 1 h, and none for CCCP.
- The sections of this notebook after the $F$ comparison (per-gene fractions, overlaps, enrichment, features, translation) are not yet re-run with the all-genes TPM variant.
""")

if len(sys.argv) > 2:  # e.g. "# Per-gene fractions": keep only the cells before this heading, and the caveats
    stop = next(i for i, c in enumerate(cells) if c.cell_type == "markdown" and c.source.startswith(sys.argv[2]))
    cells = cells[:stop] + cells[caveat_start:]
nb = nbf.v4.new_notebook()
nb["cells"] = cells
nb["metadata"] = {
    "kernelspec": {"display_name": "Python (sanity_utils)", "language": "python", "name": "sanity_utils"},
    "language_info": {"name": "python"},
}
nbf.write(nb, sys.argv[1])
print("written", sys.argv[1], len(cells), "cells")
