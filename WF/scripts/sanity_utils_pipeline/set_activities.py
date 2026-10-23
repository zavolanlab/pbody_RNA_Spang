"""
ISMARA-style gene-set activities (Sanity_utils.fit_set_activities) on the log2 TPM of all Total and
PD samples together, replicate averages per fraction;time;condition, and contrasts:

- recruitment of a set into the PD fraction in a cohort c:   (PD;c) - (Total;c)
- differential recruitment, test t vs reference r:            [(PD;t) - (Total;t)] - [(PD;r) - (Total;r)]
- expression change in the Total fraction:                     (Total;t) - (Total;r)

Fitting both fractions with one lambda keeps the model linear in the data, so these contrasts are
exactly the activities of the corresponding contrasts of log2 values (e.g. of log2 PD/Total ratios).
"""

import argparse
import os

import pandas as pd
import sanity_utils as su
from sanity_utils.genesets import load_gene_sets, set_source_table

from common import comparisons, dump_pickle, load_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--total", required=True)
    p.add_argument("--pd", required=True)
    p.add_argument("--genesets", required=True, help="store of sanity-utils-fetch-genesets")
    p.add_argument("--genes", required=True, help="genes.tsv (gene names)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-size", type=int, default=5)
    p.add_argument("--max-size", type=int, default=500)
    p.add_argument("--averaging-method", default="reml")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    tot, pdf = load_pickle(a.total), load_pickle(a.pd)
    genes = tot["result"].log2_expression.index.intersection(pdf["result"].log2_expression.index)
    values = pd.concat([tot["result"].log2_expression.loc[genes], pdf["result"].log2_expression.loc[genes]], axis=1)
    meta = pd.concat([tot["metadata"], pdf["metadata"]])

    coll = load_gene_sets(a.genesets, sources=("GO:BP", "GO:CC", "GO:MF", "Reactome", "KEGG", "WikiPathways"),
                          genes=genes, min_size=a.min_size, max_size=a.max_size, drop_duplicates=True)
    src = set_source_table(coll)
    src["size"] = src["set_id"].map(coll.sizes)
    names = pd.read_csv(a.genes, sep="\t", index_col=0)["gene_name"]
    src["genes"] = [",".join(sorted(names.get(g, g) for g in coll[k])) for k in src["set_id"]]
    src.to_csv(os.path.join(a.outdir, "gene_sets.tsv"), sep="\t", index=False)
    print(f"{len(coll)} gene sets on {len(genes)} genes", flush=True)

    act = su.fit_set_activities(values, coll, center="both", ridge="cv", cv_repeats=5)
    print(f"lambda = {act.lambda_:.4g}; median explained variance {act.explained_variance.median():.3f}", flush=True)
    avg = act.average_conditions(meta, sample_col="sample", cond_col="group", method=a.averaging_method)

    act.activities.to_csv(os.path.join(a.outdir, "activities_per_sample.tsv"), sep="\t", index_label="set_id")
    act.errors.to_csv(os.path.join(a.outdir, "activity_errors_per_sample.tsv"), sep="\t", index_label="set_id")
    act.cv_curve.to_csv(os.path.join(a.outdir, "cv_curve.tsv"), sep="\t", index=False)
    pd.DataFrame({"chi2": act.chi2, "explained_variance": act.explained_variance}).to_csv(
        os.path.join(a.outdir, "fit_per_sample.tsv"), sep="\t", index_label="sample")
    avg.means.to_csv(os.path.join(a.outdir, "activities_per_group.tsv"), sep="\t", index_label="set_id")
    avg.errors.to_csv(os.path.join(a.outdir, "activity_errors_per_group.tsv"), sep="\t", index_label="set_id")

    cohorts = sorted(meta["cohort"].unique())
    rec = []
    for c in cohorts:
        r = su.activity_contrast(avg, {f"PD;{c}": 1.0, f"Total;{c}": -1.0})
        r.insert(0, "cohort", c)
        rec.append(r.reset_index())
    pd.concat(rec).to_csv(os.path.join(a.outdir, "recruitment_per_cohort.tsv"), sep="\t", index=False)

    diff, expr = [], []
    for _, row in comparisons().iterrows():
        t, r = row["test"], row["reference"]
        d = su.activity_contrast(avg, {f"PD;{t}": 1.0, f"Total;{t}": -1.0, f"PD;{r}": -1.0, f"Total;{r}": 1.0})
        d.insert(0, "comparison", row["label"])
        d.insert(1, "kind", row["kind"])
        diff.append(d.reset_index())
        e = su.activity_contrast(avg, {f"Total;{t}": 1.0, f"Total;{r}": -1.0})
        e.insert(0, "comparison", row["label"])
        e.insert(1, "kind", row["kind"])
        expr.append(e.reset_index())
    pd.concat(diff).to_csv(os.path.join(a.outdir, "differential_recruitment.tsv"), sep="\t", index=False)
    pd.concat(expr).to_csv(os.path.join(a.outdir, "differential_expression_total.tsv"), sep="\t", index=False)
    dump_pickle({"activities": act, "averages": avg, "collection": coll}, os.path.join(a.outdir, "set_activities.pkl"))


if __name__ == "__main__":
    main()
