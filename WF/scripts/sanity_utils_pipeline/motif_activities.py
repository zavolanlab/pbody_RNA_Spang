"""
ISMARA-style motif activities (Sanity_utils.fit_set_activities with a motif count matrix N in place of
gene-set membership) on the log2 TPM of all Total and PD samples, for the genes with the region annotated.

Model per sample s: E_gs = c~_s + c_g + sum_m N_gm A_ms + noise, N_gm = number of sites of motif m in the
region of gene g (or log(1 + count) / presence, --transform), plus the column log10(region length) as a
covariate so that motif activities are not proxies of region length. Ridge penalty by cross-validation.
Because E is centred per sample, the per-cohort contrast PD - Total of a motif activity does not depend
on the global P-body fraction F.

Outputs (outdir): activities_per_sample.tsv, activity_errors_per_sample.tsv, activities_per_group.tsv,
activity_errors_per_group.tsv, cv_curve.tsv, fit_per_sample.tsv, recruitment_per_cohort.tsv
(PD - Total per cohort), differential_recruitment.tsv ((PD - Total) test - reference per comparison),
differential_expression_total.tsv (Total test - reference), motif_activities.pkl.
"""

import argparse
import os

import numpy as np
import pandas as pd
import sanity_utils as su

from common import comparisons, dump_pickle, load_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--total", required=True)
    p.add_argument("--pd", required=True)
    p.add_argument("--matrix", required=True, help="genes x motifs counts (motif_matrices.py)")
    p.add_argument("--lengths", required=True, help="region_lengths.tsv (motif_matrices.py)")
    p.add_argument("--region", required=True, choices=["5UTR", "CDS", "3UTR"])
    p.add_argument("--transform", default="count", choices=["count", "log1p", "presence"])
    p.add_argument("--min-genes", type=int, default=20, help="drop motifs present in fewer genes")
    p.add_argument("--outdir", required=True)
    p.add_argument("--averaging-method", default="reml")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    tot, pdf = load_pickle(a.total), load_pickle(a.pd)
    N = pd.read_csv(a.matrix, sep="\t", index_col=0)
    length = pd.read_csv(a.lengths, sep="\t", index_col=0)[a.region].dropna()
    genes = (tot["result"].log2_expression.index.intersection(pdf["result"].log2_expression.index)
             .intersection(N.index).intersection(length.index[length > 0]))
    N = N.loc[genes]
    N = N.loc[:, (N > 0).sum() >= a.min_genes]
    if a.transform == "log1p":
        N = np.log1p(N)
    elif a.transform == "presence":
        N = (N > 0).astype(float)
    N["log10_region_length"] = np.log10(length.loc[genes])

    values = pd.concat([tot["result"].log2_expression.loc[genes], pdf["result"].log2_expression.loc[genes]], axis=1)
    meta = pd.concat([tot["metadata"], pdf["metadata"]])
    print(f"{a.region}: {len(genes)} genes x {N.shape[1]} motifs ({a.transform}); {values.shape[1]} samples", flush=True)

    act = su.fit_set_activities(values, N, center="both", ridge="cv", cv_repeats=5)
    print(f"lambda = {act.lambda_:.4g}; median explained variance {act.explained_variance.median():.3f}", flush=True)
    avg = act.average_conditions(meta, sample_col="sample", cond_col="group", method=a.averaging_method)

    act.activities.to_csv(os.path.join(a.outdir, "activities_per_sample.tsv"), sep="\t", index_label="motif")
    act.errors.to_csv(os.path.join(a.outdir, "activity_errors_per_sample.tsv"), sep="\t", index_label="motif")
    act.cv_curve.to_csv(os.path.join(a.outdir, "cv_curve.tsv"), sep="\t", index=False)
    pd.DataFrame({"chi2": act.chi2, "explained_variance": act.explained_variance}).to_csv(
        os.path.join(a.outdir, "fit_per_sample.tsv"), sep="\t", index_label="sample")
    avg.means.to_csv(os.path.join(a.outdir, "activities_per_group.tsv"), sep="\t", index_label="motif")
    avg.errors.to_csv(os.path.join(a.outdir, "activity_errors_per_group.tsv"), sep="\t", index_label="motif")

    rec = []
    for c in sorted(meta["cohort"].unique()):
        r = su.activity_contrast(avg, {f"PD;{c}": 1.0, f"Total;{c}": -1.0})
        r.insert(0, "cohort", c)
        rec.append(r.rename_axis("motif").reset_index())
    pd.concat(rec).to_csv(os.path.join(a.outdir, "recruitment_per_cohort.tsv"), sep="\t", index=False)

    diff, expr = [], []
    for _, row in comparisons().iterrows():
        t, r = row["test"], row["reference"]
        d = su.activity_contrast(avg, {f"PD;{t}": 1.0, f"Total;{t}": -1.0, f"PD;{r}": -1.0, f"Total;{r}": 1.0})
        d.insert(0, "comparison", row["label"])
        d.insert(1, "kind", row["kind"])
        diff.append(d.rename_axis("motif").reset_index())
        e = su.activity_contrast(avg, {f"Total;{t}": 1.0, f"Total;{r}": -1.0})
        e.insert(0, "comparison", row["label"])
        e.insert(1, "kind", row["kind"])
        expr.append(e.rename_axis("motif").reset_index())
    pd.concat(diff).to_csv(os.path.join(a.outdir, "differential_recruitment.tsv"), sep="\t", index=False)
    pd.concat(expr).to_csv(os.path.join(a.outdir, "differential_expression_total.tsv"), sep="\t", index=False)
    dump_pickle({"activities": act, "averages": avg, "matrix": N}, os.path.join(a.outdir, "motif_activities.pkl"))


if __name__ == "__main__":
    main()
