"""
Final tables of the fracSanity pipeline:

- fracSanity_global_F.tsv: per cohort the global fraction F of mRNA in the PD fraction with 95%
  credible interval and fit diagnostics;
- fracSanity_alpha_per_cohort.tsv: per gene and cohort alpha_g (posterior median, 95% interval),
  posterior mean and variance of log2 alpha_g, the input log2 PD/Total ratio D_g and its variance;
- fracSanity_alpha_matrix.tsv: genes x cohorts posterior medians of alpha_g (wide);
- fracSanity_differential_recruitment.tsv: all comparisons, long format;
- fracSanity_delta_global_F.tsv: change of log2 F per comparison.
"""

import argparse
import glob
import os

import numpy as np
import pandas as pd

from common import comparisons, load_pickle, safe_name


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fits", nargs="+", required=True)
    p.add_argument("--diff-tsv", nargs="+", required=True)
    p.add_argument("--diff-meta", nargs="+", required=True)
    p.add_argument("--genes", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    names = pd.read_csv(a.genes, sep="\t", index_col=0)["gene_name"]

    glob_rows, alpha_rows = [], []
    for f in a.fits:
        fit = load_pickle(f)
        c = fit["cohort"]
        t, cond = c.split(";")
        lo, hi = fit.get("F_profile_CI", (np.nan, np.nan))
        glob_rows.append(dict(cohort=c, time_point=t, condition=cond, F=fit["F"], F_CI_lower=fit["F_CI_lower"],
                              F_CI_upper=fit["F_CI_upper"], F_MAP=fit.get("F_MAP"), F_MLE=fit.get("F_MLE"),
                              F_profile_CI_lower=lo, F_profile_CI_upper=hi, var_log2_F=fit.get("var_log2_F"),
                              F_identifiable=fit.get("F_identifiable"), a=fit.get("a"), b=fit.get("b"),
                              error_model=fit.get("error_model"), tau=fit.get("tau"),
                              pi_outlier=fit.get("pi_outlier"), ks_stat=fit.get("ks_stat"),
                              ks_pval=fit.get("ks_pval"), n_genes=len(fit["gene_ids"])))
        inp = fit["inputs"]
        alpha_rows.append(pd.DataFrame({
            "gene_id": fit["gene_ids"], "cohort": c, "time_point": t, "condition": cond,
            "alpha": fit["alpha_g"], "alpha_CI_lower": fit["alpha_g_CI_lower"],
            "alpha_CI_upper": fit["alpha_g_CI_upper"], "log2_alpha_mean": fit["log2_alpha_mean"],
            "log2_alpha_var": fit["log2_alpha_var"], "alpha_plugin": fit["alpha_g_plugin"],
            "outlier_prob": (np.full(len(fit["gene_ids"]), np.nan) if fit.get("outlier_prob") is None
                             else fit["outlier_prob"]),
            "log2_PD_over_Total": inp["ltq_pd"] - inp["ltq_total"], "var_D": inp["var_pd"] + inp["var_total"],
            "log2_Total": inp["ltq_total"], "log2_PD": inp["ltq_pd"],
        }))
    glob_df = pd.DataFrame(glob_rows).sort_values(["time_point", "condition"])
    glob_df.to_csv(os.path.join(a.outdir, "fracSanity_global_F.tsv"), sep="\t", index=False)
    alpha = pd.concat(alpha_rows)
    alpha.insert(1, "gene_name", alpha["gene_id"].map(names))
    alpha.to_csv(os.path.join(a.outdir, "fracSanity_alpha_per_cohort.tsv"), sep="\t", index=False)
    wide = alpha.pivot(index="gene_id", columns="cohort", values="alpha")
    wide.insert(0, "gene_name", wide.index.map(names))
    wide.to_csv(os.path.join(a.outdir, "fracSanity_alpha_matrix.tsv"), sep="\t")

    diff = pd.concat([pd.read_csv(f, sep="\t") for f in a.diff_tsv])
    kinds = comparisons().set_index("label")["kind"]
    diff.insert(3, "kind", diff["comparison"].map(kinds))
    diff.to_csv(os.path.join(a.outdir, "fracSanity_differential_recruitment.tsv"), sep="\t", index=False)
    for kind, sub in diff.groupby("kind"):
        sub.to_csv(os.path.join(a.outdir, f"fracSanity_differential_recruitment.{kind}.tsv"), sep="\t", index=False)

    rows = []
    for f in a.diff_meta:
        m = load_pickle(f)
        ci = m.get("delta_log2_F_CI", (np.nan, np.nan))
        rows.append(dict(comparison=m["label"], kind=kinds.get(m["label"]), cohort_test=m["cohort_test"],
                         cohort_reference=m["cohort_reference"], F_test=m["F_Stress"], F_reference=m["F_UT"],
                         delta_log2_F=m.get("delta_log2_F"), var_delta_log2_F=m.get("var_delta_log2_F"),
                         delta_log2_F_CI_lower=ci[0], delta_log2_F_CI_upper=ci[1],
                         prob_delta_log2_F_positive=m.get("prob_delta_log2_F_positive")))
    pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "fracSanity_delta_global_F.tsv"), sep="\t", index=False)
    print(glob_df[["cohort", "F", "F_CI_lower", "F_CI_upper", "F_identifiable"]].to_string(index=False))


if __name__ == "__main__":
    main()
