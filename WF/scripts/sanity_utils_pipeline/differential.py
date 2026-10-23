"""Differential recruitment (fracSanity) between two cohorts, from precomputed cohort fits."""

import argparse

import pandas as pd
import sanity_utils as su

from common import dump_pickle, load_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fit-test", required=True)
    p.add_argument("--fit-ref", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--genes", required=True, help="genes.tsv of prepare_inputs.py (names)")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-meta", required=True)
    a = p.parse_args()

    fit_t, fit_r = load_pickle(a.fit_test), load_pickle(a.fit_ref)
    if fit_t["gene_ids"] != fit_r["gene_ids"]:
        raise ValueError("Cohort fits are on different genes.")
    it, ir = fit_t["inputs"], fit_r["inputs"]
    res, meta = su.calculate_differential_recruitment(
        ltq_total_ut=ir["ltq_total"], var_total_ut=ir["var_total"], ltq_pd_ut=ir["ltq_pd"], var_pd_ut=ir["var_pd"],
        ltq_total_stress=it["ltq_total"], var_total_stress=it["var_total"], ltq_pd_stress=it["ltq_pd"],
        var_pd_stress=it["var_pd"], gene_ids=fit_t["gene_ids"], fit_ut=fit_r, fit_stress=fit_t)
    names = pd.read_csv(a.genes, sep="\t", index_col=0)["gene_name"]
    res.insert(1, "gene_name", res["gene_id"].map(names))
    res = res.rename(columns=lambda c: c.replace("_UT", "_reference").replace("_Stress", "_test"))
    res.insert(2, "comparison", a.label)
    res["padj"] = su.differential._bh(res["p_value"])
    res.to_csv(a.out_tsv, sep="\t", index=False)
    keep = {k: v for k, v in meta.items() if k not in ("fit_UT", "fit_Stress")}
    keep.update(label=a.label, cohort_test=fit_t["cohort"], cohort_reference=fit_r["cohort"])
    dump_pickle(keep, a.out_meta)
    print(a.label, "delta log2 F =", meta.get("delta_log2_F"), "P(>0) =", meta.get("prob_delta_log2_F_positive"))


if __name__ == "__main__":
    main()
