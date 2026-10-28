"""
Translation efficiency (TE) of one public Ribo-seq dataset with Sanity_utils fracSanity: Ribo-seq takes the
place of the pull-down and RNA-seq that of Total. Sanity normalizes the Ribo and the RNA samples
separately (counts in the trimmed CDS, length-corrected: log2 TPM over the CDS length used), replicates are
averaged per condition (ISMARA, reml), fracSanity fits D_g = <A_ribo> - <A_rna> = log2 alpha_g - log2 F per
condition, and calculate_differential_recruitment gives the change of log2 alpha_g (= log2 TE, up to the
global factor F) of each stress vs its control.

Outputs (outdir): te_per_condition.tsv (D_g = relative log2 TE and alpha_g per condition), te_global_F.tsv,
te_differential.tsv (stress vs control), sanity_{Ribo,RNA}.pkl, fits/<condition>.pkl.
"""

import argparse
import os
import pickle

import numpy as np
import pandas as pd
import sanity_utils as su


def dump(obj, path):
    with open(path, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--samples", required=True, help="samples table of the public Ribo-seq data")
    p.add_argument("--dataset", required=True)
    p.add_argument("--counts-dir", required=True, help="folder with <sample>.cds_counts.tsv")
    p.add_argument("--comparisons", required=True, help="test:reference[,test:reference...]")
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()
    os.makedirs(os.path.join(a.outdir, "fits"), exist_ok=True)

    meta = pd.read_csv(a.samples, sep="\t")
    meta = meta.loc[meta["dataset"] == a.dataset].copy()
    tabs = {s: pd.read_csv(os.path.join(a.counts_dir, f"{s}.cds_counts.tsv"), sep="\t", index_col=0) for s in meta["sample"]}
    counts = pd.DataFrame({s: t["count"] for s, t in tabs.items()}).fillna(0).astype(int)
    lengths = next(iter(tabs.values()))["cds_length_used"]

    res, avg = {}, {}
    for assay in ("Ribo", "RNA"):
        m = meta.loc[meta["assay"] == assay]
        c = counts[m["sample"]]
        c = c.loc[c.sum(axis=1) > 0]
        res[assay] = su.sanity_normalization(c, v_method="MAP", n_threads=a.threads, gene_lengths=lengths.loc[c.index],
                                             output_dir=os.path.join(a.outdir, f"sanity_{assay}"), overwrite=True)
        avg[assay] = su.average_conditions(res[assay], m, sample_col="sample", cond_col="condition", method="reml")
        dump({"result": res[assay], "averages": avg[assay], "metadata": m}, os.path.join(a.outdir, f"sanity_{assay}.pkl"))

    genes = avg["RNA"].means.index.intersection(avg["Ribo"].means.index)
    fits, rows, frows = {}, [], []
    for cond in sorted(meta["condition"].unique()):
        inp = su.prepare_frac_sanity_inputs(avg["RNA"], avg["Ribo"], cond, genes=genes)
        fit = su.fit_frac_sanity(inp["ltq_total"], inp["var_total"], inp["ltq_pd"], inp["var_pd"], n_threads=a.threads)
        fit.update(cohort=cond, gene_ids=list(inp["gene_ids"]), inputs=inp)
        fits[cond] = fit
        dump(fit, os.path.join(a.outdir, "fits", f"{cond}.pkl"))
        D = inp["ltq_pd"] - inp["ltq_total"]
        rows.append(pd.DataFrame({"gene_id": inp["gene_ids"], "condition": cond, "log2_TE_relative": D - np.median(D),
                                  "var_D": inp["var_pd"] + inp["var_total"], "alpha": fit["alpha_g"],
                                  "log2_alpha_mean": fit["log2_alpha_mean"], "log2_Ribo": inp["ltq_pd"],
                                  "log2_RNA": inp["ltq_total"]}))
        frows.append({"dataset": a.dataset, "condition": cond, **{k: fit.get(k) for k in
                      ("F", "F_CI_lower", "F_CI_upper", "F_identifiable", "pi_outlier", "a", "b")}})
        print(f"{a.dataset} {cond}: F = {fit['F']:.3g} [{fit['F_CI_lower']:.3g}, {fit['F_CI_upper']:.3g}]", flush=True)
    pd.concat(rows).assign(dataset=a.dataset).to_csv(os.path.join(a.outdir, "te_per_condition.tsv"), sep="\t", index=False)
    pd.DataFrame(frows).to_csv(os.path.join(a.outdir, "te_global_F.tsv"), sep="\t", index=False)

    diffs = []
    for pair in a.comparisons.split(","):
        t, r = pair.split(":")
        it, ir = fits[t]["inputs"], fits[r]["inputs"]
        d, meta_d = su.calculate_differential_recruitment(
            ltq_total_ut=ir["ltq_total"], var_total_ut=ir["var_total"], ltq_pd_ut=ir["ltq_pd"], var_pd_ut=ir["var_pd"],
            ltq_total_stress=it["ltq_total"], var_total_stress=it["var_total"], ltq_pd_stress=it["ltq_pd"],
            var_pd_stress=it["var_pd"], gene_ids=fits[t]["gene_ids"], fit_ut=fits[r], fit_stress=fits[t])
        d = d.rename(columns=lambda c: c.replace("_UT", "_reference").replace("_Stress", "_test"))
        # F-independent change: difference of the centered D_g (relative log2 TE)
        Dt, Dr = it["ltq_pd"] - it["ltq_total"], ir["ltq_pd"] - ir["ltq_total"]
        d["delta_log2_TE_relative"] = (Dt - np.median(Dt)) - (Dr - np.median(Dr))
        d["delta_log2_TE_relative_SE"] = np.sqrt(it["var_pd"] + it["var_total"] + ir["var_pd"] + ir["var_total"])
        d["padj"] = su.differential._bh(d["p_value"])
        d.insert(1, "comparison", f"{t}_vs_{r}")
        d.insert(1, "dataset", a.dataset)
        diffs.append(d)
        print(f"{t} vs {r}: delta log2 F = {meta_d.get('delta_log2_F')}", flush=True)
    pd.concat(diffs).to_csv(os.path.join(a.outdir, "te_differential.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
