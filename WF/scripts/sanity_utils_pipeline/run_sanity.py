"""
Sanity (C++) normalisation of one fraction (Total or PD) with length correction (log2 TPM), and
ISMARA replicate averages per cohort (time;condition).
"""

import argparse
import os

import pandas as pd
import sanity_utils as su

from common import dump_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", required=True, help="folder of prepare_inputs.py")
    p.add_argument("--fraction", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--v-method", default="MAP")
    p.add_argument("--averaging-method", default="reml")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    counts = pd.read_csv(os.path.join(a.inputs, "counts.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(a.inputs, "metadata.tsv"), sep="\t")
    lengths = su.read_gene_lengths(os.path.join(a.inputs, "gene_lengths.tsv"))
    meta_f = meta.loc[meta["fraction"] == a.fraction]
    sub = counts[meta_f["sample"].tolist()]
    sub = sub.loc[sub.sum(axis=1) > 0]

    res = su.sanity_normalization(sub, v_method=a.v_method, n_threads=a.threads, gene_lengths=lengths,
                                  output_dir=os.path.join(a.outdir, "sanity_run"), overwrite=True)
    avg = su.average_conditions(res, meta_f, sample_col="sample", cond_col="cohort", method=a.averaging_method)

    dump_pickle({"result": res, "averages": avg, "metadata": meta_f}, os.path.join(a.outdir, "sanity.pkl"))
    res.log2_expression.to_csv(os.path.join(a.outdir, "log2_expression.tsv"), sep="\t", index_label="gene_id")
    res.log2_sd.to_csv(os.path.join(a.outdir, "log2_sd.tsv"), sep="\t", index_label="gene_id")
    vg = res.vg if isinstance(res.vg, pd.DataFrame) else res.vg.to_frame("v_g")
    pd.concat([pd.DataFrame({"gene_mean": res.gene_mean, "gene_mean_sd": res.gene_mean_sd}),
               vg.add_prefix("vg_")], axis=1).to_csv(
        os.path.join(a.outdir, "gene_stats.tsv"), sep="\t", index_label="gene_id")
    avg.means.to_csv(os.path.join(a.outdir, "cohort_means.tsv"), sep="\t", index_label="gene_id")
    avg.errors.to_csv(os.path.join(a.outdir, "cohort_errors.tsv"), sep="\t", index_label="gene_id")
    avg.tau.to_csv(os.path.join(a.outdir, "cohort_tau.tsv"), sep="\t", index_label="gene_id")
    print(a.fraction, res.log2_expression.shape, res.params.get("units"), res.params.get("sanity_commit"))


if __name__ == "__main__":
    main()
