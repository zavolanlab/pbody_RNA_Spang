"""
Per-gene table of DMS-seq structure metrics (output of structure_features.py features) for the
recruitment analyses: one row per gene, columns <condition>_<region>_<metric> for the conditions'
own metrics and <comparison>_<region>_<metric> for the within-gene differences between conditions.

Values of a gene/region are kept only where it passes the coverage filter (mean RT stops per A/C
>= --min-reads in structure_features.py) in the condition(s) involved. The recommended structure
measure is in_vivo_minus_denatured_<region>_delta_mean_window_gini: the in vivo window Gini corrected
for the uneven DMS reactivity of unstructured RNA of the same sequence (Rouskin et al. 2014).
"""

import argparse

import pandas as pd

METRICS = ["gini", "mean_window_gini", "mean_norm_reactivity", "n_AC", "mean_reads_per_AC"]
DIFF_METRICS = ["delta_gini", "delta_mean_window_gini"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--features", required=True, help="long TSV of structure_features.py features")
    p.add_argument("--diff", required=True, help="Gini difference TSV of structure_features.py features")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    st = pd.read_csv(a.features, sep="\t")
    st = st.loc[st["passes_coverage"].astype(str) == "True"]
    wide = st.pivot_table(index="gene_id", columns=["condition", "region"], values=METRICS, aggfunc="first")
    wide.columns = [f"{c}_{r}_{m}" for m, c, r in wide.columns]

    diff = pd.read_csv(a.diff, sep="\t")
    dwide = diff.pivot_table(index="gene_id", columns=["comparison", "region"], values=DIFF_METRICS, aggfunc="first")
    dwide.columns = [f"{c}_{r}_{m}" for m, c, r in dwide.columns]

    out = wide.join(dwide, how="outer").sort_index(axis=1)
    out.to_csv(a.out, sep="\t", index_label="gene_id")
    n = out.filter(like="in_vivo_minus_denatured").notna().sum()
    print(f"{len(out)} genes x {out.shape[1]} columns; genes with in vivo - denatured values:\n{n.to_string()}")


if __name__ == "__main__":
    main()
