"""
Validation: Sanity_utils fit_frac_sanity on the pySanity inputs of pbody_RNA_Spang.current.ipynb
(fracSanity_v2_inputs.pkl), one cohort, for comparison with the zavolab_pyutils frac_sanity_v2 fits.
"""

import argparse
import os
import pickle
import time

import sanity_utils as su

from common import dump_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--inputs", required=True, help="fracSanity_v2_inputs.pkl of the pySanity analysis")
    p.add_argument("--cohort-index", type=int, required=True, help="index into the sorted cohorts")
    p.add_argument("--outdir", required=True)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()

    with open(a.inputs, "rb") as fh:
        inputs = pickle.load(fh)
    cohort = sorted(inputs["cohorts"])[a.cohort_index]
    d = inputs["cohorts"][cohort]
    t0 = time.time()
    fit = su.fit_frac_sanity(d["ltq_total"], d["var_total"], d["ltq_pd"], d["var_pd"], prior_F=d["prior_F"],
                             n_threads=a.threads, **inputs["fit_kwargs"])
    fit["cohort"] = cohort
    fit["gene_ids"] = inputs["common_genes"]
    os.makedirs(a.outdir, exist_ok=True)
    dump_pickle(fit, os.path.join(a.outdir, cohort.replace(";", "__") + ".pkl"))
    print(f"{cohort}: F = {fit['F']:.4g} [{fit['F_CI_lower']:.4g}, {fit['F_CI_upper']:.4g}], "
          f"identifiable = {fit['F_identifiable']}, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
