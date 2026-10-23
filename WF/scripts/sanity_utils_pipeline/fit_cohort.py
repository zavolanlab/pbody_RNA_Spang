"""fracSanity fit of one cohort (time;condition): global fraction F and per-gene alpha_g."""

import argparse
import os
import time

import sanity_utils as su

from common import dump_pickle, load_pickle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--total", required=True, help="sanity.pkl of the Total fraction")
    p.add_argument("--pd", required=True, help="sanity.pkl of the PD fraction")
    p.add_argument("--cohort", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()

    tot, pdf = load_pickle(a.total)["averages"], load_pickle(a.pd)["averages"]
    genes = tot.means.index.intersection(pdf.means.index)
    inp = su.prepare_frac_sanity_inputs(tot, pdf, a.cohort, genes=genes)
    t0 = time.time()
    fit = su.fit_frac_sanity(inp["ltq_total"], inp["var_total"], inp["ltq_pd"], inp["var_pd"],
                             n_threads=a.threads, verbose=True)
    fit["cohort"] = a.cohort
    fit["gene_ids"] = inp["gene_ids"]
    fit["inputs"] = inp
    dump_pickle(fit, a.out)
    print(f"{a.cohort}: F = {fit['F']:.4g} [{fit['F_CI_lower']:.4g}, {fit['F_CI_upper']:.4g}], "
          f"identifiable = {fit['F_identifiable']}, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "8")
    main()
