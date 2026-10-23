"""
Counts matrix, gene lengths and sample metadata for the Sanity_utils fracSanity pipeline.

Gene universe: nuclear protein-coding genes of the annotation (rRNA, tRNA, snoRNA etc. and
mitochondrially encoded genes are excluded, so that fractions and the global P-body fraction F refer
to cytoplasmic mRNA). Gene lengths: featureCounts 'Length' (summed exon length, the region in which
reads were counted), used for TPM-like length correction.
"""

import argparse
import glob
import os

import pandas as pd

from common import sample_metadata


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start-samples", required=True)
    p.add_argument("--featurecounts-dir", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    meta = sample_metadata(pd.read_csv(a.start_samples, sep="\t"))

    counts, lengths = {}, None
    for s in meta["sample"]:
        f = os.path.join(a.featurecounts_dir, s + ".txt")
        if not os.path.exists(f):
            raise FileNotFoundError(f)
        t = pd.read_csv(f, sep="\t", comment="#")
        t = t.set_index("Geneid")
        counts[s] = t.iloc[:, -1]
        if lengths is None:
            lengths = t["Length"]
    counts = pd.DataFrame(counts).fillna(0.0)

    # gene annotation from the GTF 'gene' records
    rows = []
    with open(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(x.strip().split(" ", 1) for x in f[8].strip().strip(";").split(";") if x.strip())
            attrs = {k: v.strip('"') for k, v in attrs.items()}
            rows.append((attrs.get("gene_id"), attrs.get("gene_name", attrs.get("gene_id")),
                         attrs.get("gene_biotype"), f[0], int(f[3]), int(f[4]), f[6]))
    genes = pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_biotype", "chrom", "start", "end",
                                        "strand"]).drop_duplicates("gene_id").set_index("gene_id")
    universe = genes.index[(genes["gene_biotype"] == "protein_coding") & (genes["chrom"] != "Mito")]
    universe = universe.intersection(counts.index)
    # counted fractional (featureCounts --fraction): round to integers for the Poisson model
    sel = counts.loc[universe].round().astype(int)
    sel = sel.loc[sel.sum(axis=1) > 0]

    genes.loc[sel.index].assign(length=lengths.reindex(sel.index)).to_csv(
        os.path.join(a.outdir, "genes.tsv"), sep="\t")
    lengths.reindex(sel.index).rename("length").to_csv(os.path.join(a.outdir, "gene_lengths.tsv"), sep="\t",
                                                       index_label="gene_id")
    sel.to_csv(os.path.join(a.outdir, "counts.tsv"), sep="\t", index_label="gene_id")
    meta.to_csv(os.path.join(a.outdir, "metadata.tsv"), sep="\t", index=False)

    excluded = counts.drop(index=universe, errors="ignore")
    frac_excl = excluded.sum() / counts.sum()
    frac_excl.rename("fraction_of_counts_excluded").to_csv(
        os.path.join(a.outdir, "excluded_counts_fraction.tsv"), sep="\t", index_label="sample")
    print(f"{sel.shape[0]} genes x {sel.shape[1]} samples; median fraction of counts outside the "
          f"universe: {frac_excl.median():.3f}")


if __name__ == "__main__":
    main()
