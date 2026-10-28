"""
Motif occurrence matrices (genes x motifs) per mRNA region, for motif-activity and enrichment analyses
of P-body recruitment.

Regions: 5UTR, CDS, 3UTR sequences written by gene_features.py (from CDS_and_UTRs.gtf: UTRs are the
parts of the TIF-seq longest full-ORF transcripts, YPD and galactose merged, outside the CDS; the CDS
includes the stop codon). Sequences are on the transcript strand, DNA alphabet (T = U).

Matrices (counts of overlapping occurrences on the transcript strand):
    kmer<k>.<region>.tsv.gz   all k-mers (k from --k), e.g. 4096 columns for k = 6
    known.<region>.tsv        consensus motifs of yeast RNA-binding proteins (CONSENSUS below) and sites
                              of the yeast position weight matrices in --meme (FIMO, p < --fimo-p)
    region_lengths.tsv        length of each region per gene (nt), as a covariate: motif counts grow
                              with region length
Genes without an annotated region (e.g. no 5'UTR) are absent from that region's matrices.
"""

import argparse
import gzip
import itertools
import os
import re
import subprocess

import numpy as np
import pandas as pd

# consensus motifs (DNA alphabet, transcript strand); lookahead regex counts overlapping sites
CONSENSUS = {
    "Puf3_UGUAAAUA": "TGTAAATA",           # Gerber et al. PLoS Biol 2004
    "Puf4_UGUAHAUA": "TGTA[ACT]ATA",        # Gerber et al. PLoS Biol 2004
    "Whi3_UGCAU": "TGCAT",                  # Colomina et al. EMBO J 2008; Cai & Futcher PLoS One 2013
    "Ssd1_CNYUCNYU": "C[ACGT][CT]TC[ACGT][CT]T",  # Hogan et al. PLoS Biol 2008; Bayne et al. NAR 2022
    "Nab3_UCUU": "TCTT",                    # Carroll et al. MCB 2004 (short, frequent)
    "Nrd1_GUAR": "GTA[AG]",                 # Carroll et al. MCB 2004 (short, frequent)
    "Hrp1_UAUAUA": "TATATA",                # Kessler et al. Genes Dev 1997
    "A8": "A{8}",                           # poly(A) stretch (Pab1)
    "U8": "T{8}",                           # poly(U) stretch (Pub1, Ngr1)
}
REGIONS = ("5UTR", "CDS", "3UTR")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items() if v}


def kmer_counts(seqs, k):
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=k)]
    idx = {km: i for i, km in enumerate(kmers)}
    M = np.zeros((len(seqs), len(kmers)), dtype=np.int32)
    for r, s in enumerate(seqs.values()):
        for i in range(len(s) - k + 1):
            j = idx.get(s[i:i + k])
            if j is not None:
                M[r, j] += 1
    return pd.DataFrame(M, index=list(seqs), columns=[km.replace("T", "U") for km in kmers])


def consensus_counts(seqs):
    return pd.DataFrame({m: [len(re.findall(f"(?={pat})", s)) for s in seqs.values()] for m, pat in CONSENSUS.items()},
                        index=list(seqs))


def fimo_counts(fasta, meme, pvalue, fimo_bin, exclude):
    # sites per gene and motif on the transcript strand only (--norc)
    out = subprocess.run([fimo_bin, "--norc", "--thresh", str(pvalue), "--text", "--skip-matched-sequence",
                          meme, fasta], capture_output=True, text=True, check=True).stdout
    rows = [l.split("\t") for l in out.splitlines() if l and not l.startswith(("motif_id", "#"))]
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows).iloc[:, :3]
    df.columns = ["motif_id", "motif_alt_id", "sequence_name"]
    if exclude:
        df = df.loc[~df["motif_alt_id"].str.contains(exclude, regex=True)]
    df["motif"] = "PWM_" + df["motif_alt_id"].str.replace(r"[^\w.]+", "_", regex=True) + "_" + df["motif_id"]
    return df.groupby(["sequence_name", "motif"]).size().unstack(fill_value=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta-dir", required=True, help="sequences/ of gene_features.py (5UTR.fa, CDS.fa, 3UTR.fa)")
    p.add_argument("--meme", action="append", default=[], help="MEME-format PWMs (DNA-encoded), repeatable")
    p.add_argument("--k", type=int, action="append", help="k-mer sizes (default 6)")
    p.add_argument("--fimo-p", type=float, default=1e-4)
    p.add_argument("--fimo-bin", default="fimo")
    p.add_argument("--pwm-exclude", default="Drosophila|Homo_sapiens|Mus_musculus",
                   help="regex on motif_alt_id: drop PWMs inferred from other species' RBPs")
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ks = a.k or [6]

    lengths = {}
    for region in REGIONS:
        fasta = os.path.join(a.fasta_dir, f"{region}.fa")
        seqs = read_fasta(fasta)
        lengths[region] = pd.Series({g: len(s) for g, s in seqs.items()})
        for k in ks:
            with gzip.open(os.path.join(a.outdir, f"kmer{k}.{region}.tsv.gz"), "wt") as fh:
                kmer_counts(seqs, k).to_csv(fh, sep="\t", index_label="gene_id")
        known = consensus_counts(seqs)
        for meme in a.meme:
            known = known.join(fimo_counts(fasta, meme, a.fimo_p, a.fimo_bin, a.pwm_exclude), how="left")
        known.fillna(0).astype(int).to_csv(os.path.join(a.outdir, f"known.{region}.tsv"), sep="\t", index_label="gene_id")
        print(f"{region}: {len(seqs)} genes, known motifs {known.shape[1]}, "
              f"median sites per gene: {known.median().to_dict()}", flush=True)
    pd.DataFrame(lengths).to_csv(os.path.join(a.outdir, "region_lengths.tsv"), sep="\t", index_label="gene_id")


if __name__ == "__main__":
    main()
