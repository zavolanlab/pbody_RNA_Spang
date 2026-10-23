"""
Per-gene sequence features of yeast mRNAs for the analysis of P-body recruitment.

Regions: CDS (spliced, including the stop codon) from the Ensembl GTF; 5' and 3' UTRs from the
TIF-seq-based annotation (CDS_and_UTRs.gtf: '<tx>___5UTR' / '<tx>___3UTR' records, longest full-ORF
transcripts of Pelechano et al. 2013 in YPD/galactose, so the 5' ends are the most upstream TSSs, not
necessarily the dominant ones).

Features:
- lengths and GC content of 5'UTR, CDS, 3'UTR; GC3 (third codon positions);
- tAI (tRNA adaptation index; dos Reis et al. 2004) from the nuclear tRNA gene copy numbers of the
  GTF (anticodon in the gene name, e.g. tR(ACG)D), with wobble weights s(G:U) = 0.41, s(I:C) = 0.28,
  s(I:A) = 0.9999, s(U:G) = 0.68; A34 is read as inosine; w of codons without any tRNA = geometric
  mean of the others; Met and stop codons excluded;
- 5' end: first nucleotide, pyrimidine fraction of the first 10 nt, 'TOP-like' (C at +1 followed by
  >= 4 pyrimidines, the mammalian 5' TOP definition), length of the 5'-terminal pyrimidine tract;
- minimum free energies (RNAfold, ViennaRNA) of the 5'UTR, of the window -50..+50 nt around the start
  codon, and of the 3'UTR (first 1000 nt), also per nucleotide;
- counts of RNA-binding-protein motifs in the 3'UTR and 5'UTR: Puf3 UGUAAAUA, Puf4 UGUAHAUA,
  Whi3 UGCAU, Ssd1 CNYUCNYU.
"""

import argparse
import itertools
import math
import os
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

GENETIC_STOP = {"TAA", "TAG", "TGA"}
S_WOBBLE = {"G:T": 0.41, "I:C": 0.28, "I:A": 0.9999, "T:G": 0.68}
MOTIFS = {"Puf3": "TGTAAATA", "Puf4": "TGTA[ACT]ATA", "Whi3": "TGCAT", "Ssd1": "C[ACGT][CT]TC[ACGT][CT]T"}


def parse_attrs(s):
    out = {}
    for x in s.strip().strip(";").split(";"):
        x = x.strip()
        if x:
            k, v = x.split(" ", 1)
            out[k] = v.strip('"')
    return out


def read_gtf(path, features):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in features:
                continue
            a = parse_attrs(f[8])
            rows.append((f[0], f[2], int(f[3]), int(f[4]), f[6], a.get("gene_id"), a.get("transcript_id"),
                         a.get("gene_name"), a.get("gene_biotype")))
    return pd.DataFrame(rows, columns=["chrom", "feature", "start", "end", "strand", "gene_id", "transcript_id",
                                       "gene_name", "gene_biotype"])


def seq_of(genome, chrom, parts, strand):
    s = "".join(str(genome[chrom].seq[a - 1:b]) for a, b in sorted(parts))
    return str(Seq(s).reverse_complement()) if strand == "-" else s


def gc(s):
    return (s.count("G") + s.count("C")) / len(s) if s else np.nan


def tai_weights(trna_genes):
    """Relative adaptiveness w of the 61 sense codons from tRNA gene copy numbers."""
    copies = {}
    for name in trna_genes:
        m = re.match(r"^t[A-Z]\(([ACGU]{3})\)", name)
        if m:
            ac = m.group(1).replace("U", "T")
            copies[ac] = copies.get(ac, 0) + 1
    W = {}
    for codon in ("".join(p) for p in itertools.product("ACGT", repeat=3)):
        if codon in GENETIC_STOP:
            continue
        total = 0.0
        for ac, n in copies.items():
            wc = str(Seq(ac).reverse_complement())  # Watson-Crick codon of this anticodon
            if wc[:2] != codon[:2]:
                continue
            a34, c3 = ac[0], codon[2]
            if a34 == "G":
                s = {"C": 0.0, "T": S_WOBBLE["G:T"]}.get(c3)
            elif a34 == "A":  # inosine
                s = {"T": 0.0, "C": S_WOBBLE["I:C"], "A": S_WOBBLE["I:A"]}.get(c3)
            elif a34 == "T":
                s = {"A": 0.0, "G": S_WOBBLE["T:G"]}.get(c3)
            else:  # C34
                s = {"G": 0.0}.get(c3)
            if s is not None:
                total += (1.0 - s) * n
        W[codon] = total
    W.pop("ATG", None)
    wmax = max(W.values())
    w = {c: v / wmax for c, v in W.items() if v > 0}
    gmean = math.exp(np.mean([math.log(v) for v in w.values()]))
    return {c: (w[c] if c in w else gmean) for c in W}


def tai(cds, w):
    logs = [math.log(w[cds[i:i + 3]]) for i in range(0, len(cds) - 2, 3) if cds[i:i + 3] in w]
    return math.exp(np.mean(logs)) if logs else np.nan


def rnafold(seqs, threads):
    """MFE (kcal/mol) per sequence with RNAfold; empty or N-containing sequences give NaN."""
    ids = [k for k, s in seqs.items() if s and len(s) >= 8 and "N" not in s]
    out = {k: np.nan for k in seqs}
    if not ids:
        return out
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, "in.fa")
        with open(fa, "w") as fh:
            for k in ids:
                fh.write(f">{k}\n{seqs[k]}\n")
        res = subprocess.run(["RNAfold", "--noPS", f"--jobs={threads}", "-i", fa], capture_output=True, text=True,
                             check=True, cwd=tmp)
    cur = None
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
        elif cur is not None and line.rstrip().endswith(")") and " (" in line:
            out[cur] = float(line.rstrip().rsplit("(", 1)[1].rstrip(")"))
            cur = None
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--utr-gtf", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--genes", required=True, help="genes.tsv of prepare_inputs.py (universe)")
    p.add_argument("--out", required=True)
    p.add_argument("--fasta-dir", required=True)
    p.add_argument("--threads", type=int, default=8)
    a = p.parse_args()
    os.makedirs(a.fasta_dir, exist_ok=True)

    genome = SeqIO.to_dict(SeqIO.parse(a.genome, "fasta"))
    universe = pd.read_csv(a.genes, sep="\t", index_col=0)
    ens = read_gtf(a.gtf, {"CDS", "stop_codon", "gene"})
    trnas = ens[(ens["feature"] == "gene") & (ens["gene_biotype"] == "tRNA") & (ens["chrom"] != "Mito")]["gene_id"]
    w = tai_weights(trnas)

    cds_parts = ens[ens["feature"].isin(["CDS", "stop_codon"]) & ens["gene_id"].isin(universe.index)]
    utr = read_gtf(a.utr_gtf, {"exon"})
    utr = utr[utr["gene_id"].isin(universe.index) & utr["transcript_id"].str.contains("___")]
    utr["region"] = utr["transcript_id"].str.split("___").str[-1]

    rows, fasta = [], {"5UTR": {}, "CDS": {}, "3UTR": {}, "start_window": {}}
    utr_by_gene = {g: d for g, d in utr.groupby("gene_id")}
    empty = utr.iloc[:0]
    for g, sub in cds_parts.groupby("gene_id"):
        chrom, strand = sub["chrom"].iloc[0], sub["strand"].iloc[0]
        cds = seq_of(genome, chrom, list(zip(sub["start"], sub["end"])), strand)
        u = utr_by_gene.get(g, empty)
        u5 = u[u["region"] == "5UTR"]
        u3 = u[u["region"] == "3UTR"]
        s5 = seq_of(genome, chrom, list(zip(u5["start"], u5["end"])), strand) if len(u5) else ""
        s3 = seq_of(genome, chrom, list(zip(u3["start"], u3["end"])), strand) if len(u3) else ""
        cds_ok = len(cds) % 3 == 0 and cds[:3] == "ATG"
        pyr = 0
        for ch in s5:
            if ch in "CT":
                pyr += 1
            else:
                break
        rows.append(dict(
            gene_id=g, len_5utr=len(s5) if len(u5) else np.nan, len_cds=len(cds), len_3utr=len(s3) if len(u3) else np.nan,
            gc_5utr=gc(s5), gc_cds=gc(cds), gc_3utr=gc(s3),
            gc3=gc("".join(cds[i + 2] for i in range(0, len(cds) - 2, 3))) if cds_ok else np.nan,
            tai=tai(cds[:-3], w) if cds_ok else np.nan, cds_valid=cds_ok,
            first_nt_5utr=s5[0] if s5 else "", pyr_frac_first10=(sum(c in "CT" for c in s5[:10]) / min(10, len(s5))) if s5 else np.nan,
            pyr_tract_5end=pyr if s5 else np.nan,
            top_like=bool(s5) and s5[0] == "C" and all(c in "CT" for c in s5[1:5]) and len(s5) >= 5,
            **{f"n_{k}_3utr": len(re.findall(f"(?={v})", s3)) for k, v in MOTIFS.items()},
            **{f"n_{k}_5utr": len(re.findall(f"(?={v})", s5)) for k, v in MOTIFS.items()},
        ))
        if s5:
            fasta["5UTR"][g] = s5
        if s3:
            fasta["3UTR"][g] = s3
        fasta["CDS"][g] = cds
        win5 = s5[-50:] if s5 else ""
        fasta["start_window"][g] = win5 + cds[:50]
    feat = pd.DataFrame(rows).set_index("gene_id")

    for region, key, maxlen in (("5UTR", "mfe_5utr", None), ("start_window", "mfe_start_window", None),
                                ("3UTR", "mfe_3utr", 1000)):
        seqs = {g: (s[:maxlen] if maxlen else s).replace("T", "U") for g, s in fasta[region].items()}
        mfe = rnafold(seqs, a.threads)
        feat[key] = pd.Series(mfe)
        lens = pd.Series({g: len(s) for g, s in seqs.items()})
        feat[key + "_per_nt"] = feat[key] / lens
    for region, seqs in fasta.items():
        with open(os.path.join(a.fasta_dir, f"{region}.fa"), "w") as fh:
            for g, s in seqs.items():
                fh.write(f">{g}\n{s}\n")
    feat.insert(0, "gene_name", universe["gene_name"].reindex(feat.index))
    feat.to_csv(a.out, sep="\t")
    print(feat.describe().T.to_string())
    print("tAI weights:", {k: round(v, 3) for k, v in sorted(w.items())})


if __name__ == "__main__":
    main()
