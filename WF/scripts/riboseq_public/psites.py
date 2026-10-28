"""
P-sites of ribosome footprints (mode 'ribo') or read positions of RNA-seq fragments (mode 'rna') from a
coordinate-sorted, uniquely mapped BAM.

ribo: per read length L, the P-site offset is the distance from the read 5' end to the first nucleotide
      of annotated start codons that is most frequent among reads whose 5' end lies 6-18 nt upstream of a
      start codon (P-site on the AUG of initiating ribosomes). Lengths with fewer than --min-reads such reads,
      or whose most frequent offset holds less than --min-peak-frac of them, are not used. P-sites:
      read 5' end + offset (strand aware). Outputs P-site bigWigs per strand in reads per million P-sites
      (RPM, for IGV), offsets with the 3-nt periodicity of the P-sites in CDSs (fraction in frame 0), and
      P-site counts per gene in the CDS without its first 15 and last 5 codons (initiation and termination
      peaks of cycloheximide-treated samples).
rna:  position = read midpoint; counts per gene in the same trimmed CDS, so that Ribo/RNA ratios refer
      to the same region. Library strandedness (sense or antisense to the CDS) is detected from the reads.
"""

import argparse
import collections
import os

import numpy as np
import pandas as pd
import pyBigWig
import pysam


def parse_gtf_cds(gtf):
    cds = collections.defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "stop_codon"):
                continue
            gid = f[8].split('gene_id "')[1].split('"')[0]
            cds[gid].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))   # 0-based half-open
    return cds


def cds_position_index(cds, chrom_sizes, trim5=45, trim3=15):
    # per chrom and strand: gene index (or -1) and frame of every CDS nucleotide, after trimming the ends
    idx = {(c, s): np.full(n, -1, dtype=np.int32) for c, n in chrom_sizes.items() for s in "+-"}
    frame = {(c, s): np.full(n, -1, dtype=np.int8) for c, n in chrom_sizes.items() for s in "+-"}
    genes, lengths = [], []
    for gid, parts in cds.items():
        chrom, strand = parts[0][0], parts[0][3]
        if (chrom, strand) not in idx:
            continue
        pos = np.concatenate([np.arange(s, e) for c, s, e, _ in sorted(parts, key=lambda p: p[1])])
        pos = np.unique(pos)
        if strand == "-":
            pos = pos[::-1]                                          # transcript order
        n = len(pos)
        fr = np.arange(n) % 3
        keep = slice(trim5, max(trim5, n - trim3))
        if n - trim5 - trim3 < 30:
            continue
        gi = len(genes)
        genes.append(gid)
        lengths.append(n - trim5 - trim3)
        a = idx[(chrom, strand)]
        free = a[pos[keep]] == -1                                    # overlapping CDSs: first gene wins
        a[pos[keep][free]] = gi
        frame[(chrom, strand)][pos[keep]] = fr[keep]
    return idx, frame, genes, np.array(lengths)


def read_starts(bed):
    plus, minus = collections.defaultdict(set), collections.defaultdict(set)
    for line in open(bed):
        f = line.split("\t")
        c, s, e, strand = f[0], int(f[1]), int(f[2]), f[5].strip()
        (plus if strand == "+" else minus)[c].add(s if strand == "+" else e - 1)
    return plus, minus


def five_prime(read):
    return read.reference_end - 1 if read.is_reverse else read.reference_start


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--mode", choices=["ribo", "rna"], required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--chrom-sizes", required=True, help="genome .fai")
    p.add_argument("--start-codons", help="BED of start codons (ribo)")
    p.add_argument("--lengths", default="25-35")
    p.add_argument("--min-reads", type=int, default=200)
    p.add_argument("--min-peak-frac", type=float, default=0.3)
    p.add_argument("--out-prefix", required=True)
    a = p.parse_args()
    os.makedirs(os.path.dirname(a.out_prefix), exist_ok=True)

    sizes = dict(l.split("\t")[:2] for l in open(a.chrom_sizes))
    sizes = {c: int(n) for c, n in sizes.items()}
    gidx, gframe, genes, glen = cds_position_index(parse_gtf_cds(a.gtf), sizes)
    bam = pysam.AlignmentFile(a.bam)
    lo, hi = (int(x) for x in a.lengths.split("-"))

    if a.mode == "ribo":
        plus, minus = read_starts(a.start_codons)
        hist = collections.defaultdict(collections.Counter)
        for r in bam.fetch(until_eof=True):
            L = r.query_alignment_length
            if r.is_unmapped or r.is_secondary or not lo <= L <= hi:
                continue
            p5 = five_prime(r)
            for d in range(6, 19):
                if (not r.is_reverse and p5 + d in plus[r.reference_name]) or \
                   (r.is_reverse and p5 - d in minus[r.reference_name]):
                    hist[L][d] += 1
        rows = []
        for L in range(lo, hi + 1):
            n = sum(hist[L].values())
            d, k = hist[L].most_common(1)[0] if n else (np.nan, 0)
            rows.append({"length": L, "n_reads_near_start": n, "offset": d, "peak_frac": k / n if n else np.nan,
                         "used": bool(n >= a.min_reads and k / n >= a.min_peak_frac)})
        off = pd.DataFrame(rows)
        offsets = {int(r.length): int(r.offset) for r in off.itertuples() if r.used}
        if not offsets:
            raise RuntimeError("no read length with a clear P-site offset")
    else:
        offsets = {}

    cov = {(c, s): np.zeros(n, dtype=np.float32) for c, n in sizes.items() for s in "+-"}
    frame_counts = collections.defaultdict(lambda: np.zeros(3))
    bam = pysam.AlignmentFile(a.bam)   # second pass over the file
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary:
            continue
        s = "-" if r.is_reverse else "+"
        if a.mode == "ribo":
            L = r.query_alignment_length
            if L not in offsets:
                continue
            p5 = five_prime(r)
            pos = p5 - offsets[L] if r.is_reverse else p5 + offsets[L]
        else:
            pos = (r.reference_start + r.reference_end - 1) // 2
        c = r.reference_name
        if 0 <= pos < sizes[c]:
            cov[(c, s)][pos] += 1
            if a.mode == "ribo":
                f = gframe[(c, s)][pos]
                if f >= 0:
                    frame_counts[r.query_alignment_length][f] += 1

    # per-gene counts; RNA: sense or antisense library, whichever gives more CDS reads
    def gene_counts(flip):
        cnt = np.zeros(len(genes))
        for (c, s), v in cov.items():
            g = gidx[(c, {"+": "-", "-": "+"}[s] if flip else s)]
            m = (g >= 0) & (v > 0)
            np.add.at(cnt, g[m], v[m])
        return cnt
    counts, strand_note = gene_counts(False), "sense"
    if a.mode == "rna":
        anti = gene_counts(True)
        if anti.sum() > counts.sum():
            counts, strand_note = anti, "antisense"
        print(f"RNA library: {strand_note} (CDS reads sense {gene_counts(False).sum():.0f}, antisense {anti.sum():.0f})")
    pd.DataFrame({"gene_id": genes, "count": counts.astype(int), "cds_length_used": glen}).to_csv(
        a.out_prefix + ".cds_counts.tsv", sep="\t", index=False)

    if a.mode == "ribo":
        off["frame0_frac"] = off["length"].map(lambda L: frame_counts[L][0] / frame_counts[L].sum()
                                                if frame_counts[L].sum() else np.nan)
        off.to_csv(a.out_prefix + ".psite_offsets.tsv", sep="\t", index=False)
        total = sum(v.sum() for v in cov.values())
        for s, name in (("+", "plus"), ("-", "minus")):
            bw = pyBigWig.open(f"{a.out_prefix}.psites.{name}.rpm.bw", "w")
            bw.addHeader(sorted(sizes.items()))
            for c in sorted(sizes):
                v = cov[(c, s)]
                nz = np.nonzero(v)[0]
                if len(nz):
                    bw.addEntries([c] * len(nz), nz.tolist(), ends=(nz + 1).tolist(),
                                  values=(v[nz] * 1e6 / total).astype(float).tolist())
            bw.close()
        print(f"offsets {offsets}; P-sites {total:.0f}, in CDS {counts.sum():.0f}; frame0 "
              f"{off.loc[off['used'], 'frame0_frac'].round(2).tolist()}")


if __name__ == "__main__":
    main()
