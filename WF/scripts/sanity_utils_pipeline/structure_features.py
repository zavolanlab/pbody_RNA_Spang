"""
Per-gene mRNA secondary-structure metrics from yeast DMS-seq (Rouskin et al. 2014), in R64-1-1 coordinates.

Data
----
Rouskin S, Zubradt M, Washietl S, Kellis M, Weissman JS (2014). Genome-wide probing of RNA structure
reveals active unfolding of mRNA structures in vivo. Nature 505:701-705. doi:10.1038/nature12894.
GEO GSE45803. The processed yeast files are per-strand variableStep WIGs (span=1) of reverse-transcriptase
stop counts; "the position number directly indicates the DMS modified position starting at 1", replicates
pooled, genome build SGD R61/R62 (UCSC sacCer2). Conditions used here:

    in_vivo      GSM1115749-52  Feb13_VivoAllextra_1_15_{PLUS,Minus}.wig   (BY4741, log phase, DMS in culture)
    in_vitro     GSM1115756-61  Feb13_VitroAll_1_15_{Plus,Minus}.wig       (polyA RNA refolded at 30 C, then DMS)
    denatured    GSM1115753-55  Feb13_DenaturedAll_1_15_{Plus,MINUS}.wig   (DMS at 95 C; 'unstructured' control)
    atp_depleted GSM1297509     Yeast_noATP_{plus,minus}.wig               (S288C, ATP-depleted, DMS in vivo)

No untreated (no-DMS) control is provided as processed data in GSE45803.

Two sub-commands
----------------
1. ``lift``      WIG (sacCer2) -> per-strand bedGraph.gz in R64-1-1 / sacCer3 coordinates with Ensembl
                 chromosome names (chrI -> I, chrM/chrMito -> Mito), using UCSC ``liftOver`` and
                 ``sacCer2ToSacCer3.over.chain.gz``. Writes a QC table with the fraction of reads (and of
                 positions) that fall on A/C of the transcribed strand in (i) the source genome, (ii) the
                 target genome *before* liftover and (iii) the target genome *after* liftover.
2. ``features``  per gene x condition x region metrics (below) -> long TSV, plus a companion table of
                 Gini differences between conditions for gene/regions passing coverage in both.

Regions (built in transcript 5'->3' order from the TIF-seq based GTF; exon pieces carry transcript_id
suffixes ___5UTR, ___CDS (includes the stop codon), ___3UTR; CDS introns are spliced out):

    5UTR, CDS, 3UTR      the annotated pieces (a row is emitted only if the piece is annotated)
    start_region         5'UTR + first ``--start-cds-nt`` (default 50) nt of the CDS; only if a 5'UTR exists
    mRNA                 5'UTR + CDS + 3'UTR as annotated (= CDS only for genes without annotated UTRs)

Only A and C positions (base of the transcribed strand in the R64-1-1 genome) are used, since DMS
methylates unpaired A (N1) and C (N3) and only those are read out as RT stops.

Mappability mask (``--mappability-k``, default 25; 0 disables). Rouskin aligned the first 25 nt of each
read with bowtie and kept only unique, perfect matches, so RT stops whose read k-mer occurs more than
once in the genome are invisible (counts 0 irrespective of structure). This affects mainly paralogous
gene pairs, e.g. ribosomal-protein genes: only ~14% of the RPL1A/RPL1B CDS positions are unique and their
CDS Gini drops from ~0.85 to ~0.32-0.37 once masked. An A/C position p is kept only if the k-mer of the read
it generates is unique in the R64-1-1 genome (both strands, no N): genomic [p+1, p+k] for '+' strand
transcripts and [p-k, p-1] for '-' strand transcripts (the read starts 1 nt 3' of the modified base in the
RNA; this geometry gave the lowest fraction of in vivo reads at non-unique k-mers, 2.3% vs 5.6% for the
opposite orientation). n_AC therefore counts mappable A/C positions.

Metrics, for a region with A/C raw counts x_1..x_n (n = n_AC):

    mean_reads_per_AC     = (1/n) * sum_i x_i                       (raw RT-stop counts; coverage)
    passes_coverage       = mean_reads_per_AC >= --min-reads (15). Rouskin 2014: "Raw data was filtered for
                            messages that have at least 15 reads on average per A or C position"; here the
                            criterion is applied to each region separately.
    gini                  Gini index over all A/C positions of the region. With x sorted ascending,
                              G = sum_{i=1..n} (2i - n - 1) x_(i) / (n * sum_i x_i)
                          (0 = uniform reactivity, i.e. unstructured; -> 1 = signal concentrated in a few
                          positions, i.e. structured/protected stretches). G is scale invariant, so it is
                          identical for raw counts and for any per-gene linear normalisation (Rouskin
                          normalised each window to its most reactive base); it is computed on raw counts.
    mean_window_gini      mean of G over consecutive non-overlapping windows of --window (50) A/C positions,
                          tiled from the 5' end of the region; the trailing incomplete window is dropped
                          and windows without any reads are skipped (NaN if the region has < 50 A/C).
                          Rouskin 2014 computed Gini (and r vs denatured) on "windows containing a total of
                          50 A/C nucleotides" and called a window structured at Gini > 0.14 (with r < 0.55),
                          unstructured at Gini < 0.08 (with r > 0.70). Note that Rouskin tiled windows across
                          the whole message from the 5'UTR; here windows are tiled within each region.
    mean_norm_reactivity  mean over the region's A/C of the per-gene normalised reactivity
                              r_i = min(x_i, P95) / P95,
                          P95 = 95th percentile of the A/C counts over the full mRNA of that gene in that
                          condition (i.e. 90% winsorisation + scaling, as Rouskin 2014 did for rRNA; if
                          P95 = 0 the maximum is used). Range 0..1; higher = more accessible on average
                          relative to the gene's own reactive positions.

Difference table (``--diff-out``): for each comparison A:B (default in_vivo:in_vitro, in_vivo:denatured,
in_vitro:denatured, atp_depleted:in_vivo, when present) and gene/region passing coverage in both
conditions, delta_gini = gini_A - gini_B and delta_mean_window_gini likewise. Rouskin 2014 report mRNAs
to be much less structured in vivo than in vitro (delta < 0 for in_vivo - in_vitro) and ATP depletion to
increase structure (delta > 0 for atp_depleted - in_vivo). Because DMS reactivity of an unstructured
sequence is itself uneven (A >> C, sequence context), the absolute Gini of the denatured control is far
from 0 (~0.3 for CDS windows); X - denatured is the denatured-adjusted structure of condition X.

Caveats: the absolute Gini increases with low coverage (Poisson noise) and with stretches without signal
(e.g. annotated UTR parts upstream/downstream of the ends of the dominant isoforms in the DMS-seq
samples), so cross-gene comparisons of UTR / start_region Gini should control for mean_reads_per_AC and
n_AC or use the within-gene differences to the denatured control. With the TIF-seq GTF, the DMS signal
of the log-phase samples is essentially confined to ~50-80 nt upstream of the AUG and ~100 nt downstream
of the stop codon, while the annotated UTRs are longer (median 5'UTR 117 nt, 3'UTR 196 nt); the optional
--max-5utr-nt / --max-3utr-nt caps restrict the UTRs to the codon-proximal part. Positions shared by overlapping genes
on the same strand count for each gene; RT-stop counts are affected by fragmentation/polyA selection
biases (Rouskin normalised within 50-200 nt windows for that reason; window Gini is robust to it,
whole-region Gini and mean_norm_reactivity less so); the WIGs of in_vivo/in_vitro/denatured contain no
mitochondrial chromosome, so Mito genes are skipped for those conditions.

Run that produced external_data/Rouskin2014_DMSseq/structure_features.R64.tsv
---------------------------------------------------------------------------
(python = sanity_utils env, run with PYTHONNOUSERSITE=1; liftOver from the py312 env; the sacCer2 FASTA,
UCSC bigZips/chromFa.tar.gz, is optional and only used for the native A/C rate in the QC table)

    D=/scicore/home/zavolan/GROUP/pbody_RNA_spang_zavolan_v2/external_data/Rouskin2014_DMSseq
    R=$D/raw; L=$D/lifted_R64
    GENOME=/scicore/home/zavolan/GROUP/Genomes/saccharomyces_cerevisiae/S288C_reference_genome_R64-1-1/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
    GTF=/scicore/home/zavolan/GROUP/pbody_RNA_spang_zavolan_v2/wf_runs/v1/modified_annotation/CDS_and_UTRs.gtf

    python structure_features.py lift \\
        --wig in_vivo:+:$R/GSE45803_Feb13_VivoAllextra_1_15_PLUS.wig.gz \\
        --wig in_vivo:-:$R/GSE45803_Feb13_VivoAllextra_1_15_Minus.wig.gz \\
        --wig in_vitro:+:$R/GSE45803_Feb13_VitroAll_1_15_Plus.wig.gz \\
        --wig in_vitro:-:$R/GSE45803_Feb13_VitroAll_1_15_Minus.wig.gz \\
        --wig denatured:+:$R/GSE45803_Feb13_DenaturedAll_1_15_Plus.wig.gz \\
        --wig denatured:-:$R/GSE45803_Feb13_DenaturedAll_1_15_MINUS.wig.gz \\
        --wig atp_depleted:+:$R/GSM1297509_Yeast_noATP_plus.wig.gz \\
        --wig atp_depleted:-:$R/GSM1297509_Yeast_noATP_minus.wig.gz \\
        --chain $D/liftover/sacCer2ToSacCer3.over.chain.gz \\
        --liftover-bin /scicore/home/zavolan/mirono0000/miniforge3_Dec2024/envs/py312/bin/liftOver \\
        --target-genome $GENOME [--source-genome sacCer2.fa] --outdir $L

    python structure_features.py features \\
        --track in_vivo:$L/in_vivo.plus.R64.bedGraph.gz:$L/in_vivo.minus.R64.bedGraph.gz \\
        --track in_vitro:$L/in_vitro.plus.R64.bedGraph.gz:$L/in_vitro.minus.R64.bedGraph.gz \\
        --track denatured:$L/denatured.plus.R64.bedGraph.gz:$L/denatured.minus.R64.bedGraph.gz \\
        --track atp_depleted:$L/atp_depleted.plus.R64.bedGraph.gz:$L/atp_depleted.minus.R64.bedGraph.gz \\
        --gtf $GTF --genome $GENOME --out $D/structure_features.R64.tsv
"""

import argparse
import collections
import gzip
import os
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam

REGIONS = ("5UTR", "CDS", "3UTR", "start_region", "mRNA")
DEFAULT_COMPARISONS = ("in_vivo:in_vitro", "in_vivo:denatured", "in_vitro:denatured", "atp_depleted:in_vivo")


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def to_ensembl_chrom(chrom):
    """UCSC/SGD-style chromosome name -> Ensembl R64-1-1 name (chrIV -> IV, chrM/chrMito -> Mito)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return "Mito" if c in ("M", "Mito", "mt", "MT") else c


def to_ucsc_chrom(chrom):
    """Chromosome name as used in the UCSC sacCer2 chain (chrMito -> chrM)."""
    c = to_ensembl_chrom(chrom)
    return "chrM" if c == "Mito" else "chr" + c


# ----------------------------------------------------------------------------------------------------
# lift
# ----------------------------------------------------------------------------------------------------

def read_wig(path):
    """variableStep WIG -> DataFrame(chrom, pos1, value). fixedStep is not used by GSE45803."""
    rows, chrom, span = [], None, 1
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(("track", "#", "browser")) or not line.strip():
                continue
            if line.startswith("variableStep"):
                kv = dict(t.split("=") for t in line.split()[1:])
                chrom, span = kv["chrom"], int(kv.get("span", 1))
                if span != 1:
                    raise ValueError(f"{path}: span={span} not supported")
                continue
            if line.startswith("fixedStep"):
                raise ValueError(f"{path}: fixedStep not supported")
            p, v = line.split()
            rows.append((chrom, int(p), float(v)))
    return pd.DataFrame(rows, columns=["chrom", "pos1", "value"])


def ac_fraction(df, strand, fasta, chrom_fn):
    """Fraction of reads and of positions on A/C of the given strand. df: chrom, pos1 (1-based), value."""
    reads_ac = reads = pos_ac = npos = 0.0
    for chrom, sub in df.groupby("chrom", sort=False):
        c = chrom_fn(chrom)
        if c not in fasta.references:
            continue
        seq = np.frombuffer(fasta.fetch(c).upper().encode(), dtype="S1")
        q = sub["pos1"].to_numpy() - 1
        ok = (q >= 0) & (q < seq.size)
        b = seq[q[ok]]
        is_ac = np.isin(b, [b"A", b"C"]) if strand == "+" else np.isin(b, [b"T", b"G"])
        v = sub["value"].to_numpy()[ok]
        reads_ac += v[is_ac].sum()
        reads += v.sum()
        pos_ac += is_ac.sum()
        npos += ok.sum()
    return (reads_ac / reads if reads else np.nan, pos_ac / npos if npos else np.nan)


def cmd_lift(a):
    os.makedirs(a.outdir, exist_ok=True)
    target = pysam.FastaFile(a.target_genome)
    source = pysam.FastaFile(a.source_genome) if a.source_genome else None
    qc = []
    for spec in a.wig:
        cond, strand, path = spec.split(":", 2)
        if strand not in "+-":
            raise ValueError(f"strand must be + or -: {spec}")
        sname = "plus" if strand == "+" else "minus"
        log(f"[lift] {cond} {strand} {path}")
        wig = read_wig(path)
        wig = wig.groupby(["chrom", "pos1"], as_index=False, sort=False)["value"].sum()

        with tempfile.TemporaryDirectory(dir=a.tmpdir) as tmp:
            bed_in, bed_out, bed_un = (os.path.join(tmp, f) for f in ("in.bed", "out.bed", "unmapped.bed"))
            pd.DataFrame({"chrom": wig["chrom"].map(to_ucsc_chrom), "start": wig["pos1"] - 1,
                          "end": wig["pos1"], "name": wig["value"]}) \
                .to_csv(bed_in, sep="\t", header=False, index=False)
            subprocess.run([a.liftover_bin, bed_in, a.chain, bed_out, bed_un], check=True)
            names = ["chrom", "start", "end", "value"]
            lifted = pd.read_csv(bed_out, sep="\t", header=None, names=names)
            un_lines = [l.rstrip("\n").split("\t") for l in open(bed_un) if l.strip() and not l.startswith("#")]
            unmapped = pd.DataFrame(un_lines, columns=names).astype({"value": float}) if un_lines \
                else pd.DataFrame(columns=names, dtype=float)
        lifted["chrom"] = lifted["chrom"].map(to_ensembl_chrom)
        lifted = lifted.groupby(["chrom", "start", "end"], as_index=False)["value"].sum()
        order = {c: i for i, c in enumerate(target.references)}
        lifted = lifted.sort_values(["chrom", "start"], key=lambda s: s.map(order) if s.name == "chrom" else s)
        out = os.path.join(a.outdir, f"{cond}.{sname}.R64.bedGraph.gz")
        lifted.to_csv(out, sep="\t", header=False, index=False, compression="gzip", float_format="%.6g")

        lifted_df = pd.DataFrame({"chrom": lifted["chrom"], "pos1": lifted["end"], "value": lifted["value"]})
        rec = dict(condition=cond, strand=strand, wig=os.path.basename(path), bedgraph=os.path.basename(out),
                   n_positions=len(wig), total_reads=wig["value"].sum(),
                   n_unmapped=len(unmapped), reads_unmapped=unmapped["value"].sum(),
                   chroms=",".join(sorted(wig["chrom"].unique())))
        if source is not None:
            rec["fracReads_AC_source"], rec["fracPos_AC_source"] = ac_fraction(wig, strand, source, to_ucsc_chrom)
        rec["fracReads_AC_target_before"], rec["fracPos_AC_target_before"] = \
            ac_fraction(wig, strand, target, to_ensembl_chrom)
        rec["fracReads_AC_target_after"], rec["fracPos_AC_target_after"] = \
            ac_fraction(lifted_df, strand, target, to_ensembl_chrom)
        log("   " + ", ".join(f"{k}={v:.4f}" for k, v in rec.items() if k.startswith("frac")))
        qc.append(rec)
    qc_out = a.qc_out or os.path.join(a.outdir, "liftover_AC_qc.tsv")
    pd.DataFrame(qc).to_csv(qc_out, sep="\t", index=False, float_format="%.6g")
    log(f"[lift] QC -> {qc_out}")


# ----------------------------------------------------------------------------------------------------
# features
# ----------------------------------------------------------------------------------------------------

def parse_gtf(path):
    """Genes with a ___CDS piece -> {gene_id: dict(chrom, strand, pieces={piece: [(start1, end1), ...]})}."""
    genes = {}
    attr = re.compile(r'(\S+) "([^"]*)"')
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            at = dict(attr.findall(f[8]))
            tid = at.get("transcript_id", "")
            if "___" not in tid:
                continue
            piece = tid.rsplit("___", 1)[1]
            if piece not in ("5UTR", "CDS", "3UTR"):
                continue
            g = genes.setdefault(at["gene_id"], dict(chrom=f[0], strand=f[6], pieces=collections.defaultdict(list)))
            if g["chrom"] != f[0] or g["strand"] != f[6]:
                raise ValueError(f"inconsistent chrom/strand for {at['gene_id']}")
            g["pieces"][piece].append((int(f[3]), int(f[4])))
    return {k: v for k, v in genes.items() if "CDS" in v["pieces"]}


def unique_kmers(fasta, k):
    """{chrom: bool array u} with u[i] True iff genome[i:i+k] (no N) occurs exactly once in the genome,
    counting both strands. 2-bit encoded k-mers (k <= 32) in uint64; ~1 GB peak memory for yeast."""
    if not 0 < k <= 32:
        raise ValueError("mappability k must be in 1..32")
    lut = np.full(256, 255, dtype=np.uint8)
    for i, b in enumerate(b"ACGT"):
        lut[b] = i
    fwd, rev, bad = {}, {}, {}
    for c in fasta.references:
        s = lut[np.frombuffer(fasta.fetch(c).upper().encode(), dtype=np.uint8)]
        n = s.size - k + 1
        if n <= 0:
            fwd[c] = rev[c] = np.zeros(0, np.uint64)
            bad[c] = np.zeros(0, bool)
            continue
        bad[c] = np.convolve((s == 255).astype(np.int32), np.ones(k, np.int32), "valid") > 0
        s = np.where(s == 255, 0, s).astype(np.uint64)
        f, r = np.zeros(n, np.uint64), np.zeros(n, np.uint64)
        for j in range(k):
            f = (f << np.uint64(2)) | s[j:j + n]
            r |= (np.uint64(3) - s[j:j + n]) << np.uint64(2 * j)  # reverse complement
        fwd[c], rev[c] = f, r
    allk = np.concatenate([x[~bad[c]] for c in fasta.references for x in (fwd[c], rev[c])])
    u, cnt = np.unique(allk, return_counts=True)
    del allk
    return {c: (cnt[np.searchsorted(u, fwd[c])] == 1) & ~bad[c] if fwd[c].size else np.zeros(0, bool)
            for c in fasta.references}


def build_transcripts(genes, fasta, start_cds_nt, uniq=None, k=0, max_5utr=None, max_3utr=None):
    """For each gene: 0-based genomic positions of (mappable) A/C nucleotides in transcript 5'->3' order and,
    per region, the indices into that A/C array. uniq/k: see unique_kmers (None = no mappability mask).
    max_5utr / max_3utr: keep only the 5'UTR nt next to the start codon / 3'UTR nt next to the stop codon."""
    seqs = {c: np.frombuffer(fasta.fetch(c).upper().encode(), dtype="S1") for c in fasta.references}
    out = {}
    for gid, g in genes.items():
        if g["chrom"] not in seqs:
            continue
        minus = g["strand"] == "-"
        pos, lab = [], []
        for li, piece in enumerate(("5UTR", "CDS", "3UTR")):
            for s, e in sorted(g["pieces"].get(piece, []), reverse=minus):
                p = np.arange(s - 1, e)
                pos.append(p[::-1] if minus else p)
                lab.append(np.full(p.size, li, dtype=np.int8))
        pos, lab = np.concatenate(pos), np.concatenate(lab)
        keep = np.ones(pos.size, bool)
        if max_5utr:
            keep[np.flatnonzero(lab == 0)[:-max_5utr]] = False
        if max_3utr:
            keep[np.flatnonzero(lab == 2)[max_3utr:]] = False
        pos, lab = pos[keep], lab[keep]
        base = seqs[g["chrom"]][pos]  # genomic plus-strand base; A/C on the minus strand = T/G on plus
        is_ac = np.isin(base, [b"T", b"G"] if minus else [b"A", b"C"])
        if uniq is not None:  # read k-mer: genomic [p+1, p+k] ('+') or [p-k, p-1] ('-')
            u = uniq[g["chrom"]]
            st = pos - k if minus else pos + 1
            ok = (st >= 0) & (st < u.size)
            mappable = np.zeros(pos.size, bool)
            mappable[ok] = u[st[ok]]
            is_ac &= mappable

        cds_idx = np.flatnonzero(lab == 1)
        in_start = (lab == 0)
        in_start[cds_idx[:start_cds_nt]] = True
        masks = {"5UTR": lab == 0, "CDS": lab == 1, "3UTR": lab == 2, "mRNA": np.ones(lab.size, bool)}
        if "5UTR" in g["pieces"]:
            masks["start_region"] = in_start
        ac_pos = pos[is_ac]
        region_idx = {r: np.flatnonzero(m[is_ac]) for r, m in masks.items()
                      if (r in ("start_region", "mRNA") or r in g["pieces"])}
        out[gid] = dict(chrom=g["chrom"], strand=g["strand"], ac_pos=ac_pos, region_idx=region_idx)
    return out


def load_track(plus_path, minus_path, chrom_sizes):
    """Per-strand bedGraph(.gz) -> {(chrom, strand): float array over the chromosome}."""
    track = {}
    for strand, path in (("+", plus_path), ("-", minus_path)):
        with open_text(path) as fh:  # skip optional 'track'/'browser' header lines
            n_header = 0
            for line in fh:
                if not line.startswith(("track", "browser", "#")):
                    break
                n_header += 1
        df = pd.read_csv(path, sep="\t", header=None, skiprows=n_header,
                         names=["chrom", "start", "end", "value"], dtype={"chrom": str})
        for chrom, sub in df.groupby("chrom", sort=False):
            chrom = to_ensembl_chrom(chrom)
            if chrom not in chrom_sizes:
                continue
            arr = np.zeros(chrom_sizes[chrom], dtype=np.float64)
            st, en, v = sub["start"].to_numpy(), sub["end"].to_numpy(), sub["value"].to_numpy()
            if np.all(en - st == 1):
                np.add.at(arr, st, v)
            else:
                for s, e, x in zip(st, en, v):
                    arr[s:e] += x
            track[(chrom, strand)] = arr
    return track


def gini(x):
    """Gini index of a 1-D non-negative vector; NaN if empty or all zero."""
    x = np.sort(np.asarray(x, dtype=float))
    n, s = x.size, x.sum()
    if n == 0 or s <= 0:
        return np.nan
    i = np.arange(1, n + 1)
    return float(np.sum((2 * i - n - 1) * x) / (n * s))


def window_gini(x, w):
    """Gini of consecutive non-overlapping windows of w values (trailing remainder dropped)."""
    k = x.size // w
    if k == 0:
        return np.array([])
    m = np.sort(x[: k * w].reshape(k, w), axis=1)
    s = m.sum(axis=1)
    i = np.arange(1, w + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        g = (m * (2 * i - w - 1)).sum(axis=1) / (w * s)
    return g[s > 0]


def gene_metrics(counts, region_idx, a):
    """Metrics for one gene x condition. counts: raw counts at the gene's A/C positions (5'->3')."""
    p = np.percentile(counts, a.norm_percentile) if counts.size else 0.0
    denom = p if p > 0 else (counts.max() if counts.size else 0.0)
    norm = np.minimum(counts, denom) / denom if denom > 0 else np.full(counts.size, np.nan)
    rows = []
    for region in REGIONS:
        if region not in region_idx:
            continue
        idx = region_idx[region]
        if idx.size == 0:
            continue
        x = counts[idx]
        wg = window_gini(x, a.window)
        mean_reads = float(x.mean())
        rows.append(dict(region=region, n_AC=int(idx.size), mean_reads_per_AC=mean_reads, gini=gini(x),
                         mean_window_gini=float(wg.mean()) if wg.size else np.nan,
                         mean_norm_reactivity=float(np.nanmean(norm[idx])) if denom > 0 else np.nan,
                         passes_coverage=bool(mean_reads >= a.min_reads)))
    return rows


def cmd_features(a):
    fasta = pysam.FastaFile(a.genome)
    chrom_sizes = dict(zip(fasta.references, fasta.lengths))
    genes = parse_gtf(a.gtf)
    uniq = None
    if a.mappability_k > 0:
        uniq = unique_kmers(fasta, a.mappability_k)
        log(f"[features] {a.mappability_k}-mers unique in the genome: "
            f"{np.mean(np.concatenate(list(uniq.values()))):.4f} of positions")
    tx = build_transcripts(genes, fasta, a.start_cds_nt, uniq, a.mappability_k, a.max_5utr_nt, a.max_3utr_nt)
    log(f"[features] {len(tx)} genes with CDS; "
        f"{sum('5UTR' in g['pieces'] for g in genes.values())} with 5'UTR, "
        f"{sum('3UTR' in g['pieces'] for g in genes.values())} with 3'UTR")

    rows, conditions = [], []
    for spec in a.track:
        cond, plus, minus = spec.split(":", 2)
        conditions.append(cond)
        track = load_track(plus, minus, chrom_sizes)
        present = sorted({c for c, _ in track})
        skipped = 0
        for gid, t in tx.items():
            arr = track.get((t["chrom"], t["strand"]))
            if arr is None:
                skipped += 1
                continue
            for r in gene_metrics(arr[t["ac_pos"]], t["region_idx"], a):
                rows.append(dict(gene_id=gid, condition=cond, **r))
        log(f"[features] {cond}: chromosomes in track {','.join(present)}; "
            f"{skipped} genes skipped (chromosome/strand absent from track)")
        del track

    cols = ["gene_id", "condition", "region", "n_AC", "mean_reads_per_AC", "gini", "mean_window_gini",
            "mean_norm_reactivity", "passes_coverage"]
    res = pd.DataFrame(rows, columns=cols)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    res.to_csv(a.out, sep="\t", index=False, float_format="%.6g", na_rep="NA")
    log(f"[features] {len(res)} rows -> {a.out}")

    comparisons = a.compare or [c for c in DEFAULT_COMPARISONS
                                if all(x in conditions for x in c.split(":"))]
    diff_out = a.diff_out or re.sub(r"(\.tsv)?$", ".gini_diff.tsv", a.out, count=1)
    ok = res[res["passes_coverage"]]
    diffs = []
    for comp in comparisons:
        ca, cb = comp.split(":")
        m = ok[ok["condition"] == ca].merge(ok[ok["condition"] == cb], on=["gene_id", "region"],
                                            suffixes=("_a", "_b"))
        diffs.append(pd.DataFrame({
            "gene_id": m["gene_id"], "region": m["region"], "comparison": f"{ca}_minus_{cb}",
            "n_AC": m["n_AC_a"], "gini_a": m["gini_a"], "gini_b": m["gini_b"],
            "delta_gini": m["gini_a"] - m["gini_b"],
            "mean_window_gini_a": m["mean_window_gini_a"], "mean_window_gini_b": m["mean_window_gini_b"],
            "delta_mean_window_gini": m["mean_window_gini_a"] - m["mean_window_gini_b"]}))
    if diffs:
        d = pd.concat(diffs, ignore_index=True)
        d.to_csv(diff_out, sep="\t", index=False, float_format="%.6g", na_rep="NA")
        log(f"[features] {len(d)} difference rows ({', '.join(comparisons)}) -> {diff_out}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    pl = sub.add_parser("lift", help="WIG (sacCer2) -> per-strand bedGraph.gz (R64-1-1) + A/C QC")
    pl.add_argument("--wig", action="append", required=True, metavar="COND:STRAND:PATH",
                    help="per-strand WIG, STRAND is + or - (repeat)")
    pl.add_argument("--chain", required=True, help="sacCer2ToSacCer3.over.chain.gz")
    pl.add_argument("--liftover-bin", default="liftOver")
    pl.add_argument("--target-genome", required=True, help="R64-1-1 FASTA (Ensembl chromosome names), indexed")
    pl.add_argument("--source-genome", help="optional sacCer2 FASTA (UCSC names) for the native A/C rate")
    pl.add_argument("--outdir", required=True)
    pl.add_argument("--qc-out", help="default: <outdir>/liftover_AC_qc.tsv")
    pl.add_argument("--tmpdir", default=None, help="directory for temporary BED files")
    pl.set_defaults(func=cmd_lift)

    pf = sub.add_parser("features", help="per gene/region structure metrics")
    pf.add_argument("--track", action="append", required=True, metavar="COND:PLUS_BEDGRAPH:MINUS_BEDGRAPH",
                    help="lifted per-strand bedGraph(.gz) of one condition, R64 coordinates (repeat)")
    pf.add_argument("--gtf", required=True, help="GTF with ___5UTR/___CDS/___3UTR transcript_id suffixes")
    pf.add_argument("--genome", required=True, help="R64-1-1 genome FASTA (indexed)")
    pf.add_argument("--out", required=True, help="output TSV (long format)")
    pf.add_argument("--diff-out", help="Gini difference TSV; default: <out>.gini_diff.tsv")
    pf.add_argument("--compare", action="append", metavar="COND_A:COND_B",
                    help=f"difference A - B (repeat); default: {', '.join(DEFAULT_COMPARISONS)} if present")
    pf.add_argument("--window", type=int, default=50, help="window size in A/C positions (default 50)")
    pf.add_argument("--min-reads", type=float, default=15.0, help="min mean reads per A/C (default 15)")
    pf.add_argument("--start-cds-nt", type=int, default=50, help="CDS nt included in start_region (default 50)")
    pf.add_argument("--norm-percentile", type=float, default=95.0,
                    help="per-gene winsorisation/scaling percentile of A/C counts (default 95)")
    pf.add_argument("--max-5utr-nt", type=int, default=None,
                    help="optional: use only the N 5'UTR nt upstream of the start codon (default: full annotation)")
    pf.add_argument("--max-3utr-nt", type=int, default=None,
                    help="optional: use only the N 3'UTR nt downstream of the stop codon (default: full annotation)")
    pf.add_argument("--mappability-k", type=int, default=25,
                    help="drop A/C positions whose read k-mer is not unique in the genome (default 25 = "
                         "Rouskin's bowtie seed length; 0 = no mask)")
    pf.set_defaults(func=cmd_features)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
