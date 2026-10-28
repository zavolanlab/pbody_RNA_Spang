"""
BED12 of the transcripts of a GTF (exon features), named by gene_id, for calculate-tin.py; and a table of the
genes with biotype, chromosome and exonic length.
"""

import argparse
import collections


def attr(field, key):
    return field.split(f'{key} "')[1].split('"')[0] if f'{key} "' in field else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gtf", required=True)
    p.add_argument("--bed", required=True)
    p.add_argument("--genes", required=True)
    a = p.parse_args()

    exons, info = collections.defaultdict(list), {}
    with open(a.gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = attr(f[8], "transcript_id")
            exons[tid].append((int(f[3]) - 1, int(f[4])))              # 0-based half-open
            info[tid] = (f[0], f[6], attr(f[8], "gene_id"), attr(f[8], "gene_name"), attr(f[8], "gene_biotype"))

    with open(a.bed, "w") as bed, open(a.genes, "w") as genes:
        genes.write("gene_id\tgene_name\tgene_biotype\tchrom\texonic_length\n")
        for tid, ex in sorted(exons.items(), key=lambda kv: (info[kv[0]][0], min(kv[1]))):
            chrom, strand, gid, gname, biotype = info[tid]
            ex = sorted(ex)
            start, end = ex[0][0], ex[-1][1]
            sizes = ",".join(str(e - s) for s, e in ex) + ","
            starts = ",".join(str(s - start) for s, _ in ex) + ","
            bed.write(f"{chrom}\t{start}\t{end}\t{gid}\t0\t{strand}\t{start}\t{end}\t0\t{len(ex)}\t{sizes}\t{starts}\n")
            genes.write(f"{gid}\t{gname or gid}\t{biotype}\t{chrom}\t{sum(e - s for s, e in ex)}\n")


if __name__ == "__main__":
    main()
