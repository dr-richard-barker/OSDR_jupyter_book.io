"""Inspect OSD-766 16S reads before DADA2: which primer they start with, read lengths,
and where quality drops, so truncation lengths are chosen from the data.

Usage:  python inspect_reads.py <reads_dir> [n_files]
"""
import collections
import gzip
import itertools
import os
import re
import statistics
import sys

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[GC]",
         "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
         "V": "[ACG]", "N": "[ACGT]"}
PRIMERS = {  # widely used 16S primers, checked against the read starts
    "515F (V4)": "GTGYCAGCMGCCGCGGTAA", "806R (V4)": "GGACTACNVGGGTWTCTAAT",
    "341F (V3-V4)": "CCTACGGGNGGCWGCAG", "785R (V3-V4)": "GACTACHVGGGTATCTAATCC",
    "805R (V3-V4)": "GACTACHVGGGTATCTAATCC", "799F (V5-V7)": "AACMGGATTAGATACCCKG",
    "1193R (V5-V7)": "ACGTCATCCCCACCTTCC", "27F (V1-V3)": "AGAGTTTGATCMTGGCTCAG",
    "534R (V1-V3)": "ATTACCGCGGCTGCTGG", "338F": "ACTCCTACGGGAGGCAGCAG",
}


def rx(seq):
    return re.compile("".join(IUPAC[b] for b in seq))


def reads(path, n):
    with gzip.open(path, "rt") as fh:
        for rec in itertools.islice(zip(*[fh] * 4), n):
            yield rec[1].strip(), rec[3].strip()


def main():
    d = sys.argv[1]
    n_files = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    r1 = sorted(f for f in os.listdir(d) if f.endswith("_R1_raw.fastq.gz"))[:n_files]
    for tag in ("R1", "R2"):
        hits, starts, lens = collections.Counter(), collections.Counter(), []
        qpos = collections.defaultdict(list)
        n = 0
        for f in r1:
            for seq, qual in reads(os.path.join(d, f.replace("_R1_", f"_{tag}_")), 4000):
                n += 1
                lens.append(len(seq))
                starts[seq[:20]] += 1
                for name, p in PRIMERS.items():
                    if rx(p).match(seq):
                        hits[name] += 1
                if n % 4 == 0:
                    for i in range(0, len(qual), 10):
                        qpos[i].append(ord(qual[i]) - 33)
        print(f"== {tag}: {n} reads from {len(r1)} samples; length median {statistics.median(lens)} "
              f"(min {min(lens)}, max {max(lens)})")
        print("   primer matches at read start:",
              {k: f"{100 * v / n:.1f}%" for k, v in hits.most_common()} or "none")
        print("   commonest 20-mers at read start:", starts.most_common(3))
        print("   median quality by position:",
              " ".join(f"{i}:{statistics.median(v):.0f}" for i, v in sorted(qpos.items())))


if __name__ == "__main__":
    main()
