"""Turn the DADA2 genus counts into the tables the pipeline chapter reads.

Writes
  samples.csv                       one row per 16S sample: OSDR spaceflight + light factors,
                                    tissue, read depth, and whether it passes the depth filter
  ../tomato_genus_by_sample.csv     relative abundance (% of bacterial reads), genera x samples,
                                    for samples with >= MIN_READS bacterial reads

Flight/ground and light come from OSDR's own metadata, not from parsing names. Leaves are
nearly all chloroplast reads (about 99%), so most fail the depth filter.

Usage:  python make_genus_table.py      (run inside data_tomato_16s/)
"""
import io
import re

import pandas as pd
import requests

MIN_READS = 1000
API = "https://visualization.osdr.nasa.gov/biodata/api/v2"


def tissue(name):
    n = name.lower()
    return "leaf" if "leaf" in n else ("adventitious root" if "advroot" in n else "root")


def main():
    counts = pd.read_csv("genus_counts.csv", index_col=0)
    track = pd.read_csv("read_tracking.csv").set_index("sample")
    q = {"id.accession": "OSD-766", "id.sample name": "", "study.factor value": "", "format": "csv"}
    meta = (pd.read_csv(io.StringIO(requests.get(f"{API}/query/metadata/", params=q, timeout=120).text))
            .drop_duplicates("id.sample name").set_index("id.sample name"))
    rows = []
    for s in counts.columns:
        plant = re.search(r"SN-?0*(\d+)", s)
        rows.append({
            "sample": s,
            "condition": {"Space Flight": "Flight", "Ground Control": "Ground"}[meta.loc[s, "study.factor value.spaceflight"]],
            "light": meta.loc[s, "study.factor value.light regimen"],
            "tissue": tissue(s),
            "plant": int(plant.group(1)) if plant else None,
            "input_reads": int(track.loc[s, "input"]),
            "organelle_reads": int(track.loc[s, "organelle"]),
            "bacterial_reads": int(track.loc[s, "bacterial"]),
        })
    samples = pd.DataFrame(rows)
    samples["kept"] = samples["bacterial_reads"] >= MIN_READS
    samples.to_csv("samples.csv", index=False)

    keep = samples.loc[samples["kept"], "sample"]
    rel = counts[keep] / counts[keep].sum() * 100
    rel.round(4).to_csv("../tomato_genus_by_sample.csv")
    print(f"{len(keep)} of {len(samples)} samples pass >= {MIN_READS} bacterial reads; "
          f"wrote ../tomato_genus_by_sample.csv ({rel.shape[0]} genera)")
    print(samples.groupby(["condition", "tissue"])["kept"].agg(["size", "sum"]).to_string())


if __name__ == "__main__":
    main()
