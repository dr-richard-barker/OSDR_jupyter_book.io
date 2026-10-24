"""Download the OSD-766 16S amplicon reads used for the tomato genus table.

Selects every 16S R1/R2 FASTQ from a per-plant root, adventitious-root or leaf sample
(OSD-766 also has fruit, soil, wick and swab samples; those are not used here) and saves
them under the target folder. Re-running skips files that are already complete.

Usage:  python download_reads.py <target_dir>
"""
import concurrent.futures as cf
import json
import os
import re
import sys
import urllib.request

FILES_API = "https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-766/files/"
DOWNLOAD = "https://osdr.nasa.gov/geode-py/ws/studies/OSD-766/download?source=datamanager&file="
# per-plant samples only (SN##), so water-locker "root-mat" and pooled samples are excluded;
# the lookahead is needed because "\b" does not split "root" from "_S96" ("_" is a word char)
WANTED = re.compile(r"16S-.*SN-?\d+-(root|leaf|AdvRoots?)(?=[-_]).*_raw\.fastq\.gz$")


def list_files():
    d = json.load(urllib.request.urlopen(FILES_API, timeout=120))
    files = d[next(iter(d))].get("files", {})
    return sorted(f for f in files if WANTED.search(f))


def fetch(name, target):
    out = os.path.join(target, name)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return name, "kept"
    tmp = out + ".part"
    with urllib.request.urlopen(DOWNLOAD + name, timeout=600) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    os.replace(tmp, out)
    return name, "downloaded"


def main():
    target = sys.argv[1]
    os.makedirs(target, exist_ok=True)
    names = list_files()
    print(f"{len(names)} files selected")
    with cf.ThreadPoolExecutor(max_workers=6) as pool:
        for name, status in pool.map(lambda n: fetch(n, target), names):
            print(status, name)


if __name__ == "__main__":
    main()
