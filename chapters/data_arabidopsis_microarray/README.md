# Arabidopsis microarray chapter — data provenance

Data backing `OSDR_arabidopsis_microarray.ipynb`.

- **`arab_table.json`** — as-conducted metadata (tissue, ecotype, age, light,
  hardware, growth medium, factors, release year) for every *Arabidopsis
  thaliana* DNA-microarray dataset in OSDR, captured from the OSDR biodata
  metadata API (`/query/metadata/` and `/dataset/<acc>/`).
- **`de_embed.json`** — pre-computed differential-expression meta-analysis for
  the six spaceflight datasets with a clean wild-type Space-Flight vs
  Ground-Control contrast (OSD-44, 121, 147, 205, 213, 469; OSD-213 and OSD-469
  no longer resolve in the OSDR APIs as of 2026-10-22). Per-study DE
  counts, cross-study locus recurrence, a log2FC heatmap matrix for the top
  recurrent loci, and the consistent-direction core.
- **`collect_de.py`** — the collector that produced `de_embed.json`. It
  downloads each GeneLab DE table (20–800 MB; too large to fetch at book-build
  time), selects the wild-type flight-vs-ground contrast, collapses multi-mapping
  probes to the strongest per locus, and calls loci at FDR < 0.05 & |log2FC| ≥ 1.

Regenerate the notebook after editing data or the generator:

```
python chapters/_make_arabidopsis_microarray_notebook.py
jupyter nbconvert --to notebook --execute --inplace chapters/OSDR_arabidopsis_microarray.ipynb
```
