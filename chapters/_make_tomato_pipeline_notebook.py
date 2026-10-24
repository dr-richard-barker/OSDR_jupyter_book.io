"""Generate `OSDR_tomato_microbiome_pipeline.ipynb` — turns OSD-766's raw 16S reads into a
genus abundance table and sets it beside the OSD-767 host transcriptome. The heavy step
(DADA2 + SILVA) is run with the scripts in data_tomato_16s/, whose outputs are committed;
the book reads those outputs and runs the light steps live.
Run: python _make_tomato_pipeline_notebook.py
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md(r"""# 🧫 Processing the tomato microbiome → host integration

The [tomato chapter](OSDR_tomato_microbiome.ipynb) paired the host transcriptome
(OSD-767) with the microbiome study (OSD-766) but hit a wall: **OSDR serves
OSD-766 as raw 16S/ITS reads — no processed taxonomy table**. This chapter builds
that table from the raw 16S reads (section 2) and sets it beside the host data.

```{admonition} What runs where — please read
:class: important
Denoising and taxonomy (DADA2 + the SILVA database) need R, about 330 MB of reads and
a 137 MB reference, so they are **not executed inside this book**. They were run with
the four scripts in
[`data_tomato_16s/`](https://github.com/dr-richard-barker/OSDR_jupyter_book.io/tree/main/chapters/data_tomato_16s),
and their outputs are committed next to them. The cells below read those outputs; the
live cells check the data path and build the host side.
```

```{admonition} An honest design caveat
:class: warning
OSD-767 (host) and OSD-766 (microbes) are **different sample sets**, so host and
microbes can only be compared at the **experimental-group level**, not sample by
sample. Almost all leaf reads also turned out to come from the tomato's own chloroplasts
and mitochondria (section 2),
which leaves the root tissues as the only usable microbiome samples. Read the
host–microbe comparison as descriptive and hypothesis-generating.
```""")

code('''import io, re, gzip, os
import requests
import numpy as np
import pandas as pd

API = "https://visualization.osdr.nasa.gov/biodata/api/v2"
GEODE = "https://osdr.nasa.gov/geode-py/ws/studies"''')

md("## 1. The raw amplicon data is there (live check)\n\nWe enumerate OSD-766's raw read files straight from OSDR.")

code('''# Full file listing (the OSD file API lists more than the biodata files endpoint)
listing = requests.get("https://osdr.nasa.gov/osdr/data/osd/files/766", timeout=60).json()
def collect(o, acc):
    if isinstance(o, dict):
        if "file_name" in o: acc.append(o)
        for v in o.values(): collect(v, acc)
    elif isinstance(o, list):
        for x in o: collect(x, acc)
    return acc
fobjs = collect(listing, [])
fq = [o["file_name"] for o in fobjs if o["file_name"].endswith(".fastq.gz")]
print(f"OSD-766 raw FASTQ files: {len(fq)}")
print(f"  16S (bacteria): {sum('16S' in n for n in fq)}   |   ITS (fungi): {sum('ITS' in n for n in fq)}")''')

code('''# Download ONE 16S sample and confirm what the reads look like (real, ~few MB)
one = next(o["file_name"] for o in fobjs
           if o["file_name"].endswith("R1_raw.fastq.gz") and "16S" in o["file_name"]
           and "root" in o["file_name"] and float(o.get("file_size", 0)) > 1e6)
blob = requests.get(f"{GEODE}/OSD-766/download?source=datamanager&file={one}", timeout=120).content
lines = gzip.decompress(blob).decode("ascii", "replace").splitlines()
print(f"file : {one}")
print(f"reads: {len(lines)//4:,}   read length: {len(lines[1])} bp")
print(f"example read: {lines[1][:60]}...")''')

md(r"""## 2. From reads to a genus table (DADA2 + SILVA)

The table was built with the scripts in `data_tomato_16s/`, run in this order:

```bash
python download_reads.py ~/OSD-766-16S/reads     # 44 per-plant root, adventitious-root and leaf samples (88 FASTQ)
python inspect_reads.py  ~/OSD-766-16S/reads     # primers? read length? where does quality drop?
Rscript run_dada2.R ~/OSD-766-16S/reads silva_nr99_v138.1_train_set.fa.gz ~/OSD-766-16S/dada2 .
python make_genus_table.py                       # sample sheet + relative-abundance table
```

What the inspection and the run found, and the settings used:

- **No primers on the reads.** None of the common 16S primers occurs at the start of the
  deposited reads, so there is no trimming step. Most merged amplicons are 253 bp (the
  cell below), the length of the V4 region between the 515F and 806R primers.
- **DADA2**: `truncLen = c(220, 180)`, `maxEE = c(2, 2)`, consensus chimera removal.
- **Taxonomy**: SILVA 138.1 ([Zenodo 10.5281/zenodo.4587955](https://doi.org/10.5281/zenodo.4587955)),
  with DADA2's default rule that a rank needs 50% bootstrap support to be named.
- **Host organelles removed**: chloroplast and mitochondrial ASVs are dropped before the
  reads are summed by genus.
- **Depth filter**: a sample needs at least 1,000 bacterial reads to enter the analysis.

The genus table holds each kept sample's relative abundance (% of its bacterial reads),
so every column sums to 100 by construction.""")

code('''samples = pd.read_csv("data_tomato_16s/samples.csv").set_index("sample")
track = pd.read_csv("data_tomato_16s/read_tracking.csv").set_index("sample")
asv = pd.read_csv("data_tomato_16s/asv_taxonomy.csv").set_index("asv")

tot = track[["input", "filtered", "merged", "nonchim", "organelle", "bacterial"]].sum()
print(f"Reads kept through DADA2: filtered {tot.filtered / tot.input:.0%}, merged {tot.merged / tot.input:.0%}, "
      f"non-chimeric {tot.nonchim / tot.input:.0%} of {int(tot.input):,} input reads")
print(f"Chloroplast/mitochondrial share of non-chimeric reads: {tot.organelle / tot.nonchim:.0%}")
lens = asv["sequence"].str.len()
by_len = asv.groupby(lens)["total_reads"].sum().sort_values(ascending=False)
print(f"Merged amplicon length carrying most bacterial reads: {by_len.index[0]} bp "
      f"({by_len.iloc[0] / by_len.sum():.0%} of them)")

samples["organelle %"] = 100 * samples["organelle_reads"] / (samples["organelle_reads"] + samples["bacterial_reads"])
summary = samples.groupby(["condition", "tissue"]).agg(
    samples=("kept", "size"), kept=("kept", "sum"),
    median_bacterial_reads=("bacterial_reads", "median"), median_organelle_pct=("organelle %", "median"))
summary.round(1)''')

md("""## 3. Map microbiome samples to experimental groups

OSDR's metadata records each 16S sample's spaceflight and light factors;
`make_genus_table.py` saved them in `data_tomato_16s/samples.csv`, and the tissue comes
from the sample name. These are the samples that pass the depth filter:""")

code('''kept = samples[samples["kept"]]
print(kept.groupby(["tissue", "condition", "light"]).size().rename("samples").reset_index().to_string(index=False))''')

md(r"""## 4. Join to the host transcriptome

The host side runs live from OSD-767: per-group expression and a **spaceflight
gene-module** (the 50 root genes that differ most between flight and ground). The
microbiome side is the genus table from section 2. With leaves excluded, the microbes
can be compared only within the two root tissues, so the cells below ask which genera
differ between flight and ground roots.""")

code('''# Host (OSD-767) expression, summarised per group  -- runs live
counts = pd.read_csv(f"{GEODE}/OSD-767/download?source=datamanager"
                     "&file=GLDS-709_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv", index_col=0)
cpm = counts / counts.sum() * 1e6
cpm = cpm[(cpm > 1).sum(axis=1) >= 3]
def parse(c):
    return ("Flight" if "Flt" in c else "Ground"), ("Root" if "Root" in c else "Leaf")
groups = pd.DataFrame({c: parse(c) for c in cpm.columns}, index=["condition", "organ"]).T
host_group_mean = cpm.T.join(groups).groupby(["organ", "condition"]).mean(numeric_only=True)
print(f"Host group means: {host_group_mean.shape[0]} groups x {host_group_mean.shape[1]} genes")

# A host "spaceflight module": the top flight-vs-ground root genes
_organs = host_group_mean.index.get_level_values("organ")
_root = host_group_mean.xs("Root", level="organ") if "Root" in _organs else host_group_mean
if {"Flight", "Ground"}.issubset(set(_root.index)):
    _lfc = (_root.loc["Flight"] - _root.loc["Ground"]).abs().sort_values(ascending=False)
    module_genes = list(_lfc.head(50).index)
else:
    module_genes = list(host_group_mean.columns[:50])''')

code('''# The genus table from section 2: relative abundance (% of bacterial reads), genera x samples
from scipy.stats import mannwhitneyu
genus = pd.read_csv("tomato_genus_by_sample.csv", index_col=0)
meta = samples.loc[genus.columns]
roots = meta[meta["tissue"].isin(["root", "adventitious root"])]
print(f"{genus.shape[0]} genera x {genus.shape[1]} samples; root tissues: {len(roots)} samples "
      f"(flight {int((roots.condition == 'Flight').sum())}, ground {int((roots.condition == 'Ground').sum())})")

# The genera the VEG-05 microbiome preprint names, flight vs ground, in each root tissue (SILVA names)
NAMED = {"Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium": "Rhizobium clade",
         "Burkholderia-Caballeronia-Paraburkholderia": "Burkholderia clade",
         "Azospirillum": "Azospirillum", "Dyadobacter": "Dyadobacter", "Sphingomonas": "Sphingomonas",
         "Methylobacterium-Methylorubrum": "Methylobacterium", "Massilia": "Massilia",
         "Curtobacterium": "Curtobacterium", "Herbaspirillum": "Herbaspirillum"}
rows = []
for tissue in ("root", "adventitious root"):
    f = roots.index[(roots.tissue == tissue) & (roots.condition == "Flight")]
    g = roots.index[(roots.tissue == tissue) & (roots.condition == "Ground")]
    for silva, label in NAMED.items():
        a, b = genus.loc[silva, f], genus.loc[silva, g]
        rows.append({"tissue": tissue, "genus": label, "n (F vs G)": f"{len(f)} vs {len(g)}",
                     "flight %": a.mean(), "ground %": b.mean(),
                     "log2 F/G": np.log2((a.mean() + 0.01) / (b.mean() + 0.01)),
                     "p (Mann-Whitney)": mannwhitneyu(a, b).pvalue})
named = pd.DataFrame(rows)
named.round(3)''')

md(r"""The **Rhizobium clade** needs care. Flight and ground plants carry *different*
Rhizobiaceae sequence variants (ASVs), and SILVA can name the flight one only with low
confidence:""")

code('''# Rhizobiaceae ASVs: genus as assigned (50% bootstrap rule), best guess, support, and reads by condition
asv_counts = pd.read_csv("data_tomato_16s/asv_counts.csv", index_col=0)
fl = samples.index[samples.condition == "Flight"]; gr = samples.index[samples.condition == "Ground"]
rz = asv[asv["Family"] == "Rhizobiaceae"].copy()
rz["flight reads"] = asv_counts.loc[rz.index, fl].sum(axis=1)
rz["ground reads"] = asv_counts.loc[rz.index, gr].sum(axis=1)
rz[["Genus", "best_genus", "genus_bootstrap", "flight reads", "ground reads"]].sort_values("flight reads", ascending=False)''')

md(r"""## 5. Results

The **reported responders** (bar) come from the VEG-05 microbiome preprint and the
**host response** (box) from OSD-767. The **random forest** and the **strip plot** use the
genus table built in section 2. The **interaction Sankey** is a curated map of
documented plant–microbe biology.""")

code('''import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "notebook_connected"   # keep charts interactive in the static book''')

md("### Bar — which microbes respond most to spaceflight (real, from the paper)")

code('''# The VEG-05 DESeq2 result, encoded (only the two clades were quantified, at >8-fold)
responders = pd.DataFrame([
    ("Allorhizobium-Rhizobium clade", 3.0), ("Burkholderia clade", 3.0),
    ("Azospirillum", 1.0), ("Sphingomonas", 1.0), ("Dyadobacter", 1.0),
    ("Methylobacterium", 1.0), ("Massilia", 1.0), ("Curtobacterium", 1.0),
    ("Herbaspirillum", -1.0),
], columns=["genus", "reported flight response"])
fig = px.bar(responders.sort_values("reported flight response"),
             x="reported flight response", y="genus", orientation="h",
             color="reported flight response", color_continuous_scale="RdBu",
             color_continuous_midpoint=0,
             title="Which microbes respond most to spaceflight (VEG-05 DESeq2, reported)")
fig.update_layout(height=400, coloraxis_showscale=False)
fig.update_xaxes(title="ordinal: +3 = >8-fold up, +1 = up, -1 = down")
fig.show()''')

md(r"""### Box — the host spaceflight gene-module (real, OSD-767)

The host module (top flight-vs-ground root genes, section 4) scored per sample.""")

code('''mps = cpm.loc[cpm.index.intersection(module_genes)].mean(axis=0)
box_df = pd.DataFrame({"host module score": mps.values}, index=mps.index)
box_df["condition"] = ["Flight" if "Flt" in s else "Ground" for s in box_df.index]
box_df["organ"] = ["Root" if "Root" in s else "Leaf" for s in box_df.index]
fig = px.box(box_df, x="organ", y="host module score", color="condition", points="all",
             title="Host spaceflight gene-module by organ and condition (real, OSD-767)")
fig.update_layout(height=400)
fig.show()''')

md(r"""### Random forest — which genera separate flight from ground roots (real)

A random forest trained on the root-tissue samples (roots and adventitious roots
together) to tell flight from ground; the out-of-bag accuracy says how well it does on
samples it did not train on. Flight has more adventitious-root samples than ground, so
tissue and condition are partly confounded here.""")

code('''def rf_responders(genus_by_sample, labels, n_estimators=400, seed=0):
    """RF importance: which genera best separate flight vs ground."""
    from sklearn.ensemble import RandomForestClassifier
    rf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed, oob_score=True)
    rf.fit(genus_by_sample.values, np.asarray(labels))
    return rf, pd.Series(rf.feature_importances_, index=genus_by_sample.columns).sort_values(ascending=False)

X = genus[roots.index].T
rf, imp = rf_responders(X, roots["condition"])
print(f"Root tissues, flight vs ground (n = {len(roots)}): out-of-bag accuracy {rf.oob_score_:.2f}")
top = imp.head(12).rename_axis("genus").reset_index(name="importance")
fig = px.bar(top.sort_values("importance"), x="importance", y="genus", orientation="h",
             title="Genera that best separate flight from ground root tissues (random forest)")
fig.update_layout(height=440)
fig.show()''')

md("### Strip plot — the named genera in each root sample (real)")

code('''long = (genus.loc[list(NAMED), roots.index].rename(index=NAMED).T
        .join(roots[["condition", "tissue"]]).melt(id_vars=["condition", "tissue"],
                                                   var_name="genus", value_name="relative abundance (%)"))
fig = px.strip(long, x="relative abundance (%)", y="genus", color="condition", facet_col="tissue",
               title="Paper-named genera in each root-tissue sample (DADA2 + SILVA, this chapter)")
fig.update_layout(height=460)
fig.show()''')

md(r"""### Sankey — how the flight-shifted microbes connect to the host response (curated)

A knowledge-based map: flight-shifted **microbes** → their **functional traits** → the
**host pathways** enriched in spaceflight (the tomato chapter). The flows converge on the
flavonoid/phenylpropanoid and ROS/defence axes.""")

code('''micro = ["Rhizobium clade", "Azospirillum", "Burkholderia clade", "Sphingomonas",
         "Methylobacterium", "Dyadobacter", "Massilia"]
traits = ["N-fixation", "Flavonoid response", "ROS interaction", "Stress tolerance"]
hostp = ["Flavonoid / phenylpropanoid", "ROS / defence", "Stress"]
nodes = micro + traits + hostp
idx = {n: i for i, n in enumerate(nodes)}
links = [
    ("Rhizobium clade", "N-fixation"), ("Rhizobium clade", "Flavonoid response"),
    ("Azospirillum", "N-fixation"), ("Azospirillum", "Stress tolerance"),
    ("Burkholderia clade", "N-fixation"), ("Burkholderia clade", "ROS interaction"),
    ("Sphingomonas", "ROS interaction"), ("Sphingomonas", "Stress tolerance"),
    ("Methylobacterium", "Stress tolerance"), ("Dyadobacter", "ROS interaction"),
    ("Massilia", "ROS interaction"),
    ("N-fixation", "Flavonoid / phenylpropanoid"), ("Flavonoid response", "Flavonoid / phenylpropanoid"),
    ("ROS interaction", "ROS / defence"), ("Stress tolerance", "Stress"),
]
fig = go.Figure(go.Sankey(
    node=dict(label=nodes, pad=15, thickness=16,
              color=["#66bb6a"] * len(micro) + ["#90caf9"] * len(traits) + ["#ffb74d"] * len(hostp)),
    link=dict(source=[idx[a] for a, b in links], target=[idx[b] for a, b in links],
              value=[1] * len(links))))
fig.update_layout(height=470,
                  title="Flight-shifted microbes -> traits -> host spaceflight-enriched pathways")
fig.show()''')

md(r"""## 6. Discussion

Building the genus table from OSD-766's raw reads lets this chapter test the microbiome
side of the story directly. Within the root tissues, four of the VEG-05 preprint's
findings hold up (section 4):

- The ***Burkholderia* clade** is higher in flight in roots (18.7% vs 2.9% of bacterial
  reads) and adventitious roots (38.0% vs 8.8%).
- ***Azospirillum*** and ***Dyadobacter*** are higher in flight in both root tissues.
- ***Herbaspirillum***, the one genus the preprint reports as lower in flight, is absent
  from every flight root-tissue sample.

*Methylobacterium* and *Curtobacterium* lean higher in flight without reaching p < 0.05,
and *Sphingomonas* and *Massilia* show no consistent direction, at these sample sizes
(4 vs 4 roots, 9 vs 3 adventitious roots). The preprint's differential abundance used
DESeq2 (Love *et al.*, 2014); the comparisons here are simpler, a Mann-Whitney test on
relative abundance.

**The *Rhizobium* clade is the exception, and the reason is taxonomic.** Flight and ground
plants carry different Rhizobiaceae sequence variants. The ground variant (ASV7) is named
*Rhizobium* clade with 95% bootstrap support. The main flight variant (ASV4, 15,898 reads,
none on the ground) is closest to the same clade but with only 46% support, below DADA2's
50% cut-off, so it stays unnamed at genus level. At genus level the clade therefore looks
absent from flight and abundant on the ground, the reverse of the preprint, while
Rhizobiaceae as a family are plentiful in both. The preprint classified with Greengenes2
and pooled plant and non-plant samples in its differential-abundance plots, so the two
analyses need not agree on this clade. Do not read it as a reversal of the preprint.

A random forest separates flight from ground root tissues well (out-of-bag accuracy 0.95,
20 samples). Its leading genera are *Xylophilus*, the *Rhizobium* clade (subject to the
caveat above), *Chitinophaga*, *Paenarthrobacter* and *Dyadobacter*. But flight
contributes 9 adventitious-root samples against 3 for ground, so tissue and condition
are partly confounded: the forest shows that the communities differ, not which genus
drives the difference.

On the host side, the tomato root genes that change most in spaceflight (in either
direction) enrich for flavonoid and phenylpropanoid biosynthesis (the tomato chapter).
These are the molecules plants use to recruit and signal to rhizosphere bacteria (Hassan &
Mathesius, 2012). With the *Burkholderia* clade and *Azospirillum*, a well-studied plant
growth-promoting bacterium (Bashan & de-Bashan, 2010), higher on flight roots, that
flavonoid–rhizobacteria axis is the clearest candidate point of contact between the two
datasets. As in *Arabidopsis*, whose spaceflight transcriptome is remodelled organ by
organ (Paul *et al.*, 2013), the signal here is organ-specific, and the microbes can only
be read in roots. Earlier Veggie crops were screened too: lettuce grown in the same
hardware tested negative for select human pathogens (Khodadad *et al.*, 2020).

**How firmly can we call this causal?** Not yet. Host and microbiome were profiled on
**different plants**, so the link is group-level. Most leaf samples carry too few bacterial
reads to analyse, and the ground adventitious roots number only three. Spaceflight also
changes water, temperature and the plant at once, any of which could move a microbe.
Turning these associations into causes needs **matched per-sample host and microbiome
data**, and ultimately **inoculation experiments** that add or omit a genus and watch the
host.

### References

- Bashan, Y. & de-Bashan, L. E. (2010) How the plant growth-promoting bacterium
  *Azospirillum* promotes plant growth — a critical assessment. *Advances in Agronomy*
  **108**, 77–136. <https://doi.org/10.1016/S0065-2113(10)08002-8>
- Hassan, S. & Mathesius, U. (2012) The role of flavonoids in root–rhizosphere signalling:
  opportunities and challenges for improving plant–microbe interactions. *Journal of
  Experimental Botany* **63**(9), 3429–3444.
  <https://doi.org/10.1093/jxb/err430>
- Khodadad, C. L. M. *et al.* (2020) Microbiological and nutritional analysis of lettuce
  crops grown on the International Space Station. *Frontiers in Plant Science* **11**, 199.
  <https://doi.org/10.3389/fpls.2020.00199>
- Love, M. I., Huber, W. & Anders, S. (2014) Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biology* **15**, 550.
  <https://doi.org/10.1186/s13059-014-0550-8>
- Paul, A.-L., Zupanska, A. K., Schultz, E. R. & Ferl, R. J. (2013) Organ-specific
  remodeling of the *Arabidopsis* transcriptome in response to spaceflight. *BMC Plant
  Biology* **13**, 112. <https://doi.org/10.1186/1471-2229-13-112>
- VEG-05 microbiome: *The microbiome of a tomato crop grown under different lighting
  regimes on the ISS.* NASA NTRS 20240016407.
- VEG-05 transcriptome: *Stress and light spectral quality influence the transcriptome of a
  tomato crop on the ISS.* *BMC Plant Biology* (2026).
  <https://doi.org/10.1186/s12870-025-07621-4>

---

**Data:** [OSD-766 (microbiome)](https://osdr.nasa.gov/bio/repo/data/studies/OSD-766) ·
[OSD-767 (host RNA-seq)](https://osdr.nasa.gov/bio/repo/data/studies/OSD-767)""")

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
with open("OSDR_tomato_microbiome_pipeline.ipynb", "w", encoding="utf-8") as f:
    nbf.write(nb, f)
print("Wrote OSDR_tomato_microbiome_pipeline.ipynb")
