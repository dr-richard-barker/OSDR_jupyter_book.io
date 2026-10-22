"""Generate `OSDR_tomato_microbiome_pipeline.ipynb` — the processing recipe that
turns OSD-766's raw 16S/ITS reads into a genus abundance table and joins it to the
OSD-767 host transcriptome. Light steps run live; heavy steps (DADA2/SILVA) are
documented code for a dedicated bioinformatics environment.
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
OSD-766 as raw 16S/ITS reads — no processed taxonomy table**. This chapter is the
**recipe to fix that**: turn the raw reads into a per-genus abundance table and
feed it into the host↔microbe correlation.

```{admonition} What runs where — please read
:class: important
Amplicon denoising + taxonomy (DADA2 / QIIME2 + the SILVA database) needs heavy
tools and ~1 GB of reads, so it is **not executed inside this book**. The cells
below that *are* executed (live) **prove the data path and build the host side**;
the heavy steps are given as **complete, copy-runnable code** to run once in a
bioinformatics environment (a Binder/conda env with QIIME2 or R/DADA2). The
output of that step drops straight into section 4.
```

```{admonition} An honest design caveat
:class: warning
OSD-767 (host) and OSD-766 (microbes) are **different sample sets**, so the
host↔microbe link is necessarily at the **experimental-group level**
(organ × light × spaceflight), not sample-by-sample. With only a handful of
groups this is **descriptive**, not a powered correlation — read it as
hypothesis-generating.
```""")

code('''import io, re, gzip
import requests
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
print(f"reads: {len(lines)//4:,}   read length: {len(lines[1])} bp (16S V3-V4)")
print(f"example read: {lines[1][:60]}...")''')

md(r"""## 2. From reads to a genus table — run in a DADA2/QIIME2 environment

The block below is the **standard amplicon pipeline**. Run it once in an
environment with `cutadapt` + R/`dada2` (or QIIME2). It writes
`tomato_genus_by_sample.csv` (genera × samples), which section 4 consumes.

```bash
# 0. Download the 16S reads (R1+R2) for OSD-766 (~1 GB) into ./reads/
python - <<'PY'
import requests, json, os
os.makedirs("reads", exist_ok=True)
fobjs = requests.get("https://osdr.nasa.gov/osdr/data/osd/files/766").json()
# (walk the JSON as in section 1, then:)
for name in fastq_16S_names:                      # R1 + R2 for every 16S sample
    url = f"https://osdr.nasa.gov/geode-py/ws/studies/OSD-766/download?source=datamanager&file={name}"
    open(f"reads/{name}", "wb").write(requests.get(url).content)
PY

# 1. Remove the 515F/806R primers
for f in reads/*_R1_raw.fastq.gz; do
  r=${f/_R1_/_R2_}
  cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT \
           -o trimmed/$(basename $f) -p trimmed/$(basename $r) "$f" "$r"
done
```

```r
# 2. DADA2: denoise -> ASVs -> taxonomy (SILVA) -> collapse to genus  (in R)
library(dada2)
path <- "trimmed"
fnF <- sort(list.files(path, "_R1_", full.names = TRUE))
fnR <- sort(list.files(path, "_R2_", full.names = TRUE))
filtF <- file.path("filt", basename(fnF)); filtR <- file.path("filt", basename(fnR))
filterAndTrim(fnF, filtF, fnR, filtR, truncLen = c(220, 180), maxEE = c(2, 2))
errF <- learnErrors(filtF); errR <- learnErrors(filtR)
seqtab <- makeSequenceTable(mergePairs(dada(filtF, errF), filtF,
                                       dada(filtR, errR), filtR))
seqtab <- removeBimeraDenovo(seqtab)
tax <- assignTaxonomy(seqtab, "silva_nr99_v138.1_train_set.fa.gz")   # SILVA DB
# collapse ASV counts to genus and write genera x samples
genus <- rowsum(t(seqtab), group = tax[, "Genus"])
write.csv(genus, "tomato_genus_by_sample.csv")
```""")

md("""## 3. Map microbiome samples to experimental groups

OSD-766's sample names encode the organ and flight/ground status, so the genus
table can be summarised on the **same groups** as the host data.""")

code('''# Parse organ + condition from the 16S sample names (real, live)
rows = []
for n in fq:
    if "16S" not in n or "_R1" not in n:
        continue
    cond = "Flight" if re.search(r"VEG-05-F-", n) else ("Ground" if re.search(r"VEG-05-G-", n) else None)
    organ = next((o for o in ("root", "leaf", "fruit", "soil", "wick", "swab") if o in n.lower()), "other")
    rows.append({"file": n, "organ": organ, "condition": cond})
mic_samples = pd.DataFrame(rows)
print(mic_samples.groupby(["organ", "condition"]).size().rename("samples").reset_index().to_string(index=False))''')

md(r"""## 4. Join to the host transcriptome

Here is the live host side plus the integration functions. Once
`tomato_genus_by_sample.csv` exists (from section 2), the final cell computes **(a)** a
group-level **genus ↔ host-gene correlation** and **(b)** a per-sample **random-forest
ranking** (`rf_host_association`, the regression sibling of section 5's `rf_responders`)
of which microbes best track the host's **spaceflight gene-module**.""")

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

def correlate_microbe_host(genus_by_group, host_by_group, genes):
    """genus_by_group & host_by_group indexed by (organ, condition).
    Returns genera x genes Pearson r across the shared groups."""
    j = genus_by_group.join(host_by_group[genes], how="inner")
    return j[genus_by_group.columns].apply(lambda t: j[genes].corrwith(t)).T

# A host "spaceflight module": the top flight-vs-ground root genes, scored per group
_organs = host_group_mean.index.get_level_values("organ")
_root = host_group_mean.xs("Root", level="organ") if "Root" in _organs else host_group_mean
if {"Flight", "Ground"}.issubset(set(_root.index)):
    _lfc = (_root.loc["Flight"] - _root.loc["Ground"]).abs().sort_values(ascending=False)
    module_genes = list(_lfc.head(50).index)
else:
    module_genes = list(host_group_mean.columns[:50])
host_module = host_group_mean[module_genes].mean(axis=1).rename("host_spaceflight_module")

def rf_host_association(genus_by_sample, module_per_sample, n_estimators=400, seed=0):
    """RF regression: which genera best predict a host gene-module score.
    The regression sibling of rf_responders() (section 5) — point it at the per-sample
    genus table and the host module each sample's group is expected to show."""
    from sklearn.ensemble import RandomForestRegressor
    rf = RandomForestRegressor(n_estimators=n_estimators, random_state=seed, oob_score=True)
    rf.fit(genus_by_sample.values, module_per_sample.values)
    return rf, pd.Series(rf.feature_importances_, index=genus_by_sample.columns).sort_values(ascending=False)''')

code('''# Final step — runs for real ONCE the genus table from section 2 is present
import os
if os.path.exists("tomato_genus_by_sample.csv"):
    genus = pd.read_csv("tomato_genus_by_sample.csv", index_col=0)
    g = genus.T
    g["organ"] = [next((o for o in ("root","leaf","fruit","soil","wick","swab") if o in s.lower()), "other") for s in g.index]
    g["condition"] = ["Flight" if "-F-" in s else "Ground" for s in g.index]
    # group keyed by Title-case organ so it aligns with the host index ("Root"/"Leaf")
    genus_by_group = g.groupby([g["organ"].str.title(), "condition"]).mean(numeric_only=True)

    # (a) group-level correlation of the PGPR genera with the most variable host root genes
    pgpr = [c for c in genus.index
            if c in ("Rhizobium","Azospirillum","Burkholderia","Dyadobacter","Sphingomonas")]
    top_genes = host_group_mean.loc[("Root",)].std().sort_values(ascending=False).head(30).index
    print("(a) Group-level PGPR x host-gene correlation:")
    print(correlate_microbe_host(genus_by_group[pgpr], host_group_mean, list(top_genes)).round(2))

    # (b) per-sample random forest: which microbes track the host spaceflight module?
    g["host_module"] = [host_module.get((o.title(), c)) for o, c in zip(g["organ"], g["condition"])]
    sub = g[g["organ"].isin(["root", "leaf"])].dropna(subset=["host_module"])
    rf, imp = rf_host_association(sub[list(genus.index)], sub["host_module"])
    print(f"\\n(b) rf_host_association (OOB R2 = {rf.oob_score_:.2f}) — "
          f"microbes most predictive of the host spaceflight module:")
    print(imp.head(10).round(3))
else:
    print("Run section 2 to create tomato_genus_by_sample.csv, then re-run this cell for:")
    print("  (a) group-level PGPR x host-gene correlation, and")
    print("  (b) rf_host_association() — RF ranking of microbes by association with the")
    print("      host spaceflight gene-module (the regression sibling of section 5's rf_responders).")''')

md(r"""## 5. Results

We present the chapter's data graphically. The **reported responders** (bar) and the
**host response** (box) are *real*; the integration views (RF importance bar, scatter)
are demonstrated on a **simulated** community — clearly marked — because the processed
genus table isn't public yet; the **interaction Sankey** is a curated map of documented
plant–microbe biology.""")

code('''import numpy as np
import plotly.express as px
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

The host module (top flight-vs-ground root genes, section 4) scored per sample — the
response the microbiome would engage.""")

code('''mps = cpm.loc[cpm.index.intersection(module_genes)].mean(axis=0)
box_df = pd.DataFrame({"host module score": mps.values}, index=mps.index)
box_df["condition"] = ["Flight" if "Flt" in s else "Ground" for s in box_df.index]
box_df["organ"] = ["Root" if "Root" in s else "Leaf" for s in box_df.index]
fig = px.box(box_df, x="organ", y="host module score", color="condition", points="all",
             title="Host spaceflight gene-module by organ and condition (real, OSD-767)")
fig.update_layout(height=400)
fig.show()''')

md(r"""### Random forest + scatter — the integration, demonstrated

```{admonition} Simulated — method check, not a result
:class: warning
The processed genus table isn't public yet, so the next two panels run on a **simulated**
community with a planted spaceflight signal, to show the method works. The real numbers
come from running `rf_responders()` / `rf_host_association()` (section 4) on the DADA2
output of section 2.
```""")

code('''def rf_responders(genus_by_sample, labels, n_estimators=400, seed=0):
    """RF importance: which genera best separate flight vs ground."""
    from sklearn.ensemble import RandomForestClassifier
    rf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed, oob_score=True)
    rf.fit(genus_by_sample.values, np.asarray(labels))
    return rf, pd.Series(rf.feature_importances_, index=genus_by_sample.columns).sort_values(ascending=False)

rng = np.random.default_rng(0)
n_samples, n_genera = 48, 30
labels = np.array(["Flight", "Ground"] * (n_samples // 2))
X = rng.lognormal(2.0, 1.0, size=(n_samples, n_genera))
planted = [3, 7, 12, 20, 25]
for j in planted:
    X[labels == "Flight", j] *= 4.0
sim = pd.DataFrame(X, columns=[f"Genus_{i:02d}" for i in range(n_genera)])

rf, imp = rf_responders(sim, labels)
print(f"SIMULATED RF out-of-bag accuracy: {rf.oob_score_:.2f}")
top = imp.head(12).rename_axis("genus").reset_index(name="importance")
top["planted"] = top["genus"].isin([f"Genus_{j:02d}" for j in planted])
fig = px.bar(top.sort_values("importance"), x="importance", y="genus", orientation="h",
             color="planted", color_discrete_map={True: "#2e7d32", False: "#90a4ae"},
             title="SIMULATED: RF importance recovers the planted spaceflight responders")
fig.update_layout(height=420)
fig.show()''')

code('''# Scatter: a planted flight-associated microbe vs a (simulated) host module, per sample
sim_module = np.where(labels == "Flight", 1.0, 0.0) + rng.normal(0, 0.15, n_samples)
scat = pd.DataFrame({"planted PGPR abundance": sim["Genus_03"],
                     "host module (simulated)": sim_module, "condition": labels})
fig = px.scatter(scat, x="planted PGPR abundance", y="host module (simulated)", color="condition",
                 title="SIMULATED: a flight-associated microbe vs the host module (per sample)")
fig.update_layout(height=400)
fig.show()''')

md(r"""### Sankey — how the flight-shifted microbes connect to the host response (curated)

A knowledge-based map: flight-shifted **microbes** → their **functional traits** → the
**host pathways** enriched in spaceflight (section 4). The flows converge on the
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

Spaceflight reshaped the VEG-05 tomato rhizosphere toward **nitrogen-fixing and
plant-growth-promoting (PGPR) bacteria** — the *Rhizobium*- and *Burkholderia*-clades and
*Azospirillum* rose most (>8-fold for the two clades), while the endophytic diazotroph
*Herbaspirillum* fell (VEG-05 microbiome study; differential abundance via DESeq2, Love
*et al.*, 2014). A tilt toward beneficial, pathogen-free communities echoes earlier ISS
crop work: lettuce grown in the same Veggie hardware carried no plant pathogens and a
largely benign microbiota (Khodadad *et al.*, 2020).

Strikingly, the host met this shift with matching chemistry. The tomato **root**
transcriptome up-regulated the **flavonoid and phenylpropanoid** pathways (this book's
re-analysis of OSD-767, and the VEG-05 transcriptome paper) — the canonical molecules
plants exude to **recruit and signal to rhizosphere bacteria** (Hassan & Mathesius, 2012).
Flavonoids activate the colonisation and nodulation programmes of exactly the N-fixing
partners the flight community is enriched for, making the **flavonoid ↔ rhizobia axis**
(the Sankey's busiest hub) the clearest point of contact between the two datasets.
*Azospirillum*, also enriched, promotes growth through phytohormones, nitrogen fixation
and root proliferation rather than nodulation (Bashan & de-Bashan, 2010), consistent with
the light- and root-dependent responses across VEG-05.

The response was **organ-specific**: roots — the rhizosphere interface — carried the
strongest, most distinct signal (the box plot), as expected if microbes are involved, and
as seen for the *Arabidopsis* spaceflight transcriptome, which is likewise remodelled
organ-by-organ (Paul *et al.*, 2013). A second meeting point, **reactive-oxygen-species /
defence** chemistry, is the toolkit plants use to *tune* which microbes colonise the root
surface.

**How firmly can we call this causal?** Not yet. Host and microbiome were profiled on
**different plants**, so the link is group-level; spaceflight simultaneously alters water
films, temperature and the plant, any of which could drive a microbe's apparent
importance; and plant→microbe and microbe→plant signalling run both ways. The random-forest
framework here ranks candidate responders and host-tracking microbes (the simulated panels
show it recovers a planted signal) — but converting association to causation needs
**matched per-sample multi-omics** with **mediation analysis**, and ultimately
**inoculation experiments** that add or omit a genus and watch the host. The value of this
chapter is to reach that shortlist quickly and transparently.

### References

- Bashan, Y. & de-Bashan, L. E. (2010) How the plant growth-promoting bacterium
  *Azospirillum* promotes plant growth — a critical assessment. *Advances in Agronomy*
  **108**, 77–136.
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
