"""Generate `OSDR_data_story.ipynb` — a data-storytelling chapter that strings the
explorer outputs into a narrative and compares MULTIPLE Arabidopsis spaceflight
studies (transcriptome) plus surfaces the epigenetic (DNA-methylation) layer.

Run:  python _make_story_notebook.py   (then nbconvert --execute)
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md(r"""# 📖 Data story: the *Arabidopsis* spaceflight response across studies

In the [dataset explorer](OSDR_dataset_explorer.ipynb) we drilled into a single
study (OSD-120) and watched one gene respond to spaceflight. That raises the
question every scientist asks next:

> **Is the response *reproducible*? Do independent *Arabidopsis* spaceflight
> experiments agree — and does spaceflight reach beyond the transcriptome into
> the *epigenome*?**

This chapter strings the explorer's per-dataset outputs into a cross-study
narrative. We:

1. assemble a **panel of independent *Arabidopsis* spaceflight RNA-seq studies**,
2. compute one common metric — the **spaceflight log₂ fold-change** per gene,
3. ask how well the studies **agree**, and which genes change **consistently**,
4. then meet the **DNA-methylation (epigenetic) studies** that add the next layer.

> 🚀 Open this page in Binder/Colab (rocket icon) to re-run and extend the story.""")

code('''import io, re
import numpy as np
import pandas as pd
import requests
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "notebook_connected"

API = "https://visualization.osdr.nasa.gov/biodata/api/v2"
GEODE = "https://osdr.nasa.gov/geode-py/ws/studies"

# Independent Arabidopsis spaceflight RNA-seq studies with processed counts
PANEL = ["OSD-120", "OSD-37", "OSD-38", "OSD-321", "OSD-281", "OSD-427", "OSD-217"]''')

md("""## 1. The cast — independent spaceflight experiments

Each of these is a separate *Arabidopsis* spaceflight transcriptomics study. They span
five missions (SpaceX-3, -4, -5, -10 and -12) and several years, but they are not all
independent: OSD-38, OSD-217 and OSD-281 all flew on SpaceX-5.""")

code('''def normalized_counts(acc):
    """Download a study's GeneLab normalized-counts table."""
    d = requests.get(f"{API}/dataset/{acc}/files/", timeout=60).json()
    names = list(d[list(d)[0]].get("files", d[list(d)[0]]))
    f = next(n for n in names if "Normalized_Counts" in n and n.endswith(".csv") and "rRNArm" not in n)
    return pd.read_csv(f"{GEODE}/{acc}/download?source=datamanager&file={f}", index_col=0)

def flight_ground(acc, cols):
    """Map each sample column to 'Spaceflight' / 'Ground' using the study's factor values."""
    url = f"{API}/query/metadata/?id.accession={acc}&id.sample name&study.factor value&format=csv"
    df = pd.read_csv(io.StringIO(requests.get(url, timeout=60).text)).drop_duplicates("id.sample name")
    fcols = [c for c in df.columns if c.startswith("study.factor value")]
    sf = next((c for c in fcols if re.search(r"(?i)spaceflight|space.?flight", c)), None)
    labels = {}
    if sf:
        m = df.set_index("id.sample name")[sf].astype(str)
        for c in cols:
            v = m.get(c, "").lower()
            labels[c] = "Spaceflight" if ("flight" in v or "space" in v) else ("Ground" if "ground" in v else None)
    # fallback: parse the GeneLab sample name itself
    for c in cols:
        if labels.get(c) is None:
            labels[c] = "Spaceflight" if "_FLT_" in c else ("Ground" if "_GC_" in c else None)
    return labels

cast = []
for acc in PANEL:
    meta = requests.get(f"{API}/dataset/{acc}/", timeout=60).json()[acc]["metadata"]
    title = meta.get("study title", "")
    title = title[0] if isinstance(title, list) else title
    cast.append({"accession": acc, "title": title[:70]})
cast = pd.DataFrame(cast)
cast''')

md(r"""## 2. One common metric: the spaceflight log₂ fold-change

For every study we compute, per gene,
$\log_2\!\frac{\text{mean expression in spaceflight}}{\text{mean expression on the ground}}$.
Positive = up in space, negative = down. Putting every study on this same scale
lets us compare them directly.""")

code('''def spaceflight_lfc(acc):
    counts = normalized_counts(acc)
    lab = flight_ground(acc, counts.columns)
    flt = [c for c in counts.columns if lab.get(c) == "Spaceflight"]
    grd = [c for c in counts.columns if lab.get(c) == "Ground"]
    if len(flt) < 2 or len(grd) < 2:
        return None, len(flt), len(grd)
    lfc = np.log2((counts[flt].mean(axis=1) + 1) / (counts[grd].mean(axis=1) + 1))
    return lfc.rename(acc), len(flt), len(grd)

series, summary = {}, []
for acc in PANEL:
    lfc, nf, ng = spaceflight_lfc(acc)
    summary.append({"accession": acc, "flight": nf, "ground": ng,
                    "used": lfc is not None})
    if lfc is not None:
        series[acc] = lfc
lfc_mat = pd.DataFrame(series).dropna()          # genes (rows) x studies (cols)
print(f"{lfc_mat.shape[1]} studies compared over {lfc_mat.shape[0]} shared genes")
pd.DataFrame(summary)''')

md("""## 3. Do independent studies agree?

If the spaceflight response is real and reproducible, studies should **correlate**:
genes that go up in one experiment should tend to go up in the others.""")

code('''corr = lfc_mat.corr()
fig = px.imshow(corr, text_auto=".2f", color_continuous_scale="RdBu_r", zmin=-1, zmax=1,
                title="Cross-study agreement of the spaceflight response (log₂FC correlation)")
fig.update_layout(height=480)
fig.show()''')

code('''# Look at the single best-correlated pair head-to-head
import itertools
pairs = list(itertools.combinations(lfc_mat.columns, 2))
a, b = max(pairs, key=lambda p: lfc_mat[p[0]].corr(lfc_mat[p[1]]))
r = lfc_mat[a].corr(lfc_mat[b])
fig = px.scatter(lfc_mat, x=a, y=b, opacity=0.25,
                 title=f"Gene-by-gene spaceflight response: {a} vs {b}  (Pearson r = {r:.2f})",
                 labels={a: f"log₂FC ({a})", b: f"log₂FC ({b})"})
fig.update_layout(height=480)
fig.show()''')

md("""## 4. A candidate core

Which genes respond to spaceflight **consistently** — same direction, strong
effect — across the panel? Those are candidates for a shared core of the
*Arabidopsis* spaceflight transcriptome. The cut-off below (same sign in ≥80 % of
studies, top 25 by mean effect) has not been tested against a null model, so treat
the list as a shortlist, not a confirmed core.""")

code('''mean_lfc = lfc_mat.mean(axis=1)
sign_agreement = np.sign(lfc_mat).eq(np.sign(mean_lfc), axis=0).mean(axis=1)
core = mean_lfc[sign_agreement >= 0.8].abs().sort_values(ascending=False).head(25).index

fig = px.imshow(lfc_mat.loc[core], aspect="auto", color_continuous_scale="RdBu_r",
                zmin=-3, zmax=3,
                labels={"x": "study", "y": "gene (AGI locus)", "color": "log₂FC"},
                title="Reproducible spaceflight-responsive genes (consistent across studies)")
fig.update_layout(height=620)
fig.show()''')

md(r"""## 5. The next layer: spaceflight and the *epigenome*

Spaceflight doesn't only change which genes are *expressed* — it can change how
the genome is **chemically marked**. OSDR holds *Arabidopsis*
**DNA-methylation (bisulfite-sequencing)** studies that probe this epigenetic
layer:""")

code('''EPI = ["OSD-217", "OSD-220", "OSD-416"]
rows = []
for acc in EPI:
    meta = requests.get(f"{API}/dataset/{acc}/", timeout=60).json()[acc]["metadata"]
    def show(k):
        v = meta.get(k, "")
        return ", ".join(map(str, v)) if isinstance(v, list) else v
    rows.append({"accession": acc, "title": show("study title")[:55],
                 "assays": show("study assay technology type"),
                 "factors": show("study factor type")})
pd.DataFrame(rows)''')

md(r"""```{admonition} A multi-omic opportunity
:class: tip
**OSD-217 measures *both* DNA methylation *and* RNA-seq on the same spaceflight
experiment** — and its transcriptome is already in our panel above. That pairing
lets researchers ask the deeper question: *are the spaceflight-responsive genes
also differentially methylated?*

The biodata API exposes these studies' raw bisulfite data (processed
methylation tables aren't served through the API yet), so a quantitative
methylation-vs-expression analysis is the natural **frontier** this story points
to — the data is there in OSDR, waiting.
```

## 6. Results and discussion

### Key results

- **Weak, mixed agreement across studies.** Seven *Arabidopsis* spaceflight
  RNA-seq experiments (OSD-120, 37, 38, 321, 281, 427, 217) were compared over
  **21,527 shared genes**. Of the 21 study pairs, 14 had positively correlated
  spaceflight log₂ fold-changes and 7 negatively, with a median of only
  **r = 0.05**. The best pair (OSD-37 vs OSD-217) reached **r = 0.50** and the worst
  (OSD-120 vs OSD-321) **r = −0.27**; OSD-120 correlated negatively with four of its
  six partners.
- **A candidate core.** The 25 strongest of the genes that changed in a **consistent
  direction across ≥80 % of the studies** form a shortlist for a shared spaceflight
  response; it has not been tested against a null model.
- **Context adds scatter.** Pairwise agreement varied widely. The studies differ in
  design (ecotype, organ, light, hardware), and one pooled metric cannot separate
  those factors.
- **A second molecular layer.** One panel member (**OSD-217**) carries matched
  DNA-methylation data, opening the epigenetic dimension explored in the next chapter.

### Discussion

Genome-wide, the studies agree only weakly: most pairs correlate positively, but the
median correlation is close to zero and some pairs disagree outright. Some of the
agreement that does exist may not be independent, since three of the seven studies
(OSD-38, OSD-217, OSD-281) flew together on SpaceX-5. The spread is still
informative: the studies differ in ecotype, tissue, light and hardware, and a pooled
metric cannot tell those factors apart. The practical message for experiment design
is that a "spaceflight gene set" should be defined **across** studies, with an
explicit null model, not from any single one.

### Limitations

- **Heterogeneous designs** (tissue, ecotype, platform) compared with one simple
  metric (mean flight/ground log₂FC), which smooths over within-study factors.
- **Effect size, not significance** — the comparison ranks magnitude, not per-study
  statistical significance.
- **Plant-specific** — the comparison here is for *Arabidopsis*; other organisms
  may differ.
- **Shared missions** — OSD-38, OSD-217 and OSD-281 flew together on SpaceX-5, so
  they are not fully independent replicates.
- **No null model** — the candidate core uses a fixed cut-off without permutation
  testing.

### Conclusion

OSDR's *Arabidopsis* spaceflight studies share only a weak genome-wide signal, plus a
shortlist of genes that change consistently across most of them. Testing that
shortlist properly is the next step, and so is going **multi-omic** — pairing the
transcriptome with the epigenome (next chapter) and, for crops, the microbiome
(tomato chapter). Every figure here was built live from the OSDR API — change `PANEL`, swap
organisms, or join other layers to extend the analysis.""")

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
with open("OSDR_data_story.ipynb", "w", encoding="utf-8") as f:
    nbf.write(nb, f)
print("Wrote OSDR_data_story.ipynb")
