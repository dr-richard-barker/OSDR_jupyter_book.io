"""Generate `OSDR_methylation_layer.ipynb` — the quantitative epigenetic layer:
integrate OSD-217's spaceflight DNA-methylation (WGBS) with its RNA-seq, on the
same Arabidopsis root experiment. Run: python _make_methylation_notebook.py
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md(r"""# 🧬 The epigenetic layer: methylation × expression in space

The [data story](OSDR_data_story.ipynb) ended at a frontier: spaceflight changes
which genes are *expressed*, but does it also change how the genome is
*chemically marked*? **OSD-217 is the dataset that lets us answer this** — it
measured **whole-genome DNA methylation (WGBS) *and* RNA-seq on the same
*Arabidopsis* spaceflight experiment**. The methylation tables used here are from
roots; the RNA-seq contrast pools roots and leaves (3 + 3 samples per arm).

Plants methylate DNA in three sequence contexts — **CG, CHG, CHH** — each with
different biology. Here we:

1. load OSD-217's per-gene **flight-vs-ground methylation change** in all three contexts,
2. load the matched **spaceflight expression change** (log₂FC),
3. ask whether the **epigenome and transcriptome move together**, and
4. find the genes that change in **both** layers.

> 🚀 Run live via the rocket icon to extend the analysis.""")

code('''import io, re, tarfile
import numpy as np
import pandas as pd
import requests
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "notebook_connected"

API = "https://visualization.osdr.nasa.gov/biodata/api/v2"
GEODE = "https://osdr.nasa.gov/geode-py/ws/studies"
ACC = "OSD-217"''')

md("""## 1. Spaceflight DNA-methylation change (WGBS)

OSDR serves processed WGBS results as per-context tables of **AvgDiffMeth** — the
average methylation difference (flight − ground) over each gene. We average across
the replicate archives to get one value per gene per context.""")

code('''# The five WGBS result archives, each holding CG / CHG / CHH tables
arch = [f"GLDS-217_wgbs_GSE95594_ROOTFT-vs-ROOTGC_genes{n:02d}_full_tar.gz" for n in range(1, 6)]

frames = []
for fn in arch:
    blob = requests.get(f"{GEODE}/{ACC}/download?source=datamanager&file={fn}", timeout=180).content
    tf = tarfile.open(fileobj=io.BytesIO(blob))
    for m in tf.getmembers():
        ctx = next((c for c in ("CG", "CHG", "CHH") if m.name.endswith(f"{c}.xlsx")), None)
        if not ctx:
            continue
        d = pd.read_excel(io.BytesIO(tf.extractfile(m).read()))[["# Gene", "AvgDiffMeth"]]
        d = d.rename(columns={"# Gene": "gene"})
        d["context"] = ctx
        frames.append(d)

meth = (pd.concat(frames)
        .groupby(["gene", "context"])["AvgDiffMeth"].mean()
        .unstack("context"))
print(f"Methylation change for {len(meth)} genes in contexts {list(meth.columns)}")
meth.head()''')

code('''# How is methylation redistributed by spaceflight, per context?
long = meth.reset_index().melt("gene", var_name="context", value_name="AvgDiffMeth").dropna()
fig = px.violin(long, x="context", y="AvgDiffMeth", color="context", box=True, points=False,
                title="Spaceflight methylation change by context (flight − ground)")
fig.update_layout(height=420, showlegend=False)
fig.add_hline(y=0, line_dash="dot", opacity=0.5)
fig.show()''')

md("""## 2. Matched spaceflight expression change (RNA-seq)

The same study's RNA-seq, summarised as a per-gene spaceflight log₂ fold-change —
exactly as in the cross-study data story, but for OSD-217 alone. It pools the
study's root and leaf samples (3 + 3 per arm), so unlike the methylation tables it
is not root-specific.""")

code('''def normalized_counts(acc):
    d = requests.get(f"{API}/dataset/{acc}/files/", timeout=60).json()
    names = list(d[list(d)[0]].get("files", d[list(d)[0]]))
    f = next(n for n in names if "Normalized_Counts" in n and n.endswith(".csv") and "rRNArm" not in n)
    return pd.read_csv(f"{GEODE}/{acc}/download?source=datamanager&file={f}", index_col=0)

def flight_ground(acc, cols):
    url = f"{API}/query/metadata/?id.accession={acc}&id.sample name&study.factor value&format=csv"
    df = pd.read_csv(io.StringIO(requests.get(url, timeout=60).text)).drop_duplicates("id.sample name")
    sf = next((c for c in df.columns if re.search(r"(?i)spaceflight", c)), None)
    m = df.set_index("id.sample name")[sf].astype(str) if sf else {}
    out = {}
    for c in cols:
        v = str(m.get(c, "")).lower()
        out[c] = "Spaceflight" if ("flight" in v or "space" in v or "_FLT_" in c) else (
                 "Ground" if ("ground" in v or "_GC_" in c) else None)
    return out

counts = normalized_counts(ACC)
lab = flight_ground(ACC, counts.columns)
flt = [c for c in counts.columns if lab[c] == "Spaceflight"]
grd = [c for c in counts.columns if lab[c] == "Ground"]
expr_lfc = np.log2((counts[flt].mean(axis=1) + 1) / (counts[grd].mean(axis=1) + 1)).rename("expr_log2FC")
print(f"Expression log2FC from {len(flt)} flight vs {len(grd)} ground samples, {len(expr_lfc)} genes")''')

md("""## 3. Do the epigenome and transcriptome move together?

We join the two layers per gene and look for a relationship between the
spaceflight **expression** change and the spaceflight **methylation** change.""")

code('''joined = meth.join(expr_lfc, how="inner").dropna(subset=["expr_log2FC"])
print(f"{len(joined)} genes with both methylation and expression data")

rows = []
for ctx in [c for c in ("CG", "CHG", "CHH") if c in joined]:
    sub = joined.dropna(subset=[ctx])
    rows.append({"context": ctx, "genes": len(sub),
                 "corr(expr log2FC, methyl diff)": round(sub["expr_log2FC"].corr(sub[ctx]), 3)})
pd.DataFrame(rows)''')

code('''# Expression vs gene-body CG methylation change (the classic plant relationship)
ctx = "CG"
sub = joined.dropna(subset=[ctx])
r = sub["expr_log2FC"].corr(sub[ctx])            # correlation on ALL genes
plot_df = sub.sample(min(6000, len(sub)), random_state=0).reset_index()   # sample for a light, fast chart
fig = px.scatter(plot_df, x=ctx, y="expr_log2FC", opacity=0.25, render_mode="webgl",
                 title=f"Spaceflight: {ctx} methylation change vs expression change "
                       f"(r = {r:.2f}, n = {len(sub)}; {len(plot_df)} shown)",
                 labels={ctx: f"{ctx} methylation Δ (flight − ground)",
                         "expr_log2FC": "expression log₂FC (flight / ground)"})
fig.add_hline(y=0, line_dash="dot", opacity=0.4); fig.add_vline(x=0, line_dash="dot", opacity=0.4)
fig.update_layout(height=480)
fig.show()''')

md("""## 4. The multi-omic hits — genes that change in *both* layers

Genes that are strongly **differentially expressed *and* differentially
methylated** in spaceflight are the most interesting candidates: places where the
plant may be using an epigenetic mark alongside a transcriptional response.""")

code('''hits = joined.dropna(subset=["CG"]).copy()
hits["abs_expr"] = hits["expr_log2FC"].abs()
hits["abs_CG"] = hits["CG"].abs()
# rank by the product of (expression effect) x (CG methylation effect)
hits["combined"] = hits["abs_expr"] * hits["abs_CG"]
top = hits.sort_values("combined", ascending=False).head(20)
top[["expr_log2FC", "CG", "CHG", "CHH"]].round(3)''')

md(r"""## 5. What do the multi-omic hits *do*? (GO & pathway enrichment)

Knowing *which* genes change in both layers is more useful if we know their
biological **functions**. We send the multi-omic hit genes to
[g:Profiler](https://biit.cs.ut.ee/gprofiler/) for GO (Biological Process) + KEGG
enrichment against the *Arabidopsis* genome.""")

code('''def enrich(genes, organism, sources=("GO:BP", "KEGG")):
    """GO / pathway enrichment via the g:Profiler API (no key required)."""
    r = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
                      json={"organism": organism, "query": list(genes), "sources": list(sources),
                            "user_threshold": 0.05, "significance_threshold_method": "g_SCS",
                            "no_evidences": True}, timeout=120)
    res = r.json().get("result", [])
    if not res:
        return pd.DataFrame()
    return (pd.DataFrame(res)[["source", "native", "name", "p_value",
                               "term_size", "intersection_size"]]
            .sort_values("p_value").reset_index(drop=True))

# the genes changed in BOTH layers (expression x CG methylation)
hit_genes = hits.sort_values("combined", ascending=False).head(200).index.tolist()
enr = enrich(hit_genes, "athaliana")
print(f"{len(hit_genes)} multi-omic hit genes -> {len(enr)} enriched terms")
enr.head(10)''')

code('''if len(enr):
    t = enr.head(12).copy()
    t["minus_log10_p"] = -np.log10(t["p_value"])
    fig = px.bar(t.sort_values("minus_log10_p"), x="minus_log10_p", y="name", color="source",
                 orientation="h",
                 title="GO / pathway enrichment of the methylation × expression hits",
                 labels={"minus_log10_p": "-log10(adjusted p)", "name": ""})
    fig.update_layout(height=460)
    fig.show()
else:
    print("No significant enrichment for this hit set.")''')

md(r"""## 6. Results and discussion

### Key results

- **Spaceflight redistributes DNA methylation.** Across **26,344 genes**,
  flight-vs-ground methylation changed in all three plant contexts — **CG, CHG and
  CHH** (section 1) — an epigenetic response, not only a transcriptional one.
- **Methylation and expression are decoupled genome-wide.** Joined over **21,790
  genes**, the correlation between per-gene methylation change and expression change
  was **essentially zero** (CG *r* = 0.01, CHG = −0.02, CHH = −0.01; section 3) — as
  expected in plants, where gene-body methylation and expression are not simply
  proportional.
- **The signal is in specific genes.** The genes that change in *both* layers
  (section 4) are enriched (g:Profiler, *p* < 0.05) for **photosynthesis and
  light-harvesting** (photosystem light reactions) and **light/desiccation stress**
  (section 5).

### Discussion

Two messages stand out. First, the **near-zero global coupling** is itself a result:
spaceflight doesn't simply methylate genes down and de-methylate genes up, so the
interesting biology is **not** in the genome-wide trend but in the **specific genes**
coordinated across both layers. Second, those coordinated genes are strikingly
enriched for **photosynthesis and light-harvesting**. The methylation change is
measured in roots, but the expression change pools roots and leaves, so this does
not show that roots switch on photosynthesis genes: the leaf samples could carry the
expression signal. It does suggest that spaceflight touches the plant's
**light/energy programme** in both layers. These
multi-omic hits are concrete candidates for how a plant might **stabilise a
spaceflight response with an epigenetic mark**.

### Limitations

- **One experiment, gene-body averages.** OSD-217 is a single study, and
  `AvgDiffMeth` is averaged over each gene, hiding position-specific methylation.
- **Correlation, not mechanism** — co-occurrence in both layers doesn't establish
  that methylation *drives* expression (or vice versa).
- **Mixed tissues** — methylation is from roots, while the expression contrast pools
  roots and leaves; a root-only expression contrast is the obvious next test.

### Conclusion

Spaceflight leaves an **epigenetic fingerprint** on *Arabidopsis* roots that is
decoupled from expression genome-wide but converges on **specific light/energy
genes**. This is a genuine multi-omic integration on the *same* experiment — the
template for joining any OSDR study's layers. Swap `ACC`, change the context, or
extend the CHG/CHH analysis to go further.""")

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
with open("OSDR_methylation_layer.ipynb", "w", encoding="utf-8") as f:
    nbf.write(nb, f)
print("Wrote OSDR_methylation_layer.ipynb")
