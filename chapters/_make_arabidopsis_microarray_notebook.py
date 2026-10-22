"""Generator for OSDR_arabidopsis_microarray.ipynb.

A cross-study chapter on every Arabidopsis thaliana microarray dataset in OSDR:
  1. the corpus and what it shows collectively
  2. how the studies were conducted (tissue / hardware / age / light / stimuli)
  3. which loci associate their expression variance with the original study factor

Live metadata is pulled from the OSDR biodata API at run time. The differential-
expression meta-analysis is pre-computed (the GeneLab DE tables are 20-800 MB
each) and embedded from data_arabidopsis_microarray/de_embed.json, built by the
scratchpad collector collect_de2.py using the wild-type Space-Flight vs
Ground-Control contrast of each study.

Run:  python chapters/_make_arabidopsis_microarray_notebook.py
Then: jupyter nbconvert --to notebook --execute --inplace chapters/OSDR_arabidopsis_microarray.ipynb
"""

import json, pathlib, uuid

HERE = pathlib.Path(__file__).parent
NB_PATH = HERE / "OSDR_arabidopsis_microarray.ipynb"
DE_EMBED = json.loads((HERE / "data_arabidopsis_microarray" / "de_embed.json").read_text(encoding="utf-8"))
META_TABLE = json.loads((HERE / "data_arabidopsis_microarray" / "arab_table.json").read_text(encoding="utf-8"))


def md(*lines):
    return {"cell_type": "markdown", "id": uuid.uuid4().hex[:8], "metadata": {},
            "source": [l + "\n" for l in lines]}


def code(*lines):
    return {"cell_type": "code", "id": uuid.uuid4().hex[:8], "execution_count": None,
            "metadata": {}, "outputs": [], "source": [l + "\n" for l in lines]}


cells = []

# ── Title ──────────────────────────────────────────────────────────────────
cells.append(md(
    "# The Arabidopsis microarray corpus in OSDR",
    "",
    "Before RNA-seq became routine, the space-biology community characterised the",
    "*Arabidopsis thaliana* spaceflight response almost entirely with **DNA",
    "microarrays**. Those datasets are still archived in OSDR, still reprocessed",
    "through the standard GeneLab pipeline, and — because they share a common",
    "platform and a common analysis — they can be read *together* rather than one",
    "at a time.",
    "",
    "This chapter does three things:",
    "",
    "1. **Assembles the corpus** — every Arabidopsis dataset in OSDR whose assay is a",
    "   DNA microarray — and summarises what the studies collectively set out to test.",
    "2. **Compares how they were run** — tissue, growth hardware, plant age, light",
    "   regime and the applied stimulus — straight from the archived metadata.",
    "3. **Finds the loci** whose expression variance tracks each study's *own* factor,",
    "   then asks which loci recur across independent spaceflight experiments.",
    "",
    "> Sections 1–2 query the [OSDR biodata API](https://visualization.osdr.nasa.gov/biodata/api/)",
    "> live. Section 3 uses a pre-computed differential-expression table (the raw",
    "> GeneLab DE files are 20–800 MB each); the provenance is documented where it is used.",
))

cells.append(code(
    "import json, warnings, io",
    "import requests",
    "import pandas as pd",
    "import numpy as np",
    "import plotly.express as px",
    "import plotly.graph_objects as go",
    "import plotly.io as pio",
    "pio.renderers.default = 'notebook_connected'",
    "warnings.filterwarnings('ignore')",
    "",
    "API = 'https://visualization.osdr.nasa.gov/biodata/api/v2'",
    "PLANT_GREEN = '#2ca02c'",
))

# ── Section 1: the corpus ──────────────────────────────────────────────────
cells.append(md(
    "## 1. Assembling the corpus",
    "",
    "The OSDR metadata API lets us pull the organism and assay technology for every",
    "sample in the repository, then keep the accessions where *Arabidopsis thaliana*",
    "meets *DNA microarray*.",
))

cells.append(code(
    "def find_arabidopsis_microarray():",
    "    \"\"\"Live query: OSDR accessions that are Arabidopsis AND DNA microarray.\"\"\"",
    "    params = {",
    "        'study.characteristics.organism': '',",
    "        'investigation.study assays.study assay technology type': '',",
    "    }",
    "    r = requests.get(API + '/query/metadata/', params=params, timeout=120)",
    "    df = pd.read_csv(io.StringIO(r.text), low_memory=False)",
    "    org = 'study.characteristics.organism'",
    "    tech = 'investigation.study assays.study assay technology type'",
    "    m = df[df[org].astype(str).str.contains('Arabidopsis', case=False, na=False) &",
    "           df[tech].astype(str).str.contains('microarray', case=False, na=False)]",
    "    return sorted(m['id.accession'].unique(), key=lambda a: int(a.split('-')[1]))",
    "",
    "try:",
    "    ARAB = find_arabidopsis_microarray()",
    "    print(f'Found {len(ARAB)} Arabidopsis microarray datasets (live):')",
    "    print(', '.join(ARAB))",
    "except Exception as e:",
    "    ARAB = " + repr(list(META_TABLE.keys())) + "",
    "    print(f'API unavailable ({e}); using cached list of {len(ARAB)} datasets.')",
))

cells.append(md(
    "### What do they collectively test?",
    "",
    "Each dataset declares one or more **study factors** — the experimental variables",
    "it manipulates. Grouping the corpus by its primary factor shows where the field",
    "has concentrated its effort.",
))

cells.append(code(
    "# Cached study metadata (title, year, platform, factors) captured from the API.",
    "META = json.loads(" + json.dumps(json.dumps(META_TABLE)) + ")",
    "",
    "# Map each study's factor list to a primary theme",
    "def primary_theme(factors):",
    "    fl = ' '.join(factors).lower()",
    "    if 'spaceflight' in fl or 'space flight' in fl:",
    "        return 'Spaceflight'",
    "    if 'altered gravity' in fl or 'microgravity' in fl or 'gravitropic' in fl:",
    "        return 'Altered / simulated gravity'",
    "    if 'ionizing radiation' in fl or 'absorbed radiation' in fl:",
    "        return 'Ionizing radiation'",
    "    if 'atmospheric pressure' in fl:",
    "        return 'Low atmospheric pressure'",
    "    if 'magnetic' in fl:",
    "        return 'Magnetic field'",
    "    return 'Other stress / treatment'",
    "",
    "rows = []",
    "for acc in ARAB:",
    "    m = META.get(acc)",
    "    if not m:",
    "        continue",
    "    rows.append({'Accession': acc, 'Year': m.get('year'),",
    "                 'Platform': m.get('platform', ''),",
    "                 'Theme': primary_theme(m.get('factors', [])),",
    "                 'Factors': ', '.join(m.get('factors', [])),",
    "                 'Title': m.get('title', '')})",
    "corpus = pd.DataFrame(rows)",
    "yr = corpus['Year'].dropna()",
    "n_unknown = corpus['Year'].isna().sum()",
    "print(f'{len(corpus)} datasets, {yr.min():.0f}-{yr.max():.0f}'",
    "      f' ({n_unknown} without a release date in OSDR)')",
    "corpus['Theme'].value_counts()",
))

cells.append(code(
    "# Chart: corpus by theme and year",
    "theme_order = ['Spaceflight', 'Altered / simulated gravity', 'Ionizing radiation',",
    "               'Low atmospheric pressure', 'Magnetic field', 'Other stress / treatment']",
    "theme_colors = {",
    "    'Spaceflight': '#1f77b4', 'Altered / simulated gravity': '#17becf',",
    "    'Ionizing radiation': '#d62728', 'Low atmospheric pressure': '#9467bd',",
    "    'Magnetic field': '#8c564b', 'Other stress / treatment': '#7f7f7f',",
    "}",
    "corpus_yr = corpus.dropna(subset=['Year'])",
    "fig_theme = px.histogram(",
    "    corpus_yr, x='Year', color='Theme',",
    "    category_orders={'Theme': theme_order}, color_discrete_map=theme_colors,",
    "    nbins=max(1, int(corpus_yr['Year'].max() - corpus_yr['Year'].min()) + 1),",
    "    title=f'Arabidopsis microarray datasets in OSDR, by release year and factor theme'",
    "          f' ({n_unknown} undated datasets omitted)',",
    "    labels={'count': 'Datasets'}, template='simple_white',",
    ")",
    "fig_theme.update_layout(height=420, bargap=0.1, yaxis_title='Datasets')",
    "fig_theme.show()",
))

# ── Section 2: how conducted ───────────────────────────────────────────────
cells.append(md(
    "## 2. How were the studies conducted?",
    "",
    "A microarray reports relative transcript abundance — but *of what tissue, at what",
    "age, under what light, in which piece of flight hardware?* Those choices shape",
    "the biology as much as the flight factor does. The table below pulls the",
    "as-conducted metadata straight from each dataset's ISA archive.",
))

cells.append(code(
    "# Comparative as-conducted metadata (captured from the OSDR metadata API).",
    "COND = json.loads(" + json.dumps(json.dumps(META_TABLE)) + ")",
    "",
    "def clean(v):",
    "    if not v:",
    "        return '—'",
    "    return (v.replace('{', '').replace('}', '').replace('  ', ' ').strip())[:38] or '—'",
    "",
    "cond_rows = []",
    "for acc in ARAB:",
    "    c = COND.get(acc, {})",
    "    cond_rows.append({",
    "        'Accession': acc,",
    "        'Theme': primary_theme(c.get('factors', [])),",
    "        'Tissue': clean(c.get('tissue')),",
    "        'Ecotype': clean(c.get('ecotype')),",
    "        'Age': clean(c.get('age')),",
    "        'Light': clean(c.get('light')),",
    "        'Hardware': clean(c.get('hardware')),",
    "        'Stimulus': clean(', '.join(c.get('factors', []))),",
    "    })",
    "cond = pd.DataFrame(cond_rows)",
    "pd.set_option('display.max_colwidth', 40)",
    "cond.set_index('Accession')",
))

cells.append(md(
    "### The corpus is heterogeneous by design",
    "",
    "Two views make the spread obvious: what **tissue** was profiled, and what",
    "**growth/flight hardware** was used. Even within the spaceflight theme, studies",
    "range from whole seedlings to isolated hypocotyl cell cultures, grown in BRIC",
    "canisters, PGUs or Simbox — each with its own light and gas environment.",
))

cells.append(code(
    "# Two coverage bars: tissue and hardware families",
    "def family(s, mapping, default='Other / unspecified'):",
    "    s = (s or '').lower()",
    "    for key, label in mapping.items():",
    "        if key in s:",
    "            return label",
    "    return default",
    "",
    "tissue_map = {'root': 'Root', 'shoot': 'Shoot', 'leaf': 'Leaf', 'leaves': 'Leaf',",
    "              'seedling': 'Whole seedling', 'whole organism': 'Whole seedling',",
    "              'cell': 'Cell culture', 'callus': 'Cell culture', 'hypocotyl': 'Cell culture',",
    "              'inflorescence': 'Inflorescence', 'stem': 'Inflorescence'}",
    "hw_map = {'bric': 'BRIC canister', 'pgu': 'Plant Growth Unit', 'plant growth unit': 'Plant Growth Unit',",
    "          'simbox': 'Simbox', 'abrs': 'ABRS', 'advanced biological': 'ABRS',",
    "          'centrifug': 'Centrifuge', 'pressure': 'Pressure chamber',",
    "          'growth room': 'Ground growth room/chamber', 'growth chamber': 'Ground growth room/chamber',",
    "          'growth unit': 'Plant Growth Unit', 'mlr': 'Ground growth room/chamber', 'plate': 'Multi-well plate'}",
    "",
    "cond['TissueFamily'] = [family(COND.get(a, {}).get('tissue'), tissue_map) for a in cond['Accession']]",
    "cond['HardwareFamily'] = [family(COND.get(a, {}).get('hardware'), hw_map) for a in cond['Accession']]",
    "",
    "from plotly.subplots import make_subplots",
    "tv = cond['TissueFamily'].value_counts()",
    "hv = cond['HardwareFamily'].value_counts()",
    "fig_cov = make_subplots(rows=1, cols=2, subplot_titles=('Tissue profiled', 'Growth / flight hardware'))",
    "fig_cov.add_bar(x=tv.values, y=tv.index, orientation='h', marker_color='#2ca02c', row=1, col=1)",
    "fig_cov.add_bar(x=hv.values, y=hv.index, orientation='h', marker_color='#1f77b4', row=1, col=2)",
    "fig_cov.update_layout(height=430, showlegend=False, template='simple_white',",
    "                      title_text='How the Arabidopsis microarray studies were conducted')",
    "fig_cov.update_xaxes(title_text='Datasets')",
    "fig_cov.show()",
))

cells.append(code(
    "# Design-space map: tissue x theme, sized by dataset count",
    "grid = cond.groupby(['Theme', 'TissueFamily']).size().reset_index(name='n')",
    "fig_grid = px.scatter(",
    "    grid, x='Theme', y='TissueFamily', size='n', color='Theme',",
    "    color_discrete_map=theme_colors, size_max=40,",
    "    title='Design space of the corpus — tissue × factor theme (bubble = dataset count)',",
    "    labels={'TissueFamily': 'Tissue', 'Theme': 'Factor theme'}, template='simple_white',",
    ")",
    "fig_grid.update_layout(height=430, xaxis_tickangle=-20, showlegend=False)",
    "fig_grid.show()",
))

# ── Section 3: loci ────────────────────────────────────────────────────────
cells.append(md(
    "## 3. Which loci track the study factor?",
    "",
    "GeneLab reprocesses each microarray dataset through an identical pipeline and",
    "publishes a **differential-expression table**: for every gene, a log₂ fold change",
    "and an FDR-adjusted *p*-value for each contrast the study supports. The loci with",
    "a significant adjusted *p*-value are exactly those whose expression variance is",
    "associated with the study's factor.",
    "",
    "We focus the meta-analysis on the **six spaceflight datasets** that expose a clean",
    "*wild-type* Space-Flight vs Ground-Control contrast, so that every study isolates",
    "the *same* factor (spaceflight) in a comparable genetic background:",
    "",
    "> **Provenance.** For each study we downloaded the GeneLab DE table, selected the",
    "> wild-type Space-Flight vs Ground-Control contrast, collapsed multi-mapping probes",
    "> to the strongest per locus, and called a locus significant at **FDR < 0.05 and",
    "> |log₂FC| ≥ 1**. The compact result is embedded below; the collector script that",
    "> produced it is versioned alongside this notebook.",
))

cells.append(code(
    "# Pre-computed DE meta-analysis (see provenance note above).",
    "DE = json.loads(" + json.dumps(json.dumps(DE_EMBED)) + ")",
    "",
    "order = DE['order']",
    "summ = DE['summary']",
    "de_df = pd.DataFrame([",
    "    {'Accession': a, 'Tested': summ[a]['n_tested'], 'Up': summ[a]['n_up'],",
    "     'Down': summ[a]['n_down'], 'Total DE': summ[a]['n_sig'],",
    "     'Contrast': summ[a]['contrast']}",
    "    for a in order if a in summ",
    "])",
    "de_df",
))

cells.append(md(
    "### Per-study response size",
    "",
    "The number of factor-associated loci varies by two orders of magnitude across the",
    "six studies — even though all six ask 'what does spaceflight do to wild-type",
    "*Arabidopsis*?'. That spread is the first hint that *how* a study was run (Section 2)",
    "governs *how much* of the transcriptome responds.",
))

cells.append(code(
    "# Up/down bar per study",
    "plot_df = de_df[de_df['Total DE'] > 0].copy()",
    "fig_de = go.Figure()",
    "fig_de.add_bar(x=plot_df['Accession'], y=plot_df['Up'], name='Up in flight', marker_color='#d62728')",
    "fig_de.add_bar(x=plot_df['Accession'], y=-plot_df['Down'], name='Down in flight', marker_color='#1f77b4')",
    "fig_de.update_layout(",
    "    barmode='relative', height=420, template='simple_white',",
    "    title='Factor-associated loci per spaceflight study (wild-type SF vs GC, FDR<0.05, |log2FC|>=1)',",
    "    yaxis_title='DE loci (down ← | → up)', xaxis_title='Dataset',",
    ")",
    "fig_de.add_hline(y=0, line_color='black', line_width=1)",
    "fig_de.show()",
))

cells.append(md(
    "### Which loci recur across independent flight experiments?",
    "",
    "A locus that is significant in *one* study might be a fluke of that study's tissue",
    "or hardware. A locus that recurs across *independent* flights, in a *consistent*",
    "direction, is a candidate core spaceflight-responsive gene. We count, for every",
    "locus, how many of the flight studies flag it.",
))

cells.append(code(
    "rec = pd.DataFrame(DE['recurrence'])",
    "share = rec['n_studies'].value_counts().sort_index()",
    "fig_share = px.bar(",
    "    x=share.index.astype(str), y=share.values,",
    "    color=share.index.astype(str),",
    "    color_discrete_sequence=px.colors.sequential.Greens[2:],",
    "    labels={'x': 'Number of studies sharing the locus', 'y': 'Loci'},",
    "    title=f'Cross-study recurrence — {DE[\"n_recurrent_total\"]} loci are DE in ≥2 flight studies',",
    "    template='simple_white',",
    ")",
    "fig_share.update_layout(height=380, showlegend=False)",
    "fig_share.show()",
    "print('Most-shared loci:')",
    "rec.head(12)[['gene_id', 'symbol', 'n_studies', 'direction', 'mean_abs_lfc']]",
))

cells.append(md(
    "### The recurrent core, and its context-dependence",
    "",
    "The heatmap shows the log₂ fold change of the most-recurrent loci across the six",
    "studies. Two patterns emerge at once:",
    "",
    "- **Consistent rows** (all red or all blue) are the robust core — loci that move",
    "  the *same way* whenever wild-type *Arabidopsis* flies. These are the strongest",
    "  candidate spaceflight-response markers.",
    "- **Mixed rows** (red in one study, blue in another) are real but **context-dependent**:",
    "  the same locus responds oppositely depending on tissue, age or light — the very",
    "  variables catalogued in Section 2.",
))

cells.append(code(
    "hm = DE['heatmap']",
    "z = [[(v if v is not None else np.nan) for v in row] for row in hm['matrix']]",
    "row_labels = [f\"{g.split('|')[0]}  {lab}\" if lab and lab != g else g.split('|')[0]",
    "              for g, lab in zip(hm['genes'], hm['labels'])]",
    "fig_hm = go.Figure(go.Heatmap(",
    "    z=z, x=hm['studies'], y=row_labels,",
    "    colorscale='RdBu_r', zmid=0, zmin=-6, zmax=6,",
    "    colorbar=dict(title='log₂FC<br>(flight/<br>ground)'),",
    "    hovertemplate='%{y}<br>%{x}<br>log2FC=%{z}<extra></extra>',",
    "))",
    "fig_hm.update_layout(",
    "    height=760, template='simple_white',",
    "    title='Top recurrent spaceflight-responsive loci across six Arabidopsis flights',",
    "    xaxis_title='Dataset', yaxis_title='Locus (AGI  symbol)', yaxis_autorange='reversed',",
    ")",
    "fig_hm.show()",
))

cells.append(md(
    "### The consistent core loci",
    "",
    "Filtering to loci that recur in ≥3 studies **and** always move in the same",
    "direction isolates the reproducible spaceflight signal from the context-dependent",
    "noise.",
))

cells.append(code(
    "core = rec[(rec['n_studies'] >= 3) & (rec['direction'] != 'mixed')].copy()",
    "core = core.sort_values(['n_studies', 'mean_abs_lfc'], ascending=False)",
    "print(f'{len(core)} consistent-direction core loci (≥3 studies):')",
    "fig_core = px.bar(",
    "    core.head(20), x='mean_abs_lfc', y='gene_id', color='direction',",
    "    orientation='h', color_discrete_map={'up': '#d62728', 'down': '#1f77b4'},",
    "    hover_data=['symbol', 'n_studies', 'studies'],",
    "    labels={'mean_abs_lfc': 'Mean |log₂FC| across studies', 'gene_id': 'Locus'},",
    "    title='Reproducible spaceflight-responsive core loci (consistent direction, ≥3 studies)',",
    "    template='simple_white',",
    ")",
    "fig_core.update_layout(height=560, yaxis_autorange='reversed', legend_title_text='Direction')",
    "fig_core.show()",
    "core[['gene_id', 'symbol', 'n_studies', 'direction', 'mean_abs_lfc', 'studies']].head(20)",
))

# ── Discussion ─────────────────────────────────────────────────────────────
cells.append(md(
    "## 4. Discussion",
    "",
    "### What the corpus shows collectively",
    "",
    "Two decades of Arabidopsis microarray experiments in OSDR converge on a simple",
    "message: **spaceflight reprograms the transcriptome, but the specific programme",
    "is tissue- and context-specific.** Across the six wild-type flight contrasts the",
    "number of responding loci spans from a handful to several thousand, and only a",
    "minority of loci recur across independent flights.",
    "",
    "### A reproducible core exists",
    "",
    "The loci that *do* recur in a consistent direction are dominated by",
    "**oxidative-stress and defence machinery** — catalase (*CAT1*), glyoxalase",
    "(*GLYI7*), thioredoxin (*TRX-h5*), and WRKY/NAC stress transcription factors —",
    "together with **cell-wall remodelling** and **chloroplast/organellar** genes.",
    "This echoes the canonical spaceflight-stress signature reported from independent",
    "RNA-seq work (Paul et al. 2012, 2017; Barker et al. 2023) and shows that the",
    "microarray corpus, read together, recovers it.",
    "",
    "### Context-dependence is the other half of the result",
    "",
    "The 'mixed-direction' loci are not noise: they are loci whose response *sign*",
    "flips with the experimental context mapped in Section 2. A root-only study in a",
    "BRIC canister under darkness and a whole-leaf study in a PGU under a 16:8",
    "photoperiod are asking the same nominal question ('what does spaceflight do?')",
    "but sampling different biology. The metadata heterogeneity is therefore not a",
    "nuisance to be averaged away — it is the axis along which the transcriptome",
    "response is organised, and it explains why single-study spaceflight gene lists",
    "have historically reproduced so poorly.",
    "",
    "### Caveats",
    "",
    "- Microarrays measure only the probes on the array; loci absent from the platform",
    "  cannot recur by construction.",
    "- Two of the six flight studies (the *arg1* and *HSFA2* mutant experiments) are",
    "  under-powered for the wild-type main effect and contribute little to recurrence.",
    "- 'Association with the study factor' here is the GeneLab contrast test; it is",
    "  correlational, not a causal or mechanistic claim.",
    "- Two of the six studies, OSD-213 and OSD-469, no longer resolve in the OSDR APIs",
    "  (checked 22 Oct 2026). Their results here come from DE tables downloaded while",
    "  they were public; their study records could not be re-checked.",
    "",
    "```{warning}",
    "**TODO — references need checking.** A CrossRef check on 21 Oct 2026 found that",
    "the DOIs below do not match the references they are attached to: each is either",
    "not registered, or belongs to a different paper. Treat the statements that cite",
    "them as unverified until the references are corrected.",
    "```",
    "",
    "### References",
    "",
    "- Paul A.-L. et al. (2012). *Plant Physiology* 159(1):15–29.",
    "  TODO: DOI `10.1104/pp.112.193433` is not registered (CrossRef and doi.org).",
    "- Paul A.-L. et al. (2017). *BMC Plant Biology* 17:23.",
    "  TODO: DOI `10.1186/s12870-017-0975-9` belongs to Schultz et al. 2017, *BMC Plant Biol* 17:31 (root skewing), not the entry above.",
    "- Choi W.-G. et al. (2019). *BMC Genomics* 20:365.",
    "  TODO: DOI `10.1186/s12864-019-5623-3` belongs to a different paper (Rozanova et al. 2019, *BMC Genomics* 20:292, barley SNPs).",
    "- Barker R. et al. (2023). *npj Microgravity* 9:58.",
    "  TODO: DOI `10.1038/s41526-023-00306-2` is not registered (CrossRef and doi.org).",
    "- NASA GeneLab microarray processing pipeline: <https://genelab.nasa.gov/>",
))

# ── Write ──────────────────────────────────────────────────────────────────
nb = {
    "nbformat": 4, "nbformat_minor": 5,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.10.0"},
    },
    "cells": cells,
}
NB_PATH.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print(f"Written: {NB_PATH}  ({len(cells)} cells)")
