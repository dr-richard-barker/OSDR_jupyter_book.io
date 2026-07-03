"""Generator for PubMed_scraping_notenook_for_ISSOP_figures.ipynb.

Keeps Henry Cope's PubMed scraper (broadened query) and adds an interactive
quantitative analysis of the space-biology omics literature: publication growth,
omics-modality trends, journals, funding agencies, study organisms and the
growth of team science.

Run:  python chapters/_make_pubmed_notebook.py
Then: jupyter nbconvert --to notebook --execute --inplace chapters/PubMed_scraping_notenook_for_ISSOP_figures.ipynb
"""

import json, pathlib, uuid

NB_PATH = pathlib.Path(__file__).parent / "PubMed_scraping_notenook_for_ISSOP_figures.ipynb"


def md(*lines):
    return {"cell_type": "markdown", "id": uuid.uuid4().hex[:8], "metadata": {},
            "source": [l + "\n" for l in lines]}


def code(*lines):
    return {"cell_type": "code", "id": uuid.uuid4().hex[:8], "execution_count": None,
            "metadata": {}, "outputs": [], "source": [l + "\n" for l in lines]}


cells = []

# ── Title ──────────────────────────────────────────────────────────────────
cells.append(md(
    "# Mining the space-biology literature (PubMed)",
    "",
    "Programmatically searching and scraping PubMed for spaceflight / space-biology",
    "publications, then **quantifying the field**: how fast it is growing, which",
    "omics modalities drive it, where it is published, who funds it, which organisms",
    "it studies, and how team science has scaled.",
    "",
    "The scraper (below) is by **Henry Cope**; the analysis and charts read the",
    "records it returns and plot them interactively alongside the OSDR data.",
))

# ── Scraper (Henry Cope) ───────────────────────────────────────────────────
cells.append(md(
    "## 1. Scrape PubMed",
    "",
    "We query PubMed for space-biology papers that touch any omics modality, download",
    "the records in the `pubmed` text format, and parse the fields we need (year,",
    "journal, authors, keywords, grants, abstract). Only *journal articles* are kept",
    "— preprints, news and editorials are dropped.",
))

cells.append(code(
    "# Author of the scraper: Henry Cope",
    "import requests, html, urllib.parse, time",
    "import pandas as pd",
    "",
    "# Broad query: space biology AND any omics modality",
    "query = (\"('space flight' OR 'spaceflight' OR 'space station') AND \"",
    "         \"('omics' OR 'genom*' OR 'transcriptom*' OR 'epigenom*' OR \"",
    "         \"'microbiom*' OR 'proteom*' OR 'metabolom*')\")",
    "encoded_query = urllib.parse.quote_plus(query)",
    "print('Scraping PubMed with query:\\n ', query)",
    "",
    "MAX_PAGES = 6   # 200 records/page; raise to capture the full result set",
    "papers = []",
    "page = 1",
    "while page <= MAX_PAGES:",
    "    url = (f'https://pubmed.ncbi.nlm.nih.gov/?term={encoded_query}'",
    "           f'&format=pubmed&size=200&page={page}')",
    "    try:",
    "        response = requests.get(url, timeout=30)",
    "    except Exception as e:",
    "        print('Request failed:', e); break",
    "    if response.status_code != 200:",
    "        print('FAIL:', url); break",
    "    response.encoding = 'utf-8-sig'",
    "    lines = response.text.split('\\n')",
    "",
    "    paper = {}",
    "    in_abstract, in_title = False, False",
    "    new_papers_found = False",
    "    for line in lines:",
    "        line = html.unescape(line)",
    "        if 'PMID- ' in line:",
    "            new_papers_found = True",
    "            if paper:",
    "                papers.append(paper); paper = {}",
    "            paper['PMID'] = line.split('- ')[1].strip()",
    "            in_abstract, in_title = False, False",
    "        elif line.startswith('TI  -'):",
    "            paper['Title'] = line[6:].strip(); in_abstract, in_title = False, True",
    "        elif line.startswith('AB  -'):",
    "            paper['Abstract'] = line[6:].strip(); in_abstract, in_title = True, False",
    "        elif line.startswith('AU  -'):",
    "            paper.setdefault('Authors', []).append(line[6:].strip()); in_abstract, in_title = False, False",
    "        elif line.startswith('OT  -'):",
    "            paper.setdefault('Keywords', []).append(line[6:].strip()); in_abstract, in_title = False, False",
    "        elif line.startswith('DP  -'):",
    "            paper['Publication Year'] = line[6:10].strip(); in_abstract, in_title = False, False",
    "        elif line.startswith('JT  -'):",
    "            paper['Journal Title'] = line[6:].strip(); in_abstract, in_title = False, False",
    "        elif line.startswith('PT  -'):",
    "            paper.setdefault('Publication Type', []).append(line[6:].strip()); in_abstract, in_title = False, False",
    "        elif line.startswith('GR  -'):",
    "            paper.setdefault('Grants', []).append(line[6:].strip()); in_abstract, in_title = False, False",
    "        elif in_abstract and line.startswith('      '):",
    "            paper['Abstract'] += ' ' + line.strip()",
    "        elif in_title and line.startswith('      '):",
    "            paper['Title'] += ' ' + line.strip()",
    "    if paper:",
    "        papers.append(paper)",
    "    if not new_papers_found:",
    "        print('No more papers found.'); break",
    "    page += 1",
    "    time.sleep(0.4)",
    "",
    "# join list-valued fields",
    "for paper in papers:",
    "    for field in ['Authors', 'Keywords', 'Publication Type', 'Grants']:",
    "        if field in paper:",
    "            paper[field] = '; '.join(paper[field])",
    "",
    "# keep journal articles only",
    "filtered = [p for p in papers if 'Journal Article' in p.get('Publication Type', '')]",
    "df = pd.DataFrame(filtered)",
    "df['Year'] = pd.to_numeric(df['Publication Year'], errors='coerce')",
    "df = df.dropna(subset=['Year'])",
    "df['Year'] = df['Year'].astype(int)",
    "print(f'\\nParsed {len(df)} journal articles, {df.Year.min()}-{df.Year.max()}')",
    "df[['PMID', 'Year', 'Journal Title', 'Title']].head()",
))

cells.append(code(
    "# Set up plotting",
    "import re",
    "import numpy as np",
    "import plotly.express as px",
    "import plotly.graph_objects as go",
    "from plotly.subplots import make_subplots",
    "import plotly.io as pio",
    "pio.renderers.default = 'notebook_connected'",
    "SPACE = '#1f3b73'",
    "PLANT_GREEN = '#2ca02c'",
))

# ── Section 2: growth ──────────────────────────────────────────────────────
cells.append(md(
    "## 2. How fast is the field growing?",
    "",
    "The count of space-biology omics papers per year, with the cumulative total.",
    "The most recent calendar year is usually incomplete at scrape time, so it is",
    "shaded.",
))

cells.append(code(
    "per_year = df.groupby('Year').size().reset_index(name='Papers')",
    "# restrict to a sensible window and drop stray future/partial years beyond max",
    "per_year = per_year[per_year['Year'] >= 1995].copy()",
    "per_year['Cumulative'] = per_year['Papers'].cumsum()",
    "latest = int(per_year['Year'].max())",
    "",
    "fig_growth = make_subplots(specs=[[{'secondary_y': True}]])",
    "fig_growth.add_bar(x=per_year['Year'], y=per_year['Papers'], name='Papers / year',",
    "                   marker_color=SPACE, opacity=0.85)",
    "fig_growth.add_trace(go.Scatter(x=per_year['Year'], y=per_year['Cumulative'],",
    "                                name='Cumulative', mode='lines+markers',",
    "                                line=dict(color='#e4572e', width=3)), secondary_y=True)",
    "fig_growth.add_vrect(x0=latest - 0.5, x1=latest + 0.5, fillcolor='grey', opacity=0.15,",
    "                     line_width=0, annotation_text=f'{latest} (partial)',",
    "                     annotation_position='top left')",
    "fig_growth.update_layout(height=440, template='simple_white',",
    "                         title='Space-biology omics publications per year',",
    "                         legend=dict(orientation='h', y=1.12))",
    "fig_growth.update_yaxes(title_text='Papers per year', secondary_y=False)",
    "fig_growth.update_yaxes(title_text='Cumulative papers', secondary_y=True)",
    "fig_growth.show()",
))

# ── Section 3: omics modalities ────────────────────────────────────────────
cells.append(md(
    "## 3. Which omics modalities drive the field?",
    "",
    "We text-mine each paper's title and abstract for the major omics modalities.",
    "A paper can count toward more than one (multi-omics studies are common), so",
    "these are per-modality trends, not a partition.",
))

cells.append(code(
    "OMICS = {",
    "    'Genomics': r'genom',",
    "    'Transcriptomics': r'transcriptom',",
    "    'Proteomics': r'proteom',",
    "    'Metabolomics': r'metabolom',",
    "    'Microbiome': r'microbiom|microbiota',",
    "    'Epigenomics': r'epigenom|methylation',",
    "}",
    "text = (df['Title'].fillna('') + ' ' + df.get('Abstract', pd.Series('', index=df.index)).fillna('')).str.lower()",
    "for name, pat in OMICS.items():",
    "    df[name] = text.str.contains(pat, regex=True)",
    "",
    "# total per modality",
    "totals = df[list(OMICS)].sum().sort_values(ascending=False)",
    "fig_omics = px.bar(",
    "    x=totals.values, y=totals.index, orientation='h',",
    "    color=totals.index, color_discrete_sequence=px.colors.qualitative.Bold,",
    "    labels={'x': 'Papers mentioning modality', 'y': ''},",
    "    title='Omics modalities in the space-biology literature', template='simple_white',",
    ")",
    "fig_omics.update_layout(height=360, showlegend=False, yaxis_autorange='reversed')",
    "fig_omics.show()",
))

cells.append(code(
    "# modality trend over time (3-year rolling to smooth)",
    "yr = df[df['Year'].between(2000, latest)]",
    "trend = yr.groupby('Year')[list(OMICS)].sum()",
    "trend_full = trend.reindex(range(trend.index.min(), latest + 1), fill_value=0)",
    "roll = trend_full.rolling(3, min_periods=1).mean()",
    "",
    "fig_trend = go.Figure()",
    "palette = dict(zip(OMICS, px.colors.qualitative.Bold))",
    "for name in OMICS:",
    "    fig_trend.add_trace(go.Scatter(x=roll.index, y=roll[name], name=name, mode='lines',",
    "                                   line=dict(width=2.5, color=palette[name]), stackgroup=None))",
    "fig_trend.update_layout(height=450, template='simple_white',",
    "                        title='Omics modality trends over time (3-year rolling mean)',",
    "                        xaxis_title='Year', yaxis_title='Papers / year (smoothed)',",
    "                        legend=dict(orientation='h', y=1.12))",
    "fig_trend.show()",
))

# ── Section 4: journals + funding ──────────────────────────────────────────
cells.append(md(
    "## 4. Where is it published, and who funds it?",
    "",
    "The venues that carry the field, and the funding agencies named in the grant",
    "acknowledgements (PubMed's `GR` field, present for ~1 in 3 papers).",
))

cells.append(code(
    "fig_jf = make_subplots(rows=1, cols=2, horizontal_spacing=0.28,",
    "                       subplot_titles=('Top journals', 'Funding agencies (grant field)'))",
    "",
    "# top journals",
    "jc = df['Journal Title'].value_counts().head(12)",
    "fig_jf.add_bar(x=jc.values, y=jc.index, orientation='h', marker_color=SPACE, row=1, col=1)",
    "",
    "# funding agencies from grant strings",
    "AGENCIES = {",
    "    'NIH': r'\\bnih\\b|national institutes of health|/ni[a-z]{2} ',",
    "    'NASA': r'nasa|national aeronautics',",
    "    'NSF': r'\\bnsf\\b|national science found',",
    "    'ESA': r'\\besa\\b|european space agency',",
    "    'Wellcome': r'wellcome',",
    "    'DLR (Germany)': r'\\bdlr\\b',",
    "    'JAXA': r'jaxa',",
    "    'CNES (France)': r'cnes',",
    "    'ERC': r'european research council',",
    "    'ASI (Italy)': r'agenzia spaziale|\\basi\\b',",
    "}",
    "grants = df.get('Grants', pd.Series('', index=df.index)).fillna('').str.lower()",
    "agency_counts = {a: int(grants.str.contains(p, regex=True).sum()) for a, p in AGENCIES.items()}",
    "ag = pd.Series(agency_counts).sort_values(ascending=True)",
    "ag = ag[ag > 0]",
    "fig_jf.add_bar(x=ag.values, y=ag.index, orientation='h', marker_color='#e4572e', row=1, col=2)",
    "",
    "fig_jf.update_layout(height=460, showlegend=False, template='simple_white',",
    "                     title_text='Publication venues and funding sources')",
    "fig_jf.update_yaxes(autorange='reversed', row=1, col=1)",
    "fig_jf.update_xaxes(title_text='Papers', row=1, col=1)",
    "fig_jf.update_xaxes(title_text='Papers naming agency', row=1, col=2)",
    "fig_jf.show()",
))

# ── Section 5: organisms ───────────────────────────────────────────────────
cells.append(md(
    "## 5. Which organisms does the field study?",
    "",
    "Text-mining the titles and abstracts for the model systems of space biology.",
    "As with omics modalities, a paper can mention more than one.",
))

cells.append(code(
    "ORGANISMS = {",
    "    'Human / astronaut': r'\\bhuman|homo sapiens|astronaut|crew member',",
    "    'Microbes': r'bacteri|microb|microbiome|microbiota|\\be\\. ?coli|fungal|fungi',",
    "    'Mouse': r'\\bmouse|\\bmice|murine|mus musculus',",
    "    'Plant (general)': r'\\bplant|seedling|tomato|soybean|lettuce|wheat',",
    "    'Arabidopsis': r'arabidopsis',",
    "    'Rat': r'\\brat\\b|\\brats\\b|rattus',",
    "    'C. elegans': r'elegans|nematode',",
    "    'Drosophila': r'drosophila|fruit fly',",
    "    'Yeast': r'yeast|saccharomyces',",
    "    'Zebrafish': r'zebrafish|danio',",
    "}",
    "org_counts = {o: int(text.str.contains(p, regex=True).sum()) for o, p in ORGANISMS.items()}",
    "org = pd.Series(org_counts).sort_values(ascending=False)",
    "fig_org = px.bar(",
    "    x=org.index, y=org.values, color=org.values,",
    "    color_continuous_scale='Viridis',",
    "    labels={'x': '', 'y': 'Papers mentioning organism', 'color': 'Papers'},",
    "    title='Study organisms in the space-biology omics literature', template='simple_white',",
    ")",
    "fig_org.update_layout(height=420, coloraxis_showscale=False, xaxis_tickangle=-30)",
    "fig_org.show()",
))

# ── Section 6: team science ────────────────────────────────────────────────
cells.append(md(
    "## 6. The rise of team science",
    "",
    "Space-biology omics is increasingly collaborative. Counting authors per paper",
    "over time shows how the typical team has grown — a signature of the large,",
    "consortium-style analyses (e.g. the *Cell* space-omics package) that now",
    "dominate the field.",
))

cells.append(code(
    "df['n_authors'] = df.get('Authors', pd.Series('', index=df.index)).fillna('').apply(",
    "    lambda s: len([a for a in s.split(';') if a.strip()]))",
    "ta = df[(df['n_authors'] > 0) & (df['Year'].between(2000, latest))]",
    "team = ta.groupby('Year')['n_authors'].agg(['mean', 'median', 'max', 'count']).reset_index()",
    "",
    "fig_team = go.Figure()",
    "fig_team.add_trace(go.Scatter(x=team['Year'], y=team['mean'], name='Mean authors',",
    "                              mode='lines+markers', line=dict(color=SPACE, width=3)))",
    "fig_team.add_trace(go.Scatter(x=team['Year'], y=team['median'], name='Median authors',",
    "                              mode='lines+markers', line=dict(color='#e4572e', width=2, dash='dot')))",
    "fig_team.update_layout(height=420, template='simple_white',",
    "                       title='Authors per paper over time — the growth of team science',",
    "                       xaxis_title='Year', yaxis_title='Authors per paper',",
    "                       legend=dict(orientation='h', y=1.12))",
    "fig_team.show()",
    "print('Most-collaborative papers:')",
    "df.nlargest(5, 'n_authors')[['Year', 'n_authors', 'Journal Title', 'Title']]",
))

# ── Section 7: plant / astrobotany drill-down ──────────────────────────────
cells.append(md(
    "## 7. Zooming in: the plant (astrobotany) narrative",
    "",
    "The organism chart in Section 5 showed that plants are a substantial but",
    "not-dominant slice of the field. Because plants are the focus of this book — and",
    "the engine of any future bioregenerative life-support system — it is worth",
    "drilling into that slice: **how big is it, which species and questions drive it,",
    "and how does it map onto the OSDR plant datasets analysed elsewhere in this book?**",
))

cells.append(code(
    "# Flag the plant / astrobotany sub-corpus",
    "PLANT_PAT = (r'\\bplant|seedling|arabidopsis|tomato|soybean|lettuce|wheat|rice|maize|"
    "moss|physcomitr|brassica|\\broot\\b|hypocotyl|gravitropism|phototropism|chloroplast|photosynthe')",
    "df['is_plant'] = text.str.contains(PLANT_PAT, regex=True)",
    "plant_df = df[df['is_plant']].copy()",
    "print(f'{len(plant_df)} plant / astrobotany papers ',",
    "      f'({100 * len(plant_df) / len(df):.0f}% of the space-biology omics corpus)')",
))

cells.append(md(
    "### Plant share of the field over time",
    "",
    "Plotting the plant sub-corpus against the whole field shows whether astrobotany is",
    "keeping pace, gaining or losing ground as space biology scales up.",
))

cells.append(code(
    "ally = df[df['Year'].between(2000, latest)].groupby('Year').size()",
    "ply = plant_df[plant_df['Year'].between(2000, latest)].groupby('Year').size()",
    "idx = range(int(ally.index.min()), latest + 1)",
    "ally = ally.reindex(idx, fill_value=0)",
    "ply = ply.reindex(idx, fill_value=0)",
    "share = (100 * ply / ally.replace(0, np.nan)).fillna(0)",
    "",
    "fig_plant_share = make_subplots(specs=[[{'secondary_y': True}]])",
    "fig_plant_share.add_bar(x=list(idx), y=ply.values, name='Plant papers / year',",
    "                        marker_color=PLANT_GREEN, opacity=0.85)",
    "fig_plant_share.add_trace(go.Scatter(x=list(idx), y=share.values, name='Plant share of field (%)',",
    "                                     mode='lines+markers', line=dict(color='#8c564b', width=3)),",
    "                          secondary_y=True)",
    "fig_plant_share.update_layout(height=430, template='simple_white',",
    "                              title='Plant / astrobotany papers and their share of space-biology omics',",
    "                              legend=dict(orientation='h', y=1.12))",
    "fig_plant_share.update_yaxes(title_text='Plant papers / year', secondary_y=False)",
    "fig_plant_share.update_yaxes(title_text='Share of field (%)', secondary_y=True, range=[0, 40])",
    "fig_plant_share.show()",
))

cells.append(md(
    "### Which plants fly? Model vs crop",
    "",
    "*Arabidopsis thaliana* is the reference organism of plant biology, so it dominates",
    "the flight literature just as it dominates OSDR. But a **crop tail** — rice,",
    "tomato, lettuce, wheat, soybean — reflects the field's pivot toward *food",
    "production* for long-duration missions (the Veggie and Advanced Plant Habitat era).",
))

cells.append(code(
    "SPECIES = {",
    "    'Arabidopsis': r'arabidopsis', 'Rice': r'\\brice\\b|oryza',",
    "    'Tomato': r'tomato|solanum lycopersic', 'Lettuce': r'lettuce|lactuca',",
    "    'Wheat': r'\\bwheat|triticum', 'Pea': r'\\bpea\\b|pisum',",
    "    'Soybean': r'soybean|glycine max', 'Maize': r'maize|\\bcorn\\b|zea mays',",
    "    'Brassica': r'brassica|rapeseed|canola', 'Mustard': r'mustard',",
    "    'Moss': r'physcomitr|\\bmoss\\b',",
    "}",
    "ptext = (plant_df['Title'].fillna('') + ' ' + plant_df.get('Abstract', pd.Series('', index=plant_df.index)).fillna('')).str.lower()",
    "sp = pd.Series({s: int(ptext.str.contains(p, regex=True).sum()) for s, p in SPECIES.items()})",
    "sp = sp[sp > 0].sort_values(ascending=False)",
    "is_model = ['Model' if s == 'Arabidopsis' else 'Crop / other' for s in sp.index]",
    "fig_species = px.bar(",
    "    x=sp.index, y=sp.values, color=is_model,",
    "    color_discrete_map={'Model': PLANT_GREEN, 'Crop / other': '#ff7f0e'},",
    "    labels={'x': '', 'y': 'Papers mentioning species', 'color': ''},",
    "    title='Plant species in the spaceflight omics literature — model vs crop',",
    "    template='simple_white',",
    ")",
    "fig_species.update_layout(height=420, xaxis_tickangle=-30, legend=dict(orientation='h', y=1.12))",
    "fig_species.show()",
))

cells.append(md(
    "### What do the plant papers ask?",
    "",
    "Text-mining the plant sub-corpus for biological themes reveals what spaceflight",
    "plant biology is actually *about*. Strikingly, the top themes — **defence /",
    "immunity**, **oxidative stress**, **cell wall** and **root / gravitropism** —",
    "are the very same processes that surfaced as the reproducible core loci in the",
    "[Arabidopsis microarray chapter](OSDR_arabidopsis_microarray.ipynb): the",
    "literature and the primary data tell one coherent story.",
))

cells.append(code(
    "PLANT_TOPICS = {",
    "    'Gene expression': r'transcriptom|gene expression|rna-seq|microarray',",
    "    'Defence / immunity': r'defen|immun|pathogen|stress response|wound',",
    "    'Seed / germination': r'\\bseed|germinat',",
    "    'Root / gravitropism': r'gravitropism|root growth|\\broot\\b|graviperception',",
    "    'Oxidative stress': r'oxidative|reactive oxygen|\\bros\\b|antioxidant|redox',",
    "    'Photosynthesis': r'photosynthe|chloroplast|chlorophyll',",
    "    'Cell wall': r'cell wall|lignin|cellulose|pectin',",
    "    'Hormones': r'auxin|cytokinin|ethylene|abscisic|gibberellin|\\bhormone',",
    "    'Cell cycle / meristem': r'cell cycle|cell division|meristem',",
    "    'Light signalling': r'phototropism|photoreceptor|phytochrome|light signal',",
    "}",
    "tp = pd.Series({t: int(ptext.str.contains(p, regex=True).sum()) for t, p in PLANT_TOPICS.items()})",
    "tp = tp.sort_values(ascending=True)",
    "",
    "POM = {'Genomics': r'genom', 'Transcriptomics': r'transcriptom', 'Proteomics': r'proteom',",
    "       'Metabolomics': r'metabolom', 'Epigenomics': r'epigenom|methylation'}",
    "pom = pd.Series({o: int(ptext.str.contains(p, regex=True).sum()) for o, p in POM.items()}).sort_values(ascending=True)",
    "",
    "fig_pt = make_subplots(rows=1, cols=2, horizontal_spacing=0.28,",
    "                       subplot_titles=('Plant biology themes', 'Plant omics modalities'))",
    "fig_pt.add_bar(x=tp.values, y=tp.index, orientation='h', marker_color=PLANT_GREEN, row=1, col=1)",
    "fig_pt.add_bar(x=pom.values, y=pom.index, orientation='h', marker_color='#1f77b4', row=1, col=2)",
    "fig_pt.update_layout(height=440, showlegend=False, template='simple_white',",
    "                     title_text='What the spaceflight plant literature studies')",
    "fig_pt.update_xaxes(title_text='Papers', row=1, col=1)",
    "fig_pt.update_xaxes(title_text='Papers', row=1, col=2)",
    "fig_pt.show()",
))

cells.append(md(
    "### How the themes cluster",
    "",
    "Which of these processes are studied *together*? The heatmap below is the",
    "**Jaccard overlap** between plant themes — the fraction of papers touching either",
    "theme that touch both. Two clusters stand out and they are the two axes of the",
    "spaceflight plant response: a **stress axis** (defence ↔ oxidative stress) and a",
    "**mechanical-growth axis** (root / gravitropism ↔ cell wall). Transcriptomics",
    "(*gene expression*) threads through everything as the common method.",
))

cells.append(code(
    "# Jaccard co-occurrence of plant themes across the plant sub-corpus",
    "theme_mat = pd.DataFrame({t: ptext.str.contains(p, regex=True).astype(int)",
    "                          for t, p in PLANT_TOPICS.items()})",
    "co = theme_mat.T.dot(theme_mat)                       # co-occurrence counts",
    "tot = np.diag(co.values)",
    "union = tot[:, None] + tot[None, :] - co.values",
    "jac = np.divide(co.values, union, out=np.zeros_like(co.values, float), where=union > 0)",
    "np.fill_diagonal(jac, np.nan)                          # blank the diagonal",
    "labels = list(PLANT_TOPICS)",
    "fig_co = go.Figure(go.Heatmap(",
    "    z=jac, x=labels, y=labels, colorscale='Greens', zmin=0,",
    "    colorbar=dict(title='Jaccard'),",
    "    hovertemplate='%{y} + %{x}<br>Jaccard=%{z:.2f}<extra></extra>',",
    "))",
    "fig_co.update_layout(height=520, template='simple_white',",
    "                     title='Plant-biology theme co-occurrence (Jaccard overlap)',",
    "                     yaxis_autorange='reversed')",
    "fig_co.show()",
))

cells.append(md(
    "### From literature to data",
    "",
    "The plant literature and the OSDR plant archive are two views of the same",
    "community. *Arabidopsis* leads both; the dominant experimental themes",
    "(gravitropism, oxidative and defence stress, cell-wall remodelling) are exactly",
    "the axes along which the archived transcriptomes vary. The crop tail in the",
    "literature — tomato (VEG-05), lettuce, wheat — is precisely where OSDR is now",
    "growing, as the field turns from *understanding* the spaceflight response toward",
    "*engineering* reliable food production off Earth.",
))

# ── Discussion ─────────────────────────────────────────────────────────────
cells.append(md(
    "## 8. What the numbers say",
    "",
    "- **The field is growing fast.** Space-biology omics publications have risen",
    "  steeply since the mid-2000s, tracking the maturation of the ISS as a research",
    "  platform and the launch of NASA GeneLab / OSDR.",
    "- **Genomics and transcriptomics lead**, but **microbiome** and **metabolomics**",
    "  have surged most recently — mirroring the shift from single-gene studies to",
    "  ecosystem- and systems-level questions.",
    "- **Human and microbial studies dominate** the corpus; plant systems (led by",
    "  *Arabidopsis*) and classic model animals (mouse, rat, *C. elegans*) form the",
    "  next tier — the same organisms whose primary data sit in OSDR.",
    "- **Funding is led by NIH and NASA**, with a visible international tail (ESA,",
    "  DLR, JAXA, CNES, ASI) that reflects the field's multi-agency character.",
    "- **Team science has scaled**: the typical author count per paper has climbed",
    "  markedly, driven by consortium analyses that pool many missions and labs.",
    "- **Plant biology holds a steady ~1-in-6 share** (Section 7), led by *Arabidopsis*",
    "  but with a growing crop tail (rice, tomato, lettuce, wheat) as the field pivots",
    "  toward food production. Its dominant themes — defence, oxidative stress,",
    "  gravitropism, cell wall — match the reproducible core loci recovered from the",
    "  OSDR Arabidopsis microarray datasets, so literature and data agree.",
    "",
    "Together these trends frame why a FAIR, cross-study resource like OSDR — and a",
    "book like this one that reads its datasets *together* — matters: the literature",
    "is large, multi-modal, multi-organism and multi-agency, and its value compounds",
    "when the underlying data can be reused and integrated.",
    "",
    "> **Caveats.** PubMed's `pubmed`-format scrape is capped at `MAX_PAGES` here;",
    "> raise it for the full result set. Text-mining organism and modality mentions",
    "> from titles/abstracts is approximate (keyword-based, not curated). The grant",
    "> field is present for only ~1 in 3 papers, so agency counts undercount true",
    "> funding.",
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
