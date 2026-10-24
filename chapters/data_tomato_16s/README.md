# OSD-766 tomato 16S genus table — provenance

Built 2026-10-24 for the [tomato microbiome pipeline chapter](../OSDR_tomato_microbiome_pipeline.ipynb).
OSDR serves OSD-766 (VEG-05) as raw amplicon reads only; this folder turns the 16S reads into
a genus table.

## Inputs

- **Reads:** OSD-766 16S R1/R2 FASTQ for the 44 per-plant samples (root, adventitious root,
  leaf), downloaded from OSDR by `download_reads.py`, 88 files, about 330 MB. They are not
  committed: re-download them with the script.
- **Reference:** `silva_nr99_v138.1_train_set.fa.gz`, SILVA 138.1 formatted for DADA2,
  [Zenodo 10.5281/zenodo.4587955](https://doi.org/10.5281/zenodo.4587955), md5
  `6b41db7139834c71171f8ce5b5918fc6` (matches the Zenodo record).
- **Sample factors:** spaceflight and light regimen from OSDR's metadata API, fetched by
  `make_genus_table.py`.

## Steps

| Script | What it does |
|---|---|
| `download_reads.py <dir>` | Fetches the 88 FASTQ files (per-plant samples only, so it excludes water-locker "root-mat" and pooled samples). |
| `inspect_reads.py <dir>` | Checks read starts for common 16S primers (none found), read length and quality by position. |
| `run_dada2.R <reads> <silva> <work> .` | DADA2 1.40.0 on R 4.6.1: `truncLen = c(220, 180)`, `maxEE = c(2, 2)`, consensus chimera removal, SILVA taxonomy with a 50% bootstrap threshold, chloroplast and mitochondrial ASVs removed. |
| `make_genus_table.py` | Builds `samples.csv` and `../tomato_genus_by_sample.csv`, keeping samples with at least 1,000 bacterial reads. |

Re-running `run_dada2.R` reproduces `genus_counts.csv` byte for byte (fixed seeds).

## Outputs

| File | Contents |
|---|---|
| `read_tracking.csv` | Reads per sample at each step: input, filtered, denoised, merged, non-chimeric, organelle, bacterial. |
| `genus_counts.csv` | Bacterial read counts, 105 genera × 44 samples (`Unassigned` = no genus at 50% support). |
| `asv_taxonomy.csv` | The 206 bacterial ASVs: sequence, SILVA ranks, best-guess genus and its bootstrap, total reads. |
| `asv_counts.csv` | Bacterial read counts per ASV and sample. |
| `samples.csv` | Per sample: OSDR spaceflight and light factors, tissue, plant, read depths, and whether it passes the depth filter. |
| `../tomato_genus_by_sample.csv` | Relative abundance (% of bacterial reads), 105 genera × 21 samples that pass the filter. |

## Things to know

- **Region:** the deposited reads carry no primer sequence. Most merged amplicons are 253 bp,
  the length of the V4 region between the 515F and 806R primers.
- **Leaves:** about 99% of leaf reads are tomato chloroplast or mitochondrial DNA, so all but
  one leaf sample fail the 1,000-bacterial-read filter.
- **Rhizobium clade:** the main flight Rhizobiaceae variant (ASV4) matches the
  *Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium* clade at only 46% bootstrap support,
  so it is unnamed at genus level. Genus-level results for that clade depend on the
  taxonomy threshold.
