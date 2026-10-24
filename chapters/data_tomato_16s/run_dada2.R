# DADA2: OSD-766 16S reads -> ASVs -> SILVA taxonomy -> genus table.
#
# The deposited reads already have their primers removed (checked with inspect_reads.py),
# so there is no cutadapt step. Truncation follows the recipe in the pipeline chapter.
# Chloroplast and mitochondrial ASVs (host organelles) and non-bacterial ASVs are removed
# before collapsing to genus.
#
# Usage: Rscript run_dada2.R <reads_dir> <silva_train_set.fa.gz> <work_dir> <repo_out_dir>
args <- commandArgs(trailingOnly = TRUE)
reads <- args[1]; silva <- args[2]; work <- args[3]; out <- args[4]
suppressPackageStartupMessages(library(dada2))
cat("dada2", as.character(packageVersion("dada2")), "| R", R.version$major, R.version$minor, "\n")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

fnF <- sort(list.files(reads, "_R1_raw\\.fastq\\.gz$", full.names = TRUE))
fnR <- sort(list.files(reads, "_R2_raw\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fnF) == length(fnR), length(fnF) > 0)
# sample names as OSDR lists them, e.g. "113-VEG-05F-SN06-leaf-tom-red-rich_S111_L001"
sample <- sub("_R1_raw\\.fastq\\.gz$", "", sub("^GLDS-672_Amplicon_16S-", "", basename(fnF)))
stopifnot(identical(sample, sub("_R2_raw\\.fastq\\.gz$", "", sub("^GLDS-672_Amplicon_16S-", "", basename(fnR)))))

filtF <- file.path(work, "filtered", paste0(sample, "_F.fastq.gz"))
filtR <- file.path(work, "filtered", paste0(sample, "_R.fastq.gz"))
names(filtF) <- sample; names(filtR) <- sample

set.seed(100)
ft <- filterAndTrim(fnF, filtF, fnR, filtR, truncLen = c(220, 180), maxN = 0,
                    maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE,
                    multithread = TRUE)
rownames(ft) <- sample
keep <- ft[, "reads.out"] > 0 & file.exists(filtF)
cat(sum(keep), "of", length(keep), "samples have reads after filtering\n")

errF <- learnErrors(filtF[keep], multithread = TRUE)
errR <- learnErrors(filtR[keep], multithread = TRUE)
dadaF <- dada(filtF[keep], err = errF, multithread = TRUE)
dadaR <- dada(filtR[keep], err = errR, multithread = TRUE)
mergers <- mergePairs(dadaF, filtF[keep], dadaR, filtR[keep])
seqtab <- makeSequenceTable(mergers)
cat("merged amplicon lengths (bp: number of ASVs):\n"); print(table(nchar(getSequences(seqtab))))
nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
cat(sprintf("non-chimeric: %d ASVs, %.1f%% of merged reads\n", ncol(nochim), 100 * sum(nochim) / sum(seqtab)))

set.seed(100)
# classify once with no cutoff to keep each ASV's best guess, then apply DADA2's default
# 50% bootstrap rule: a rank below 50% support (and every rank under it) becomes NA
tx <- assignTaxonomy(nochim, silva, minBoot = 0, multithread = TRUE, outputBootstraps = TRUE)
best <- tx$tax
tax <- best
tax[tx$boot < 50] <- NA
for (j in 2:ncol(tax)) tax[is.na(tax[, j - 1]), j] <- NA
organelle <- (tax[, "Order"] %in% "Chloroplast") | (tax[, "Family"] %in% "Mitochondria")
nonbact <- !(tax[, "Kingdom"] %in% c("Bacteria", "Archaea"))
bact <- !(organelle | nonbact)
cat(sprintf("ASVs: %d chloroplast/mitochondria, %d non-bacterial, %d kept\n",
            sum(organelle), sum(nonbact & !organelle), sum(bact)))

getN <- function(x) sum(getUniques(x))
track <- data.frame(sample = sample, input = ft[, "reads.in"], filtered = ft[, "reads.out"],
                    denoisedF = NA, denoisedR = NA, merged = NA, nonchim = NA,
                    organelle = NA, bacterial = NA, row.names = sample)
s <- names(dadaF)
track[s, "denoisedF"] <- sapply(dadaF, getN)
track[s, "denoisedR"] <- sapply(dadaR, getN)
track[s, "merged"] <- sapply(mergers, getN)
track[s, "nonchim"] <- rowSums(nochim)[s]
track[s, "organelle"] <- rowSums(nochim[s, organelle, drop = FALSE])
track[s, "bacterial"] <- rowSums(nochim[s, bact, drop = FALSE])
write.csv(track, file.path(out, "read_tracking.csv"), row.names = FALSE)

genus_name <- ifelse(is.na(tax[bact, "Genus"]), "Unassigned", tax[bact, "Genus"])
counts <- rowsum(t(nochim[, bact, drop = FALSE]), group = genus_name)   # genera x samples
write.csv(counts, file.path(out, "genus_counts.csv"))
asv_id <- paste0("ASV", seq_len(sum(bact)))
asv <- data.frame(asv = asv_id, sequence = colnames(nochim)[bact],
                  tax[bact, , drop = FALSE], best_genus = best[bact, "Genus"],
                  genus_bootstrap = tx$boot[bact, "Genus"],
                  total_reads = colSums(nochim[, bact, drop = FALSE]), row.names = NULL)
write.csv(asv, file.path(out, "asv_taxonomy.csv"), row.names = FALSE)
asv_counts <- t(nochim[, bact, drop = FALSE]); rownames(asv_counts) <- asv_id   # ASVs x samples
write.csv(asv_counts, file.path(out, "asv_counts.csv"))
saveRDS(list(seqtab = nochim, tax = tax, track = track), file.path(work, "dada2_result.rds"))
cat("wrote genus_counts.csv (", nrow(counts), "genera x", ncol(counts), "samples ), read_tracking.csv, asv_taxonomy.csv, asv_counts.csv\n")
