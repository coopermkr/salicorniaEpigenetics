#'''''''''''''''''''''''''''''''''''''''''''
#' Genetic/Methylation Divergence Cross Comparison
#' @date 2025-03-27
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(data.table)

# Load genetic and methylation windows to merge
meanFst <- read_csv("6.fst/meanfst.csv") |>
  mutate(type = "snp")

kw <- read_csv(file = "6.dma/kw.csv") |>
  mutate(type = "smp")

# Filter out NAs and merge with SNP set
outcombo <- kw |>
  filter(!is.na(p)) |>
  mutate(significance = -log10(p) > 1.3) |>
  merge(meanFst)

# Plot p vs. Fst
vol <- ggplot(data = outcombo,
       mapping = aes(x = ZFst_mean,
                     y = -log10(p))) +
  geom_point() +
  geom_vline(xintercept = 2.5*sd(outcombo$ZFst_mean), linetype = "dashed") +
  geom_vline(xintercept = -2.5*sd(outcombo$ZFst_mean), linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  theme_classic() +
  labs(title = "Window-by-window Divergence",
       x = "Genetic Divergence (Z-transformed Mean Fst)",
       y = "Divergence in Methylation Density (-log10(p))") +
  guides (size = "none") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

png("6.dma/volcano.png", width = 800, height = 600)
vol
dev.off()


#### Pull out significant transcripts
methSig <- kw |>
  arrange(chrom) |>
  filter(p < 0.05) |>
  mutate(type = "DMR") |>
  select(chrom, window, type)

lowdiv <- meanFst |>
  filter(ZFst_mean < -2.5*sd(meanFst$ZFst_mean)) |>
  mutate(type = "Low Genetic Divergence") |>
  select(chrom, window, type)

hidiv <- meanFst |>
  filter(ZFst_mean > 2.5*sd(meanFst$ZFst_mean)) |>
  mutate(type = "High Genetic Divergence") |>
  select(chrom, window, type)

# Join the significant types together
sigTrans <- rbind(methSig, lowdiv, hidiv)

# rejoin with annotation bed file
genes <- fread("7.transcripts/gene.features.bed") |>
  filter(stringr::str_starts(V1, "^Sdep")) |>
  rename(chrom = V1,
         start = V2,
         end = V3,
         name = V4,
         score = V5,
         strand = V6,
         thicks = V7,
         thicke = V8,
         rgb = V9,
         bcount = V10,
         bsizes = V11,
         bstart = V12)

siganno <- function(chr, win, type) {
  winstart <- win - 10000
  winend <- win + 20000
  
  genes |>
    filter(str_starts(chrom, chr),
           start > win,
           start < winend) |>
    mutate(divType = type)
}


sigTranscripts <- pmap(.l = list(sigTrans$chrom, sigTrans$window, sigTrans$type), .f = siganno) |> list_rbind()

# Write out list for supp mat table
write_tsv(sigTranscripts, file = "7.transcripts/significantTranscripts.tsv", col_names = FALSE)

# create list of background genes for GO analysis
geneWindows <- pmap(.l = list(meanFst$chrom, meanFst$window, meanFst$type), .f = siganno) |> list_rbind()
methWindows <- pmap(.l = list(kw$chrom, kw$window, kw$type), .f = siganno) |> list_rbind()

# Join the lists and filter out the significant genes and duplicates
back <- rbind(geneWindows, methWindows) |>
  # Mark duplicated genes-- these appear in both SNP and SMP datasets
  mutate(dup = duplicated(name)) |>
  # distinct() only keeps the first row, so we need to arrange to make sure dup will show TRUE
  arrange(chrom, start, desc(dup)) |>
  # Filter out significant genes
  group_by(chrom) |>
  filter(!start %in% sigTranscripts$start) |>
  # Filter out one of each duplicate
  distinct(start, .keep_all = TRUE) |>
  ungroup()


write_tsv(back, file = "8.go/background.txt")

# remove duplicate hits by start position- sometimes the same gene can be recorded multiple times based on annotation rounds
goi <- sigTranscripts |>
  group_by(chrom) |>
  distinct(start, .keep_all = TRUE) |>
  ungroup() |>
  select(name) |>
  rename(gene = name)

write_tsv(goi, file = "8.go/goi.txt")
