#'''''''''''''''''''''''''''''''''''''''''''
#' Selection Cross Comparison
#' @date 2025-03-27
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(data.table)

# Load genetic and methylation windows to merge
meanFst <- read_csv("6.fst/meanfst.csv")

kw <- read_csv(file = "6.dma/kw.csv")

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
  labs(title = "Window-by-window selection scores",
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
  mutate(significance = "DMR") |>
  select(chrom, window, significance)

lowdiv <- meanFst |>
  filter(ZFst_mean < -2.5*sd(meanFst$ZFst_mean)) |>
  mutate(significance = "Low Genetic Divergence") |>
  select(chrom, window, significance)

hidiv <- meanFst |>
  filter(ZFst_mean > 2.5*sd(meanFst$ZFst_mean)) |>
  mutate(significance = "High Genetic Divergence") |>
  select(chrom, window, significance)

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

siganno <- function(chr, win, sig) {
  winstart <- win - 10000
  winend <- win + 20000
  
  genes |>
    filter(str_starts(chrom, chr),
           start > win,
           start < winend) |>
    mutate(sigType = sig)
}

transcripts <- pmap(.l = list(sigTrans$chrom, sigTrans$window, sigTrans$significance), .f = siganno) |> list_rbind()

# Write out list for supp mat table
write_tsv(transcripts, file = "7.transcripts/significantTranscripts.tsv", col_names = FALSE)

# create list of background genes for GO analysis
gen <- transcripts
  #filter(sigType != "DMR")

back <- data.frame(gene = genes$name) |>
  filter(!gene %in% gen$name)

write_tsv(back, file = "8.go/1background_genes.sb")

# remove duplicate hits by start position- sometimes the same gene can be recorded multiple times based on annotation rounds
goi <- gen |>
  filter(!duplicated(start)) |>
  select(name) |>
  rename(gene = name)

write_tsv(goi, file = "8.go/2genes_of_interest.sb")
