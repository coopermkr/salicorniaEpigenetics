#'''''''''''''''''''''''''''''''''''''''''''
#' Gene Ranking Analysis
#' @date 2026-02-10
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''

library(tidyverse)

# Read in annotation
anno <- read_tsv("data/correctLabels.gff3", skip = 1, col_names = FALSE) |>
  rename(Chr = X1, Feature = X3, Start = X4, End = X5, name = X9) |>
  filter(str_detect(name, "transcript_id"),
         Feature == "gene") |>
  mutate(name = str_split_i(name, pattern = "transcript_id=", 2)) |>
  distinct()

# Read in data
transcript <- read_tsv("data/allTreatTranscripts.tsv")

diverge <- read_tsv("7.transcripts/divergenceMetrics.tsv") |>
  rename("Chr" = chrom, "Start" = feat, "End" = end)

allMetrics <- merge(transcript, diverge, by = c("Chr", "Start", "End"), all = TRUE)
# no loss adding in the end coordinate

# Using CDS sequences joining without end coords: 13973
# Using mRNA sequences joining without end coords: 11760
# Using mRNA sequences joining with end coords: 11617
# Using gene sequences joining without end coords: 11143
# Using gene sequences joining with end coords: 10447

filtMetrics <- allMetrics |>
  merge(anno, by = c("Chr", "Start", "End")) |>
  select(name, logFC.treatmentA, logFC.treatmentN, logFC.treatmentU,
         logCPM, FDR, aveFst, p) |>
  rename("methP" = p) |>
  na.omit() |>
  distinct()

# How many duplicates are there?
filtMetrics |> group_by(name) |>
  count() |> filter(n > 1)
# CDS Without joining by end coordinate there are 201 duplicates
# CDS With joining by end coordinates there are 200 duplicates
# mRNA with/without joining by end coordinates there are 200 duplicates
# gene with/without joining by end coordinates there are 0 duplicates

#write_csv(filtMetrics, "8.go/allGeneSigMetrics.csv")

# total of 10,447 genes with data for DE, DG, and DM
filtMetricsLogDE <- filtMetrics |> mutate(sigDE = FDR < 0.05 & (logFC.treatmentA > 1 | logFC.treatmentA < -1 | 
                                      logFC.treatmentN > 1 | logFC.treatmentA < -1 |
                                      logFC.treatmentU > 1 | logFC.treatmentU < -1))
# 1863 DE transcripts
filtMetrics |> filter(methP < 0.05) # 5 DM transcripts with no overlap between DE or FST
filtMetrics |> filter(aveFst > (sd(aveFst)*3+mean(aveFst))) |> nrow() # 226 DG transcripts

# Is there overlap between DG and DE?
filtMetricsLogDE |> filter(aveFst > (sd(aveFst)*3+mean(aveFst)),
                      sigDE == TRUE) # Yes! Actually lots! 41 genes overlap.

# Can we do a chi-squared to see if this is an overrepresentation?
matrix(c(13, 218, 696, 9520),
             nrow = 2,
             dimnames = list(DE = c("divExp", "backExp"),
                             DG = c("divGen", "backGen"))) |>
  fisher.test()
# DG/DE: p = 0.687 (0.435, 1.47)
# DM/DE: p = 1 (0, 15.0)
# DM+DG/DE: p = 0.596 (0.425, 1.43)

#### Lollipop Graphs
# Read in Transcriptome and wild marsh GOs
growthGO <- read_csv("evoApps/SigGOtermsFDRNitrogenSource.csv") |>
  mutate(logMFP = -log10(CorrectedMFPvalue)) |>
  slice_min(order_by = CorrectedMFPvalue, n = 6)

growthGO$`GO Term` <-  factor(growthGO$`GO Term`, levels = unique(growthGO$`GO Term`[order(growthGO$logMFP)]), ordered = TRUE)


wildGO <- read_csv("evoApps/GOtermsWildMarsh.csv") |>
  mutate(logMFP = -log10(CorrectedMFPvalue))

wildGO$`GO Term` <-  factor(wildGO$`GO Term`, levels = unique(wildGO$`GO Term`[order(wildGO$logMFP)]), ordered = TRUE)

# Make plots
growthPlot <- ggplot(data = growthGO, 
       mapping = aes(x = logMFP, y = `GO Term`, size = NumGenes)) +
  geom_count() +
  theme_light(base_size = 22) +
  theme(axis.title.y=element_blank()) +
  scale_size_continuous(limits=c(1, 165), breaks=seq(0, 165, by=20), guide = "legend",
                        range = c(1,15)) +
  labs(x = "-log10(Corrected P Value)")
  #scale_x_continuous(limits = c(5, 10))

jpeg(filename= "8.go/growthGOs.jpg", width = 10, height = 6.5, units = "in", res = 100)
growthPlot
dev.off()

addline_format <- function(x,...){
  gsub(',','\n',x)
}

wildPlot <- ggplot(data = wildGO, 
       mapping = aes(x = logMFP, y = `GO Term`, size = NumGenes, color = Type)) +
  geom_count() +
  theme_light(base_size = 22) +
  scale_color_manual(values = c("#704d99", "#F78C45")) +
  theme(axis.title.y=element_blank()) +
  scale_size_continuous(limits=c(1, 165), breaks=seq(0, 165, by=40), guide = "legend",
                        range = c(1,15)) +
  labs(x = "-log10(Corrected P Value)", size = "N Genes") +
  scale_y_discrete(labels = addline_format(wildGO$`GO Term`))

jpeg(filename= "8.go/wildGOs.jpg", width = 5, height = 3, units = "in", res = 100)
wildPlot
dev.off()