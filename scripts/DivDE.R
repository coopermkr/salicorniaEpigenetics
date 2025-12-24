#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa DE and Divergent Genes
#' @date 2025-12-01
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)

# Load genetic, methylation, and RNA gene significance
geneFst <- read_csv("6.fst/geneFst.csv") |>
  filter(!is.na(aveFst))

nsigs <- geneFst |> filter(aveFst > 3*sd(aveFst)) |> nrow()

genekw <- read_csv(file = "6.dma/kw_detectedTrans.csv") |>
  filter(!is.na(p))

geneDE <- read_tsv(file = "../atlas/3.edger/shootTreatment.tsv") |>
  rowwise() |>
  # Only take start and end of gene body
  mutate(feat = min(as.numeric(str_split(Start, ";", simplify = TRUE))),
         end = max(as.numeric(str_split(End, ";", simplify = TRUE)))) |>
  # Only take the chromosome
  mutate(chrom = as.character(str_split_i(Chr, ";", 1))) |>
  # Filter to only the core chromosomes
  filter(str_starts(Chr, "Sdep"),
         !is.na(FDR)) |>
  select(chrom, feat, end, Geneid, FDR, ammonium, nitrate, urea,
         logCPM, logFC.tissueV.treatmentA) |>
  ungroup()

geneAmm <- read_tsv(file = "../atlas/3.edger/shootAmmonium.tsv") |>
  rowwise() |>
  # Only take start and end of gene body
  mutate(feat = min(as.numeric(str_split(Start, ";", simplify = TRUE))),
         end = max(as.numeric(str_split(End, ";", simplify = TRUE)))) |>
  # Only take the chromosome
  mutate(chrom = as.character(str_split_i(Chr, ";", 1))) |>
  # Filter to only the core chromosomes
  filter(str_starts(Chr, "Sdep"),
         !is.na(FDR)) |>
  select(chrom, feat, end, Geneid, FDR) |>
  ungroup()

# Take the top genes in each category
topFst <- geneFst |> filter(aveFst > 3*sd(aveFst))

topDE <- geneDE |> filter(urea != FALSE,
                          nitrate != FALSE,
                          ammonium != FALSE, FDR < 0.05) |>
  slice_min(n = 33, order_by = FDR)

topKw <- genekw |> slice_min(n = 32, order_by = p)


# Restore metadata to Fst and KW genes
detectedGenes <- read_csv("../atlas/3.edger/detectedGenes.csv") |>
  rename(chrom = Chr, feat = Start, end = End) |>
  select(chrom, feat, end)

anno <- read_tsv("7.transcripts/gene.features.bed", skip = 1, col_names = FALSE) |>
  rename(chrom= X1, feat = X2, end = X3, name = X4)

# Write out files with end digit
topFst |>
  merge(anno)
  filter(type == "gene",
         str_detect(name, "ID=gene-")) |>
  select(chrom, feat, end) |>
  write_tsv("7.transcripts/topFst.tsv", col_names = FALSE)

topKw |>
  merge(anno)
  filter(type == "gene") |>
  select(chrom, feat, end)
  write_tsv("7.transcripts/topKw.tsv", col_names = FALSE)

topDE |>
  select(chrom, feat, end) |>
  write_tsv(file = "7.transcripts/topDE.tsv", col_names = FALSE)

# Combine the three sets
sigGenes <- topDE |>
  merge(topKw, all = TRUE, by = c("chrom", "feat")) |>
  merge(topFst, all = TRUE, by = c("chrom", "feat"))

#### Independence Test
# Do the number of genes we see in each genetic/methylation category match what we would expect by random chance?
# Null Hypothesis: There independence between methylation and genetic divergence
sigFst <- geneFst |> filter(aveFst > 3*sd(aveFst)) |> nrow()

sigKw <- genekw |> filter(p < 0.05) |>nrow()

undiv <- nrow(geneFst) + nrow(genekw) - sigFst - sigKw

ftDiv <- matrix(c(0, sigFst, sigKw, undiv),
                nrow = 2,
                dimnames = list(FstSig = c("yes", "no"),
                                KWSig = c("yes", "no"))) |>
  fisher.test() # one sided bc can't be lower than 0
# odds ratio = 0, p-value = 1
