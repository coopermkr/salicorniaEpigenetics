#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Gene Fst
#' @date 2025-11-26
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(vcfR)
library(tidyverse)
library(hierfstat)
library(adegenet)

setwd("/work/pi_brook_moyers_umb_edu/salicornia/epigenetics/")

# Get chromosome from task ID in the array
task_id <- commandArgs(trailingOnly=TRUE)[1] |> as.numeric()
#loop <- paste("4.filter/Sdep", 1, "A.subset.vcf", sep = "")
loop <- paste("4.filter/Sdep", task_id, "A.recode.vcf", sep = "")
#loop <- paste("4.filter/Sdep", task_id, "B.recode.vcf", sep = "")
paste(loop)

# Load in subset VCF
sdgvcf <- read.vcfR(loop)

popList <- read.delim(file = "3.vcf/popmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Generate genInd and hierfstat objects
snpfst <- vcfR2genind(sdgvcf, pop = popList$pop) |>
  
  genind2hierfstat() |>

# Calculate SNP-wise Fst
  basic.stats()

### All-population gene Fsts
# For window Fst: Either sum all Fsts in the window or average them
# Or do SNP density: #SNPs/gene length
locFst <- as.data.frame(snpfst$perloc) |>
  rownames_to_column() %>%
  mutate(chr = str_split_i(rowname, pattern = "_", i = 1),
         pos = str_split_i(rowname, pattern = "_", i = 2)) |>
  select(!rowname)

write_csv(locFst, file = paste("6.fst/perloc/", loop, sep = ""))

# Load Nitrogen differential genes
# anno <- read_tsv("data/shootTreatment.tsv", col_names = TRUE) |>
#   rename("chrom" = Chr,
#          "start" = Start,
#          "end" = End,
#          "strand" = Strand,
#          "name" = Geneid) |>
#   mutate(chrom = substr(chrom, 1, 6),
#          start = as.numeric(str_split_i(start, ";", 1)),
#          end = as.numeric(str_extract(end, "[^;]+$"))) |>
#   select(chrom, start, end)
# 
# # Write function to calculate Fst over gene region
# geneFst <- function(featchrom, featstart, featend) {
#   len = featend-featstart
#   locFst |>
#     filter(chr == featchrom,
#            pos >= featstart,
#            pos <= featend) |>
#     group_by(chr) |>
#     summarize(chrom = featchrom,
#               feat = featstart,
#               aveFst = sum(Fst)/n(),
#               aveFstp = sum(Fstp)/n(),
#               nSNP = n(),
#               denSNP = n()/len)
# }
# 
# purrAnno <- anno |>
#   rename("featchrom" = chrom,
#          "featstart" = start,
#          "featend" = end) |>
#   as.list()
# 
# Fsts <- pmap_dfr(.l = anno, .f = geneFst)
# 
# write_csv(Fsts, file = paste("6.fst/", loop, sep = ""))
