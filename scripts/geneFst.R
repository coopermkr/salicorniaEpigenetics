#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Gene Fst
#' @date 2025-11-04
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(vcfR)
library(tidyverse)
library(hierfstat)

# Load in VCF
sdgvcf <- read.vcfR("4.filter/snp.flt.vcf")

popList <- read.delim(file = "2.vcf/popmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Generate genInd and hierfstat objects
sdgInd <- vcfR2genind(sdgvcf, pop = popList$pop)

sdgHier <- genind2hierfstat(sdgInd)

# Calculate SNP-wise Fst
snpfst <- basic.stats(sdgHier)

save(snpfst, file = "4.filter/snpfst.Rdata")

### All-population gene Fsts
# For window Fst: Either sum all Fsts in the window or average them
# Or do SNP density: #SNPs/gene length
locFst <- as.data.frame(snpfst$perloc) |>
  rownames_to_column() %>%
  mutate(chr = str_split_i(rowname, pattern = "_", i = 1),
         pos = str_split_i(rowname, pattern = "_", i = 2)) |>
  select(!rowname)

write_tsv(locFst, file = "4.filter/snpfst.tsv")

# Load Nitrogen differential genes
anno <- read_tsv("7.transcripts/shootTreatment.tsv", col_names = TRUE) |>
  rename("chrom" = Chr,
         "start" = Start,
         "end" = End,
         "strand" = Strand,
         "name" = Geneid) |>
  mutate(chrom = substr(chrom, 1, 6),
         start = as.numeric(str_split_i(start, ";", 1)),
         end = as.numeric(str_extract(end, "[^;]+$"))) |>
  select(chrom, start, end)

# Write function to calculate Fst over gene region
geneFst <- function(featchrom, featstart, featend) {
  len = featend-featstart
  locFst |>
    filter(chr == featchrom,
           pos >= featstart,
           pos <= featend) |>
    group_by(chr) |>
    summarize(chrom = featchrom,
              feat = featstart,
              aveFst = sum(Fst)/n(),
              aveFstp = sum(Fstp)/n(),
              nSNP = n(),
              denSNP = n()/len)
}

purrAnno <- anno |>
  rename("featchrom" = chrom,
         "featstart" = start,
         "featend" = end) |>
  as.list()

pmap_dfr(.l = testAnno, .f = geneFst)
