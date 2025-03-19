#''''''''''''''''''''''''''''''''''''''''''''''''
#' Epigenetic Bedfile Processing
#' @date 2025-03-05
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

library(tidyverse)
library(data.table)

#'bed1 <- fread(file = "2.modkit/ef1.bam.bed", header = FALSE)

colnames(bed1) <- c("chrom", "start", "end", "code", "score", "strand", "spos", "epos",
                    "color", "cov", "fracMod", "Nmod", "Ncan", "Noth", "Ndel", "Nfail", "Ndiff", "Nno")
#### Establish dataframes
bedm <- bed1 |>
  filter(code == "m",
         cov > 1) |>
  
  # Add columns for sample and population name
  mutate(id = "ef1")

bedh <- bed1 |>
  filter(code == "h",
         cov > 1) |>
  
  # Add columns for sample and population name
  mutate(id = "ef1")

#### Bed Processing
bedProc <- function(samp){
  print(samp)
  
  # To Read in the bed file and associate it with column names
  bed <- fread(file = paste("2.modkit/", samp, ".bam.bed", sep = ""), header = FALSE)
  
  colnames(bed) <- c("chrom", "start", "end", "code", "score", "strand", "spos", "epos",
                     "color", "cov", "fracMod", "Nmod", "Ncan", "Noth", "Ndel", "Nfail", "Ndiff", "Nno")
  
  # Filter by coverage and code
  bedFilth <- bed |>
    filter(code == "h",
           cov > 1) |>
    
    # Add column for sample name
    mutate(id = samp)
  
  bedh <<- rbind(bedFilth, bedh)
  # Save to an aggregate table
  
  # Filter by coverage and code
  bedFiltm <- bed |>
    filter(code == "m",
           cov > 1) |>
    
    # Add column for sample name
    mutate(id = samp)
  
  bedm <<- rbind(bedFiltm, bedm)
  # Save to an aggregate table
}

# Iterate
sampleList <- c("ef10", "ef11", "ef12", "ef13", "ef15", "ef3", "ef4", "ef5", "ef6", "ef8", "ef9",
                "es1", "es10", "es11", "es12", "es13", "es15", "es2", "es3", "es4", "es5", "es7", "es9",
                "et1", "et10", "et11", "et13", "et14", "et15", "et2", "et3", "et4", "et7", "et8", "et9",
                "ew2", "ew10", "ew11", "ew12", "ew13", "ew2", "ew3", "ew4", "ew6", "ew7", "ew8", "ew9")

lapply(sampleList, bedProc)

# Save outputs as m, h, and combined bed files
write_tsv(bedm, file = "2.modkit/filtered.m.bed")
write_tsv(bedh, file = "2.modkit/filtered.h.bed")
rbind(bedh, bedm) |>
  write_tsv(file = "2.modkit/filtered.all.bed")