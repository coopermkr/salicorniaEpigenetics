#'''''''''''''''''''''''''''''''''''''''''''
#' GO Analysis
#' @date 2025-04-15
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''


# Load Libraries
library(tidyverse)

# Load in background and significant gene lists
back <- read_tsv("8.go/background.ps.csv", col_names = FALSE) |>
  rename(gene = V1,
         prot = V2)

goi <- read_tsv("8.go/goi.ps.csv", col_names = FALSE) |>
  rename(gene = V1,
         prot = V2)

# Make summary tables for goi and background
back |>
  group_by()