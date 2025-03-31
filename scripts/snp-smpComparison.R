#'''''''''''''''''''''''''''''''''''''''''''
#' Selection Cross Comparison
#' @date 2025-03-27
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)

# Load genetic and methylation windows to merge
meanFst <- read_csv("5.fst/meanfst.csv")

kw <- read_csv(file = "5.dma/kw.csv")

# Filter out NAs and merge with SNP set
outcombo <- kw |>
  filter(!is.na(p)) |>
  mutate(significance = -log10(p) > 1.3) |>
  merge(meanFst)

# Plot p vs. Fst
ggplot(data = outcombo,
       mapping = aes(x = ZFst_mean,
                     y = -log10(p))) +
  geom_point() +
  geom_vline(xintercept = 3*sd(outcombo$ZFst_mean)) +
  geom_vline(xintercept = -3*sd(outcombo$ZFst_mean)) +
  geom_hline(yintercept = 1.3) +
  theme_classic()



