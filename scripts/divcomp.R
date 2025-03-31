#''''''''''''''''''''''''''''''''''''''''''''''''
#' Diversity Window Comparison
#' @date 2025-03-24
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

#### Setup
library(tidyverse)
library(data.table)

# Load data
snpDiv <- read_tsv("5.diversity/snpDiversity.tsv") |>
  rename(chrom = CHROM,
         window = BIN_START,
         pop = population)
smpDiv <- read_tsv("5.diversity/smpDiversity.tsv") |>
  rename(methSD = sd,
         methMean = mean)

# Join the tables and keep only windows we have info for
div <- merge(snpDiv, smpDiv)

# Plot methylation vs genetic diversity
div |>
  filter(methSD > 0) |>
  ggplot(mapping = aes(x = methSD,
                       y = PI, color = pop)) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10') +
    geom_point() +
    facet_wrap(vars(pop)) +
    theme_classic()

# Again with D
div |>
  filter(methSD > 0) |>
  ggplot(mapping = aes(x = methSD,
                       y = TajimaD, color = pop)) +
  scale_x_continuous(trans = 'log10') +
  geom_point() +
  facet_wrap(vars(pop)) +
  theme_classic()


