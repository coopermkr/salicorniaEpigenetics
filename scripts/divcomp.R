#''''''''''''''''''''''''''''''''''''''''''''''''
#' Diversity Window Comparison
#' @date 2025-03-24
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

#### Setup
library(tidyverse)
library(data.table)

# Load data
snpDiv <- read_tsv("6.diversity/snpDiversity.tsv") |>
  rename(chrom = CHROM,
         window = BIN_START,
         pop = population)

smpDiv <- read_tsv("6.diversity/smpDiversity.tsv") |>
  rename(methSD = sd,
         methMean = mean)

# Join the tables and keep only windows we have info for
div <- merge(snpDiv, smpDiv, by = c("chrom", "window", "pop"))

#### Test if genetic diversity predicts methylation diversity
# Make histograms of distributions
ggplot(data = div,
       mapping = aes(x = PI)) +
  geom_histogram() # Non-normal


# Make histograms of distributions
ggplot(data = div,
       mapping = aes(x = methSD)) +
  geom_histogram() # Only normal if log transformed


cor.test(x = div$PI, y = div$methSD, method = "spearman")
#rho = 0.1442, p = 2.2e-16

# Test again without zero inflation
divfilt <- div |>
  filter(PI > 0, methSD > 0)

cor.test(x = divfilt$PI, y = divfilt$methSD, method = "spearman")
#rho = -0.025, p = 0.04

# Again with a gamma

# Plot methylation vs genetic diversity
piMeth <- div |>
  mutate(across(pop, str_replace, "ef", "Folger's Marsh"),
         across(pop, str_replace, "es", "Savin Hill Cove"),
         across(pop, str_replace, "ew", "Waquoit Bay"),
         across(pop, str_replace, "et", "The Creeks")) |>
  #filter(methSD > 0) |>
  ggplot(mapping = aes(x = PI,
                       y = methSD, color = pop)) +
    #scale_y_continuous(trans = 'log10') +
    #scale_x_continuous(trans = 'log10') +
    geom_point() +
    facet_wrap(vars(pop)) +
    theme_classic() +
    labs(title = "Window-by-window diversity comparison",
       x = "Methylation diversity (standard deviation of density over 10kb window)",
       y = "Genetic diversity (Pi over 10kb window)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

png("6.diversity/piMeth.png", width = 800, height = 600)
piMeth
dev.off()

# Again with D
div |>
  #filter(methSD > 0) |>
  ggplot(mapping = aes(x = TajimaD,
                       y = methSD, color = pop)) +
  #scale_y_continuous(trans = 'log10') +
  geom_point() +
  facet_wrap(vars(pop)) +
  theme_classic()

