#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa N-Gene Divergence
#' @date 2025-12-01
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(data.table)

# Load genetic and methylation windows to merge
geneFst <- read_csv("6.fst/geneFst.csv") |>
  filter(!is.na(aveFst))

genekw <- read_csv(file = "6.dma/kw_detectedGenes.csv") |>
  filter(!is.na(p))

# Filter out NAs and merge with SNP set
outcombo <- genekw |>
  merge(geneFst, all = TRUE) |>
  # Replace NAs with 0s
  mutate(p = replace_na(p, 1),
         aveFst = replace_na(aveFst, 0),
         aveFstp = replace_na(aveFstp, 0),
         nSNP = replace_na(nSNP, 0),
         denSNP = replace_na(denSNP, 0)) |>
  mutate(aveFst = ifelse(aveFst < 0, 0, aveFst),
         aveFstp = ifelse(aveFst < 0, 0, aveFstp)) |>
  # Assign significance categories
  mutate(msig = -log10(p) > 1.3,
         gsig = aveFst > 3*sd(aveFst),
         cat = as.factor(paste(gsig, msig, sep = " ")))

# Plot distributions
ggplot(data = outcombo,
       mapping = aes(x = aveFst)) +
  geom_histogram() # Is this type of Fst supposed to go negative?

ggplot(data = outcombo,
       mapping = aes(x = p)) +
  geom_histogram()

# Plot p vs. aveFst
Fvol <- outcombo |>
  ggplot(mapping = aes(x = aveFst,
                       y = -log10(p),
                       color = cat)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 3*sd(outcombo$aveFst), linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  theme_classic() +
  labs(title = "Gene Body Divergence",
       x = "Genetic Divergence (Average SNP Fst)",
       y = "Methylation Density Divergence (-log10(p))") +
  guides(size = "none", color = "none") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_color_manual(values = c("grey", "#F78C45", "#704d99"))

jpeg("7.transcripts/geneDivergence.png", width = 700, height = 600, quality = 100)
Fvol
dev.off()


## SNP density against KW p
# Is there a relationship between SNP density and Fst
ggplot(data = outcombo,
       mapping = aes(x = aveFst, y = denSNP)) +
  geom_point()

# Kendall Rank Correlation Test
cor.test(x = outcombo$denSNP, y = outcombo$p, method = "kendall")

# Plot SNP density against kw
densCor <- ggplot(data = outcombo,
       mapping = aes(x = denSNP, y = -log10(p))) +
  geom_point() +
  theme_classic(base_size = 16) +
  labs(title = "Gene Body Methylation and SNP Correlation",
       x = "SNP Density",
       y = "Methylation Density Divergence (-log10(p))") +
  guides(size = "none", color = "none") +
  annotate("text",
           label = "tau = 0.07, p = 2.2e-16", 
           x = 0.05, y = 2, size = 8) +
  theme(plot.title = element_text(hjust = 0.5))

jpeg("7.transcripts/densityCorrelation.jpg", width = 800, height = 600, quality = 100)
densCor
dev.off()


