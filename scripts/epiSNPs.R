#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Epigenetic Methylation Calculations
#' @date 2024-11-11
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("NanoMethViz")

# Load required libraries
library(NanoMethViz)
library(dplyr)
library(ggplot2)

# This script picks up from work done on the remote server to calculate log
# methylation ratios. Here, we will subsample the dataset down to ~100k sites
# and then plot the PCA.


## Subsample dataset
lmr <- read.delim(file = "3.methViz/lmrSubset.tsv", header = TRUE,
                  sep = " ")

# Subsample lmr from 4.1 million to 100k SNPs
lmrSub <- lmr[seq(1,nrow(lmr), 42),]
# Now we have 98,741 methylation sites

#write.csv(lmrSub, file = "3.methViz/lmrResub.csv")


## PCA
lmrSub <- read.csv(file = "3.methViz/lmrResub.csv") |>
  select(!X) |>
  na.omit()

# Remove the rows that are all 0 because we removed et1
lmrFilt <- lmrSub[rowSums(lmrSub[]) != 0,]
# Now we have 93,714 sites (very similar to our genetic SNPs)

popList <- read.delim(file = "2.vcf/popmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

#' NOTE: Using the filtered dataset results in no differentiation between populations
#' I believe this is because there is so much more noise within populations
#' due to how methylation patterns are inherited, which obscures higher level 
#' structure. To overcome this noise, we need more data points, so I went with the
#' unfiltered dataset (which has still been filtered for 3x coverage).


# Plot
plot_mds(lmr, labels = NULL, groups = popList$pop) +
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Methylation MDS",
       x = "Log FC Dimension 1 (19%)",
       y = "Log FC Dimension 2 (4%)") +
  guides (size = "none") +
  theme_classic(base_size = 16) +
  stat_ellipse(type = "t", linetype = 2, mapping = aes(color = popList$pop)) +
  geom_point()


#### Diversity Calculations
lmr <- read.delim(file = "3.methViz/lmrSubset.tsv", header = TRUE,
                  sep = " ")

# Code haplotypes
lmr[lmr > 0] <- 1
lmr[lmr < 0] <- 0

# Filter by sample coverage
lmrFilt <- lmr |>
  mutate(samps = rowSums(across(where(is.numeric)))) |>
  filter(samps > 12) |>
  select(!samps)

# Pivot and add population metadata column and locus ID
lmrFiltTrans <- t(lmrFilt) |>
  cbind(popList) |>
  select(!id) |>
  relocate(pop)

write.csv(lmrFiltTrans, file = "lmrFiltTrans.csv")

# Group by locus ID and population

# Mutate/summarize mean and standard deviation per locus per population

# model sd ~ pop

# emmeans for differences between populations

# Plot x = pop, y = sd for each locus



