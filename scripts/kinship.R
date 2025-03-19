#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Epi/Genomic Kinship
#' @date 2024-11-16
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(pcadapt)
library(popkin)
library(distances)
library(tidyverse)
library(vegan)

# Set metadata
samples <- c("ef10","ef11","ef12","ef13","ef15","ef1","ef3","ef4","ef5","ef6","ef8","ef9","es10","es11","es12","es13","es15","es1","es2","es3","es4","es5","es7","es9","et10","et11","et13","et14","et15","et1","et2","et3","et4","et7","et8","et9","ew10","ew11","ew12","ew13","ew1","ew2","ew3","ew4","ew6","ew7","ew8","ew9")

pop <- str_sub(samples, start = 1, end = 2) %>% 
  str_replace_all(c("ef" = "Folgers Marsh",  "es" = "Savin Hill", "et" = "The Creeks", "ew" = "Waquoit Bay"))

##Genetic SNPs: Read in vcf and join metadata to the matrix
sdg.vcf <- read.pcadapt(input = "2.vcf/sdg.filt.vcf", type = "vcf", type.out = "matrix")

sdg.mat <- t(sdg.vcf)
colnames(sdg.mat) <- samples

# Calculate kinship matrix
sdg.kinship <- popkin(sdg.mat, subpops = pop, want_M = T)

# Drop the original objects to free up memory

##Methylation: Read in tsv and subsample
meth <- read.table("3.methViz/lmrSubset.tsv", header = TRUE, row.names = 1)

samplerows <- sample(1:nrow(meth), size = 50000) |>
  sort()

smmeth <- meth[samplerows,]
transMeth <- t(smmeth)

#' Subsample to prevent crashes when calculating distances across the whole data set
#' At this point it's a good idea to get rid of the original objects to make space

# Distance Calculations (subsample smaller if this crashes):
mdist <- distances(transMeth)

# Convert to matrix and associate sample metadata
dmat <- as.matrix(distance_matrix(mdist))

rownames(dmat) <- samples
colnames(dmat) <- samples


## Mantel Test
mant <- mantel(sdg.kinship$kinship, 1-dmat/max(dmat),
               method = "spearman", permutations = 1000, na.rm = TRUE)
mant
# Mantel statistic Spearman's rho based on 1000 permutations:
# r = 0.3008, p = 0.010989

# Plot
diag(sdg.kinship$kinship) <- NA
diag(dmat) <- NA

m <- as.vector(1-dmat/max(dmat, na.rm = TRUE))
k <- as.vector(sdg.kinship$kinship)
mk <- cbind(m,k)

ggplot(data = mk, mapping = aes(x = k, y = m)) +
  geom_point(alpha = 0.5) +
  labs(x = "Genetic Kinship", y = "Methylation Similarity") +
  theme_bw()

ggsave(filename = "kinship.png", device = "png", units = "in", width = 5, height = 5)
