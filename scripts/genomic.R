#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Genomic VCF Calculations
#' @date 2024-11-11
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(SNPfiltR)
library(vcfR)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(hierfstat)
library(adegenet)

# Load in VCF
lpvcf <- read.vcfR("2.vcf/sdep.epi.snp.flt.vcf")

# Load in popMap
popList <- read.delim(file = "2.vcf/popmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Note, we filter for depth with bcftools, so we don't need to do that step again
# So our goal with this pipeline is to investigate and filter for missingness
sdgvcf <- missing_by_snp(vcfR = lpvcf, cutoff = 0.6)

# At a 0.6 threshold, we keep 93476 SNPs

write.vcf(sdgvcf, "2.vcf/sdg.filt.vcf")


#### PCA ####
sdgvcf <- read.vcfR("2.stacks/sdg.filt.vcf")

# Generate genInd object
sdgInd <- vcfR2genind(sdgvcf, pop = popList$pop)

# Average missing data
x.sdg <- tab(sdgInd, freq = TRUE, NA.method = "mean")

# Calculate PCA
pca.sdg <- dudi.pca(x.sdg, center = TRUE, scale = FALSE, scannf = FALSE, nf = 47)

eig.perc <- 100*pca.sdg$eig/sum(pca.sdg$eig)

# Plot
par(oma = c(2, 2, 2, 2))

s.class(pca.sdg$li, fac = pop(sdgInd),
        col=c("orange", "forestgreen", "darkblue", "purple"),
        axesel = FALSE, cstar = 0, cpoint = 3)

# NOTE: adegenet is not integrated with ggplot, so we have to edit with base R
title(outer = TRUE, main = "PCA of S. depressa Genomic SNPs")

axis(1, outer = TRUE, at = 8, labels = "Principal Component 1 (3.83%)")

axis(2, outer = TRUE, at = 15, labels = "Principal Component 2 (3.50%)")

add.scatter.eig(pca.sdg$eig[1:47],3,1,2, ratio = 0.25)

#### Heterozygosity ####

# Convert to hierfstat
sdgHier <- genind2hierfstat(sdgInd)

# Calculate heterozygosity
div <- data.frame(Hs = hierfstat::Hs(sdgHier),
                  Ho = hierfstat::Ho(sdgHier))

lpGeneDiv <- hierfstat::Hs(lpHier)

write.csv(x = div, file = "2.vcf/sdgDiv.csv")

div <- read.csv(file = "2.vcf/sdgDiv.csv") |>
  rename(Pop = X)


# Graph diversity
ggplot(data = div,
       mapping = aes(x = Hs, y = Ho, 
                     color = Pop)) +
  geom_point() +
  labs(title = "L. perennis Population Heterozygosity") +
  theme_classic()
geom_segment(aes(x = 0.216, xend = 0.,
                 y = 0.23, yend = 0.255),
             linewidth = 0.1, linetype = "dashed")

