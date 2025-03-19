#'''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' Visualizing Multiple Files of Salicornia Methylation
#' @date 2024-07-11
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Install Libraries
#devtools::install_github("Shians/NanoMethViz")
#install.packages("scico") #This package may not install from the github

# Load Libraries
library(NanoMethViz)
library(dplyr)
#setwd("/scratch/workspace/c_kimballrhines001_umb_edu-ddrad/1.demux/samples")
setwd("/scratch/workspace/c_kimballrhines001_umb_edu-ddrad/3.methViz")

#### Load all data into one modbam
# Set your working directory to where you have all of your .bam and .bai files
# Create custom files with name and group using cut in bash and load them in here
#pathList <- list.files(pattern = ".sorted$")
sampList <- read.delim(file = "sampList", header = FALSE, sep = "\n")
groupList <- read.delim(file = "groupList", header = FALSE, sep = "\n")
 
# Put the sample and group info into a metadata table that matches
# the order in the pathList
meta <- data.frame(sampList, groupList) |>
    rename(sample = V1,
           group = V1.1)

# Associate all the mod files together and link to metadata
#multiMod <- ModBamResult(
#  methy = ModBamFiles(
#    samples = meta$sample,
#    paths = pathList),
#  samples = meta
#)

# Convert the object to tabix
# Set filtering options
#options("NanoMethViz.site_filter" = 3)

#modbam_to_tabix(multiMod,
#                "../../3.methViz/multi.filtered.tsv.bgz",
#                mod_code = NanoMethViz::mod_code(multiMod))
# Note: this creates one tabix file with all the samples inside it

# Load the tabix file back in (or start here if you've already created it)

#methy <- "multi.filtered.tsv.bgz" # Store path to the .tsv.bgz file

#multiRes <- NanoMethResult(methy, meta)
#multiRes

# Check for methylation in SALTY regions
#plot_region(x = multiRes, chr = "group3",
#            start = 78449000, end = 78460000)

#SALTY6 <- plot_region(x = multiRes, chr = "group6",
#            start = 74200750, end = 74400750)

#png("salty6.png", width = 1200, height = 800)
#SALTY6
#dev.off()

#### Convert the Result file into BSseq format
#multiBss <- methy_to_bsseq(multiRes)
#multiBss

#save(multiBss, file = "bsseq.filtered.RData") # Save object into .RData file

#load("bsseq.filtered.RData")

#### Generate the log-methylation-ratio matrix
#lmr <- bsseq_to_log_methy_ratio(multiBss)

#save(lmr, file = "lmr.filtered.RData")
#lmr141

load("lmr.filtered.RData")

#### MDS Plot

# Set group names
groups <- gsub('[[:digit:]]+', '', colnames(lmr))

mds <- plot_mds(lmr, labels = NULL, groups = groups)

png("mds.png", width = 600, height = 600)
mds
dev.off()

save(mds, file = "mdsGroup.RData")

#load(file = "mds.RData")

mdsColor <- mds +
  ggtitle("MDS Plot") +
  scale_color_brewer(palette = "Set1")

png("mdsColorBig.png", width = 1200, height = 1200)
mdsColor
dev.off()

#### PCA Plot
# Generate the PCA
#pca <- plot_pca(lmr) +
#  ggtitle("Filtered Methylation PCA")

#save(pca, file = "pca.filtered.RData")

# Save PCA
# png("pcaSmall.png", width = 600, height = 600)
# pca
# dev.off()
# 
# png("pcaLarge.png", width = 1200, height = 1200)
# pca
# dev.off()
# 
# #Now with color
# pca <- pca +
#   scale_colour_brewer(palette = "Set1")
# 
# png("colorSmall.png", width = 600, height = 600)
# pca
# dev.off()
# 
# png("colorLarge.png", width = 1200, height = 1200)
# pca
# dev.off()
