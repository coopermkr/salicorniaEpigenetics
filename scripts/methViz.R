# Load Libraries
library(NanoMethViz)
library(dplyr)

## Load in data

# NOTE: The methy() call differs from the documentation if you supply .bam files
# Change system.file() to path as I have done below and supply a path to the .bam
# Also, I am not supplying exons at this time.

mbr <- ModBamResult(
  methy = ModBamFiles(
    samples = "es11",
    paths = "../data/es11.bam.sorted"),
  
  samples = data.frame(
    sample = "es11",
    group = 1)
)

mbr #Note empty exon information

# Check for methylation in SALTY region
plot_region(x = mbr, chr = "group3",
            start = 78450000, end = 78500000)

plot_region(x = mbr, chr = "group1",
            start = 0, end = 1500000)

#### Load in all bam files
multibr <- ModBamResult(
  methy = ModBamFiles(
    samples = "es11",
    paths = "../data/es11.bam.sorted"),
  
  samples = data.frame(
    sample = "es11",
    group = 1)
)



## Convert modbam to tabix
modbam_to_tabix(mbr, 
                "../data/es11.tsv.bgz", 
                mod_code = NanoMethViz::mod_code(mbr))

## Load in tabix data
# Again, because we aren't accessing data from the package, we just supply a path to our tabix
methy <- "../data/es11.tsv.bgz"
sample <- "es11"
group <- "es"
sample_anno <- data.frame(sample, group, stringsAsFactors = FALSE)

# Create methy object
es11 <- NanoMethResult(methy, sample_anno)
es11

## Create PCA
# Convert methy to BSseq object
bss <- methy_to_bsseq(es11)
bss

save(bss, file = "bsseq.RData") # Save object into .Rdata file
load("bsseq.RData")

# Generate log-methylation ratio
# If you're interested in certain genes,
# you can supply a regions argument with exon info
lmr <- bsseq_to_log_methy_ratio(bss)

save(lmr, file = "lmr.RData") # Save object into .Rdata file

load("lmr.RData")

# Plot PCA
plot_pca(lmr) +
  ggtitle("PCA Plot")
