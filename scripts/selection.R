#'''''''''''''''''''''''''''''''''''''''''''
#' SNP Selection Scan
#' @date 2025-03-25
#' @author Cooper Kimball-Rhines- from Brook's file
#'''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(pcadapt)
library(qvalue)
library(tidyverse)

# Load data
bed <- read.pcadapt("4.filter/snp.miss.bed", type = "bed", type.out = "matrix")

# Set population info
samples <- c("ef10","ef11","ef12","ef13","ef15","ef1","ef3","ef4","ef5","ef6","ef8","ef9","es10","es11","es12","es13","es15","es1","es2","es3","es4","es5","es7","es9","et10","et11","et13","et14","et15","et1","et2","et3","et4","et7","et8","et9","ew10","ew11","ew12","ew13","ew1","ew2","ew3","ew4","ew6","ew7","ew8","ew9")

pop <- str_sub(samples, start = 1, end = 2) %>% 
  str_replace_all(c("ef" = "Folgers Marsh",  "es" = "Savin Hill", "et" = "The Creeks", "ew" = "Waquoit Bay"))

# use pcadapt to look at some stuff
gente <- pcadapt(input = bed, K = 10, ploidy = 4, LD.clumping = list(size = 200, thr = 0.1)) # default 0.05 MAF filter, thin by LD

str(gente)

plot(gente, option = "screeplot") + theme_bw()

# Selection scan
genk <- pcadapt(input = bed, K = 4, ploidy = 4, LD.clumping = list(size = 200, thr = 0.1))

plot(genk, option = "screeplot") + theme_bw()
plot(gente, option = "scores", pop = pop) + theme_bw()
plot(gente, option = "scores", pop = pop, i = 3, j = 4) + theme_bw()


plot(genk, option = "manhattan")
plot(genk, option = "qqplot") # Inflated--transforming to qvalues would help

summary(genk$pvalues)

# transform p to q values
q <-  qvalue(genk$pvalues)

sigtest <- data.frame(pos = 1:length(q$qvalues),
                      p = q$pvalues,
                      q = q$qvalues) |>
  filter(!is.na(q))

sigtest |>
  filter(q < 0.05) |>
  nrow()


