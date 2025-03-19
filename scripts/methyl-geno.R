# Comparing methyl to genotype in Salicornia
# Written by Brook Moyers
# Started 25 July 2024

setwd("~/methyl-geno/")

# load packages
library(pcadapt)
library(tidyverse)

# load data
# metadata
md <- read_table("groups.txt", col_names = c("star","barcode","id")) %>% 
  select(barcode:id) %>% 
  mutate(pop = str_sub(id, start = 1, end = 2))

# vcf
bed <- read.pcadapt("/project/pi_brook_moyers_umb_edu/salicornia/epigenetics/brook_scripts/snp.flt.bed", type = "bed", type.out = "matrix") 

samples <- c("ef10","ef11","ef12","ef13","ef15","ef1","ef3","ef4","ef5","ef6","ef8","ef9","es10","es11","es12","es13","es15","es1","es2","es3","es4","es5","es7","es9","et10","et11","et13","et14","et15","et1","et2","et3","et4","et7","et8","et9","ew10","ew11","ew12","ew13","ew1","ew2","ew3","ew4","ew6","ew7","ew8","ew9")

pop <- str_sub(samples, start = 1, end = 2) %>% 
  str_replace_all(c("ef" = "Folgers Marsh",  "es" = "Savin Hill", "et" = "The Creeks", "ew" = "Waquoit Bay"))

# use pcadapt to look at some stuff
gente <- pcadapt(input = bed, K = 10, ploidy = 4, LD.clumping = list(size = 200, thr = 0.1)) # default 0.05 MAF filter, thin by LD

str(gente)

plot(gente, option = "screeplot") + theme_bw()
ggsave(filename = "geneticPCAscree.png", device = "png", units = "in", width = 4, height = 2.5)

plot(gente, option = "scores", pop = pop, i = 1, j = 2) + 
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73","#CC79A7")) +
  scale_shape_manual(values = c(15, 19, 17,18)) +
  labs(color = "Marsh") +
  theme_bw()

ggsave(filename = "geneticPCA12.png", device = "png", units = "in", width = 6, height = 4)

plot(gente, option = "scores", pop = pop, i = 1, j = 3) + 
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73","#CC79A7")) +
  labs(color = "Marsh") +
  theme_bw()

ggsave(filename = "geneticPCA13.png", device = "png", units = "in", width = 6, height = 4)

# ok we choose a K of 3
genk <- pcadapt(input = bed, K = 3, ploidy = 4, LD.clumping = list(size = 200, thr = 0.1))

plot(genk, option = "manhattan")
plot(genk, option = "qqplot") # looks a bit inflated still
summary(genk$pvalues)
which(-log10(genk$pvalues) > 15) #okay let's leave this alone for now

# okay, let's make a kinship matrix
mat<-t(bed)
colnames(mat) <- samples

library(popkin)
kinship <- popkin(mat, subpops = pop, want_M = T)
kinin <- inbr_diag(kinship$kinship)

in_coef <- inbr(kinship$kinship)

inc <- cbind(in_coef, pop, samples)
inc <- data.frame(inc)
inc$in_coef <- as.numeric(inc$in_coef)

ggplot(inc, aes(y = in_coef, x = pop, color = pop)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.25, show.legend = F) +
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73","#CC79A7")) +
  labs(x = "Marsh", y = "Individual inbreeding coefficient") +
  theme_bw()

zin <- aov(in_coef ~ pop, data = inc) 
anova(zin) # highly significant difference!
TukeyHSD(zin)

ggsave(filename = "inbreeding_by_pop.png", device = "png", units = "in", width = 5, height = 4)

png("kinship_heatmap.png", units = "in", width = 5.5, height = 4, res = 600)
plot_popkin(kinin, labs = pop, labs_even = TRUE, labs_line = 1, labs_cex = 0.7, mar = 2, ylab = NULL)
dev.off()  

fst(kinship$kinship) # 0.009929493 very low overall
pairwise_fst <- pwfst(kinship$kinship)
leg_title <- expression(paste('Pairwise ', F[ST]))

png("fst_heatmap.png", units = "in", width = 5.5, height = 4, res = 600)
plot_popkin(pairwise_fst, labs = pop, labs_even = TRUE, labs_line = 1, labs_cex = 0.7, mar = 2, ylab = NULL, leg_title = leg_title)
dev.off()  

rm(bed) #for space

# can we do the same with the methyl matrix?
meth <- read.table("/project/pi_brook_moyers_umb_edu/salicornia/epigenetics/3.methViz/lmrSubset.tsv", header = T, row.names = 1)

samplerows <- sample(1:nrow(meth), size = 50000) %>% sort(.) # random subsample of rows, any larger and this session with 48G RAM crashes
smmeth <- meth[samplerows,] # actually sample
metht <- t(smmeth) # transpose

library(distances) # apparently good with large data, well maybe not that big
mdist <- distances(metht)
md <-distance_matrix(mdist)
mdm <- as.matrix(md)
rownames(mdm) <- samples
colnames(mdm) <- samples

png("methylation_heatmap.png", units = "in", width = 5.5, height = 4, res = 600)
plot_popkin(1-mdm/max(mdm), leg_title = "Methylation Distance", labs = pop, labs_even = TRUE, labs_line = 1, labs_cex = 0.7, mar = 2, ylab = NULL)
# something odd with et15, which looks unlike everyone
dev.off()  


#mantel test
library(vegan)

man <- mantel(kinship$kinship, 1-mdm/max(mdm), method = "spearman", permutations = 1000, na.rm = T)
man #significantly correlated!

diag(kinship$kinship) <- NA
diag(mdm) <- NA

m <- as.vector(1-mdm/max(mdm, na.rm = T))
k <- as.vector(kinship$kinship)
mk <- cbind(m,k)

mk <- data.frame(mk[!duplicated(mk),])

ggplot(data = mk, aes(x = k, y = m)) +
  geom_point(alpha = 0.5) +
  #geom_smooth(method='lm', color = "tomato") +
  labs(x = "Genetic Kinship", y = "Methylation Similarity") +
  theme_bw()

ggsave(filename = "kin-methyl-corr.png", device = "png", units = "in", width = 5, height = 5)

# why is et15 so weird? maybe coverage?
meth_zero <- colSums(as.matrix(meth)==0)/nrow(meth)
methzero <- cbind(meth_zero, pop, samples)
methzero <- data.frame(methzero)
methzero$meth_zero <- as.numeric(methzero$meth_zero)
# nothing special about it

ggplot(methzero, aes(y = 1 - meth_zero, x = pop, color = pop)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = 0.25, show.legend = F) +
  scale_color_manual(values = c("#E69F00","#56B4E9","#009E73","#CC79A7")) +
  labs(x = "Marsh", y = "Proportion of site methylated") +
  theme_bw()

anova(aov(meth_zero ~ pop, data = methzero)) # no significant difference

ggsave(filename = "p_methylated.png", device = "png", units = "in", width = 5, height = 4)

# figure out the positions for Fst outliers
pos <- read_table(file = "positions.txt", col_names = "position")

fst15 <- c(89056, 89066, 89068, 89178, 89197, 89199, 89226, 159419, 262229, 313062, 437838, 468762, 508944, 533883, 603724, 646914, 663826, 842594, 1001694, 1285715, 1407144, 1431862, 1776896, 1880270, 2188998)

pos_fst <- pos[fst15,]

pos_fst <- cbind(fst15,pos_fst)

