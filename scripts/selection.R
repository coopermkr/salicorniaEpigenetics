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

#### Fst using specified populations
# Load in data from python scripts
winFst <- read_csv("5.fst/snp.fst.csv") |>
  # Select Fst columns and major chromosomes
  select(scaffold, start, sites, Fst_ef_es, Fst_ef_et, Fst_ef_ew, Fst_es_et, Fst_es_ew, Fst_et_ew) |>
  filter(str_starts(scaffold, "^group")) |>
  # Get rid of all NAs so we only have windows with data for all four pairwise comparisons
  na.omit() |>
  rowwise() |>
  # Average the four Fsts as a method of detecting outliers
  mutate(Fst_mean = mean(c_across(starts_with("Fst")))) |>
  arrange(desc(Fst_mean))

# plot mean Fst values
ggplot(data = winFst,
       mapping = aes(x = Fst_mean)) +
  geom_histogram()

# Pull out top 0.1% of Fst value windows
top <- nrow(winFst)*0.001
outFst <- winFst[1:top,]

## Plot Fst
# Calculate chromosome offsets for plotting
sizes <- read_delim(file = "data/sizes.tetra.scaff.18.txt", 
                    delim = " ", col_names = c("chrom", "len")) |>
  filter(str_starts(chrom, "group")) |>
  mutate(offset = cumsum(lag(len, default = 0))) |>
  rename(scaffold = chrom)

manFst <- winFst |>
  merge(sizes) |>
  mutate(xcoord = start + offset)

# Plot
plotFst <- ggplot(data = manFst,
                 mapping = aes(x = xcoord,
                                  y = Fst_mean,
                                  color = scaffold)) +
  geom_point(size = 3) +
  theme_classic(base_size = 15) +
  labs(title = "Mean SNP-wise Fst",
       x = "10kb Window Position",
       y = "Mean of Pairwise Fst") +
  scale_x_continuous(label = unique(manFst$scaffold),
                     breaks = unique(manFst$offset)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 15, vjust = 0.5))


png("manhattanFst.png", width = 800, height = 600)
plotFst
dev.off()


