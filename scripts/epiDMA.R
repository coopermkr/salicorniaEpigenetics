#''''''''''''''''''''''''''''''''''''''''''''''''
#' Methylation Divergence Analysis
#' @date 2025-03-05
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

library(tidyverse)
library(data.table)

winsum <- read_tsv("5.dma/winsum.tsv")
winpop <- read_tsv("5.dma/winpop.tsv")

#### Model population effects on gene windows
winpop |>
  group_by(npop) |>
  summarize(n = n())

# Note, you only return windows that you have data for, otherwise the filter returns
# a zero row frame for summarize to work with. We're not interested in those anyway

# Merge together the ndistinct population info into the winsum table
winper <- winpop |> 
  select(!pop) |>
  merge(winsum) |>
  
  # filter for windows that have info from all four pops
  filter(npop > 3) |>
  mutate(pop = substr(id, 1, 2)) |>
  
  # Calculate methylation percentage across the window for each sample
  mutate(methper = winmeth/wincov,
         # and transform it so we don't end up with 0s or 1s
         methtrans = (methper*(wincov-1)+0.5)/wincov)

#### Test for the effect of population on methylation
methdens <- read_csv(file = "6.dma/methdens.csv")

methdens <- methdens |>
  mutate(grp = paste(chrom, window, sep = "_"))

kwtest <- function(grp) {
  tmp <- kw |>
    filter(chrom == strsplit(grp, "_")[[1]][1],
           window == strsplit(grp, "_")[[1]][2]) %>%
    kruskal.test(methtrans ~ pop, data = .)
  
  data.frame(chrom = strsplit(grp, "_")[[1]][1],
             window = strsplit(grp, "_")[[1]][2],
             p = tmp$p.value)
}

kw <- map(.x = unique(methdens$grp), .f = kwtest) |> list_rbind()

kw <- read_csv(file = "6.dma/kw.csv") |>
  arrange(chrom)

# Filter out NAs
kw <- kw |>
  filter(!is.na(p)) |>
  mutate(significance = -log10(p) > 1.3)

# Sum number of significant windows
kw |>
  #group_by(chrom) |>
  summarize(nsig = sum(significance == TRUE),
            total = n(),
            prop = nsig / total)
# 0.1% of windows are significant- use as a cutoff for Fst

# Calculate chromosome offsets for plotting
sizes <- read_delim(file = "data/sizes.tetra.scaff.18.txt", 
                    delim = " ", col_names = c("chrom", "len")) |>
  filter(str_starts(chrom, "Sdep")) |>
  arrange(chrom) |>
  mutate(offset = cumsum(lag(len, default = 0)))

man <- kw |>
  merge(sizes) |>
  mutate(xcoord = window + offset,
         significance = p < 0.05)

# Plot windows by position and pvalue to make a manhattan plot
manhattan <- ggplot(data = man,
                    mapping = aes(x = xcoord,
                                  y = -log10(p),
                                  color = chrom,
                                  shape = significance)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  theme_classic(base_size = 15) +
  labs(title = "Methylation Divergence by Population",
       x = "10kb Methylation Window Position",
       y = "Kruskal Wallis of Methylation Density (-log10(p))") +
  scale_x_continuous(label = unique(man$chrom),
                     breaks = unique(man$offset)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 15, vjust = 0.5)) +
  scale_colour_manual(values = rep(c("steelblue", "lightblue"), 9))



png("6.dma/manhattanMeth.png", width = 600, height = 400)
manhattan
dev.off()

# Pull out significant transcripts
sig <- kw |>
  filter(significance == TRUE)

# Calculate windows around genes
genes <- fread("data/filtered.long.bed") |>
  filter(stringr::str_starts(V1, "^group")) |>
  rename(chrom = V1,
         start = V2,
         end = V3,
         name = V4,
         score = V5,
         strand = V6,
         thicks = V7,
         thicke = V8,
         rgb = V9,
         bcount = V10,
         bsizes = V11,
         bstart = V12)

sigtrans <- function(chr, win) {
  winend <- win + 10000
  
  genes |>
    filter(str_starts(chrom, chr),
           start > win,
           start < winend)
}

transcripts <- map2(.x = sig$chrom, .y = sig$window, .f = sigtrans) |> list_rbind()
write.csv(transcripts, file = "7.transcripts/significantTranscripts.csv")
