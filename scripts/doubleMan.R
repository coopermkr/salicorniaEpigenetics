#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa Double Manhattan Plot
#' @date 2025-12-22
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
library(tidyverse)

fst <- read_csv("6.fst/transFst.csv") |>
  filter(!is.na(aveFst)) |>
  mutate(aveFst = case_when(aveFst < 0 ~ 0, TRUE ~ aveFst)) |>
  arrange(desc(aveFst))

kw <- read_csv("6.dma/kw_detectedTrans.csv") |>
  arrange(p)

combo <- merge(fst, kw, all.y = TRUE) |> # All rows should have a kw reading, but now all kw rows will have SNPs
  mutate(p = replace_na(p, 1)) |>
  filter(!is.na(aveFst)) |>
  mutate(msig = p < 0.05,
         gsig = aveFst > mean(fst$aveFst) + 3*sd(fst$aveFst),
         cat = paste(msig, gsig)) |>
  select(!aveFstp)

write_tsv(combo, "7.transcripts/divergenceMetrics.tsv")

# Double distribution plot
divPlot <- combo |> #filter(p < 1, aveFst > 0) |>
  ggplot(combo, mapping = aes(x = aveFst, y = -log10(p), color = cat)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_vline(xintercept = mean(fst$aveFst) + 3*sd(fst$aveFst), linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  theme_classic() +
  labs(title = "Gene Body Divergence",
       x = "Genetic Divergence (Mean Fst)",
       y = "Divergence in Methylation Density (-log10(p))") +
  guides(size = "none", color = "none") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_color_manual(values = c("grey", "#F78C45", "#704d99"))

jpeg("7.transcripts/transcriptDivergence.jpg", width = 800, height = 800, quality = 100)
divPlot
dev.off()


## Plot Fst
# Calculate chromosome offsets for plotting
sizes <- read_delim(file = "data/sizes.tetra.scaff.18.txt", 
                    delim = " ", col_names = c("chrom", "len")) |>
  filter(str_starts(chrom, "Sdep")) |>
  arrange(chrom) |>
  mutate(offset = cumsum(lag(len, default = 0)))

man <- combo |>
  merge(sizes) |>
  mutate(xcoord = feat + offset)

# Double manhattan
plotMan <- ggplot(data = man, mapping = aes(x = xcoord,
                                            y = aveFst,
                                            color = as.factor(chrom))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_point(mapping = aes(x = xcoord,
                           y = p,
                           color = as.factor(chrom)), size = 3, shape = 3, alpha = 0.5) +
  theme_classic(base_size = 15) +
  labs(title = "Divergence",
       x = "Transcript Start Position") +
  scale_x_continuous(label = unique(man$chrom),
                     breaks = unique(man$offset)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 15, vjust = 0.5)) +
  scale_colour_manual(values = rep(c("#704d99", "#D5CAE4"), 9)) +
  scale_y_continuous(name = "Genetic Divergence (Average Pairwise Fst)",
                     sec.axis = sec_axis(~ 1-., name = "Methylation Divergence (-p)"))



