#''''''''''''''''''''''''''''''''''''''''''''''''
#' 
#' Divergence Distribution Sampling
#' @date 2025-06-03
#' @author Cooper Kimball-Rhines
#' 
#''''''''''''''''''''''''''''''''''''''''''''''''

#### Setup
library(tidyverse)
library(data.table)
library(gridExtra)
library(grid)

# Load in real data
# fst <- read_csv("6.fst/transFst.csv") |>
#   filter(!is.na(aveFst)) |>
#   mutate(aveFst = case_when(aveFst < 0 ~ 0, TRUE ~ aveFst))
# 
# kw <- read_csv("6.dma/kw_detectedTrans.csv")
# 
# combo <- merge(fst, kw, all.y = TRUE) |> # All rows should have a kw reading, but now all kw rows will have SNPs
#   mutate(p = replace_na(p, 1)) |>
#   filter(!is.na(aveFst)) |>
#   mutate(msig = p < 0.05,
#          gsig = aveFst > mean(fst$aveFst) + 3*sd(fst$aveFst),
#          cat = paste(msig, gsig))

allsigs <- read_csv(file = "8.go/allGeneSigMetrics.csv")

filter(allsigs, FDR < 0.05) |> nrow()

sigCats <- allsigs |>mutate(ammonium = (FDR < 0.05 & (logFC.treatmentA > 1 | logFC.treatmentA < -1)),
                    nitrate = (FDR < 0.05 & (logFC.treatmentN > 1 | logFC.treatmentN < -1)),
                    urea = (FDR < 0.05 & (logFC.treatmentU > 1 | logFC.treatmentU < -1)),
                    any = ammonium == T | nitrate == T | urea == T,
                    all = ammonium == T & nitrate == T & urea == T)

# Make venn diagram
library(ggVennDiagram)
venn <- sigCats |>
  filter(FDR < 0.05) |>
  select(ammonium, nitrate, urea) |>
  mutate(ammonium = ammonium*row_number(),
         nitrate = nitrate*row_number(),
         urea = urea*row_number())

# make a numbered list of genes for each treatment
a <- venn |> filter(ammonium != 0) |> select(ammonium)
n <- venn |> filter(nitrate != 0) |> select(nitrate)
u <- venn |> filter(urea != 0) |> select(urea)

x <- list(Ammonium = a$ammonium, Nitrate = n$nitrate, Urea = u$urea)

# Make the Venn diagram
treatVennPlot <- ggVennDiagram(x) +
  theme_classic() +
  theme(legend.position = "none",
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

jpeg(filename = "../atlas/3.edger/vennTreat.jpg", width = 600, height = 600, quality = 100)
treatVennPlot
dev.off()


# Write function to randomly sample each dataset, combine, and count the categories.
divBoot <- function(n, i) {
  fstSub <- sample_n(outcombo, size = n, replace = FALSE) |>
    select(aveFst)

  kwSub <- sample_n(outcombo, size = n, replace = FALSE) |>
    select(methP)
  
  cbind(fstSub, kwSub) |>
    mutate(methSig = methP < 0.05,
           genDiv = aveFst > mean(outcombo$aveFst) + 3*sd(outcombo$aveFst),
           boot = i) |>
    group_by(methSig, genDiv, boot) |>
    summarize(counts = n())
}

# Set number of subsampled rows to 3500 and repeat the bootstrap 1000 times
bootOut <- map2_df(.x = 2500, .y = 1:1000, .f = divBoot)

# Resummarize the results across the bootstraps
bootOut |>
  group_by(methSig, genDiv) |>
  summarize(m = mean(counts),
            sdev = sd(counts))

categories <- bootOut |>
  mutate(cat = paste(methSig, genDiv)) |>
  mutate(cat = str_replace_all(cat, c("FALSE FALSE" = "No Divergence",  
                    "FALSE TRUE" = "Genetic Divergence",
                    "TRUE FALSE" = "Methylation Divergence", 
                    "TRUE TRUE" = "Methylation Divergence/Genetic Divergence")))

# Plot the categories
ggplot(data = categories,
       mapping = aes(x = cat,
                     y = counts,
                     color = cat)) +
  geom_boxplot() +
  theme_classic()

# Create four histograms to plot together on their own y-axes
# No divergence
noDiv <- categories |>
  filter(cat == "No Divergence") |>
  ggplot(mapping = aes(x = counts, color = cat, fill = cat)) +
  geom_histogram() +
  theme_light() +
  labs(title = "No Divergence", tag = "(A)") +
  labs(x = "Number of Transcripts",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "grey") +
  scale_fill_manual(values = "grey") +
  guides (color = "none", fill = "none")

# Genetic Divergence
genDiv <- categories |>
  filter(cat == "Genetic Divergence") |>
  ggplot(mapping = aes(x = counts, color = cat, fill = cat)) +
  geom_histogram() +
  theme_light() +
  labs(title = "Genetic Divergence", tag = "(B)") +
  labs(x = "Number of Transcripts",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "#704D99") +
  scale_fill_manual(values = "#704D99") +
  guides (color = "none", fill = "none")

# Genetic Conservation
# genCon <- categories |>
#   filter(cat == "Genetic Conservation") |>
#   ggplot(mapping = aes(x = counts)) +
#   geom_histogram() +
#   theme_classic() +
#   labs(title = "Genetic Conservation", tag = "(C)") +
#   labs(x = "Number of Transcripts",
#        y = "Frequency") +
#   theme(plot.title = element_text(hjust = 0.5))

# Methylation Divergence
methDiv <- categories |>
  filter(cat == "Methylation Divergence") |>
  ggplot(mapping = aes(x = counts, color = cat, fill = cat)) +
  geom_histogram() +
  theme_light() +
  labs(title = "Methylation Divergence", tag = "(C)") +
  labs(x = "Number of Transcripts",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "#F78C45") +
  scale_fill_manual(values = "#F78C45") +
  guides (color = "none", fill = "none")

# Methylation Divergence/Genetic Conservation
mDgC <- categories |>
  filter(cat == "Methylation Divergence/Genetic Conservation") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_light() +
  labs(title = "Methylation Divergence/Genetic Conservation", tag = "(E)") +
  labs(x = "Number of Transcripts",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Methylation Divergence/Genetic Divergence
mDgD <- categories |>
  filter(cat == "Methylation Divergence/Genetic Divergence") |>
  ggplot(mapping = aes(x = counts, color = "cat", fill = "cat")) +
  geom_histogram() +
  theme_light() +
  labs(title = "Methylation Divergence\nand Genetic Divergence", tag = "(D)") +
  labs(x = "Number of Transcripts",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "#b46d6f") +
  scale_fill_manual(values = "#b46d6f") +
  guides (color = "none", fill = "none") +
  coord_cartesian(xlim = c(0, 2))

# Make a combo plot
png("evoApps/bootCats.png", width = 800, height = 600)
grid.arrange(noDiv, genDiv, methDiv, mDgD, ncol = 2, 
             top=textGrob("Category Distributions",
                          gp=gpar(fontsize = 20),
                          just = "centre"))
dev.off()

# So this plot isn't very informative, but what might be more interesting is the
# distribution of the p values and Fst values themselves without assigning them
# to categories
distribLong <- outcombo |>
  select(methP, aveFst) |>
  pivot_longer(names_to = "type",
               values_to = "Significance",
               cols = c(methP, aveFst))

# Filter out the zeros since they're so over-represented
distribLong |>
  filter(Significance != 0) |>
  # make a histogram
  ggplot(mapping = aes(x = Significance)) +
  geom_histogram() +
  facet_wrap(vars(type)) +
  theme_classic()

# Or create two plots and join them together
# Methylation divergence histogram
methyl <- distribLong |>
  filter(type == "methP") |>
  # Make methylation plot
  ggplot(mapping = aes(x = Significance, color = type, fill = type)) +
  geom_histogram() +
  theme_light(base_size = 16) +
  labs(title = "Methylation Divergence", tag = "(B)") +
  labs(x = "KW Test P-Value",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "#F78C45") +
  scale_fill_manual(values = "#F78C45") +
  guides (color = "none", fill = "none") +
  geom_vline(xintercept = 0.05, linetype = "dashed")

# Genetic divergence histogram
gen <- distribLong |>
  filter(type == "aveFst",
         Significance != 0) |>
  ggplot(mapping = aes(x = Significance, color = type, fill = type)) +
  geom_histogram() +
  theme_light(base_size = 16) +
  labs(title = "Genetic Divergence", tag = "(A)") +
  labs(x = "Transcript Ave. Fst",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = "#704d99") +
  scale_fill_manual(values = "#704d99") +
  guides (color = "none", fill = "none") +
  geom_vline(xintercept = mean(outcombo$aveFst) + 3*sd(outcombo$aveFst),
             linetype = "dashed")

# Combine the two plots at subplots
png(filename = "evoApps/divTestDistribs.png", height = 600, width = 800)
grid.arrange(gen, methyl, ncol = 2, 
             top=textGrob("Distribution of Divergence Tests",
                          gp=gpar(fontsize = 20),
                          just = "centre"))
dev.off()
