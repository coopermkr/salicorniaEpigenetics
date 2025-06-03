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
meanFst <- read_csv("6.fst/meanfst.csv")

kw <- read_csv(file = "6.dma/kw.csv")

joined <- full_join(kw, meanFst) |>
  mutate(p = replace_na(p, 0),
         ZFst_mean = replace_na(ZFst_mean, 0),
         sites = replace_na(sites, 0))

# Calculate divergence and conservation cutoff
sd(joined$ZFst_mean) # This is 0.02, so use that for the sd number

# Write function to randomly sample each dataset, combine, and count the categories.
divBoot <- function(n, i) {
  fstSub <- sample_n(joined, size = n, replace = FALSE) |>
    select(ZFst_mean)

  kwSub <- sample_n(joined, size = n, replace = FALSE) |>
    select(p)
  
  cbind(fstSub, kwSub) |>
    mutate(methSig = p < 0.05,
           genDiv = ZFst_mean > 2.5*sd(meanFst$ZFst_mean),
           genCon = ZFst_mean < -2.5*sd(meanFst$ZFst_mean),
           boot = i) |>
    group_by(methSig, genDiv, genCon, boot) |>
    summarize(counts = n())
}

# Set number of subsampled rows to 100 and repeat the bootstrap 1000 times
bootOut <- map2_df(.x = 3500, .y = 1:1000, .f = divBoot)

# Resummarize the results across the bootstraps
bootOut |>
  group_by(methSig, genDiv, genCon) |>
  summarize(m = mean(counts),
            sdev = sd(counts))

categories <- bootOut |>
  mutate(cat = paste(methSig, genDiv, genCon)) |>
  mutate(cat = str_replace_all(cat, c("FALSE FALSE FALSE" = "No Divergence",  
                    "FALSE FALSE TRUE" = "Genetic Conservation", 
                    "FALSE TRUE FALSE" = "Genetic Divergence",
                    "TRUE FALSE FALSE" = "Methylation Divergence", 
                    "TRUE FALSE TRUE" = "Methylation Divergence/Genetic Conservation",
                    "TRUE TRUE FALSE" = "Methylation Divergence/Genetic Divergence")))

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
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "No Divergence", tag = "(A)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Genetic Divergence
genDiv <- categories |>
  filter(cat == "Genetic Divergence") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Genetic Divergence", tag = "(B)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Genetic Conservation
genCon <- categories |>
  filter(cat == "Genetic Conservation") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Genetic Conservation", tag = "(C)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Methylation Divergence
methDiv <- categories |>
  filter(cat == "Methylation Divergence") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Methylation Divergence", tag = "(D)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Methylation Divergence/Genetic Conservation
mDgC <- categories |>
  filter(cat == "Methylation Divergence/Genetic Conservation") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Methylation Divergence/Genetic Conservation", tag = "(E)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Methylation Divergence/Genetic Divergence
mDgD <- categories |>
  filter(cat == "Methylation Divergence/Genetic Divergence") |>
  ggplot(mapping = aes(x = counts)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Methylation Divergence/Genetic Divergence", tag = "(F)") +
  labs(x = "Number of Windows",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Make a combo plot
grid.arrange(noDiv, genDiv, genCon, methDiv, mDgC, mDgD, ncol = 3, 
             top=textGrob("Category Distributions",
                          gp=gpar(fontsize = 20),
                          just = "centre"))

# So this plot isn't very informative, but what might be more interesting is the
# distribution of the p values and Fst values themselves without assigning them
# to categories
distribLong <- joined |>
  select(p, ZFst_mean) |>
  pivot_longer(names_to = "type",
               values_to = "Significance",
               cols = c(p, ZFst_mean))

# Filter out the zeros since they're so over-represented
distribLong |>
  filter(test != 0) |>
  # make a histogram
  ggplot(mapping = aes(x = test)) +
  geom_histogram() +
  facet_wrap(vars(type)) +
  theme_classic()

# Or create two plots and join them together
# Methylation divergence histogram
methyl <- distribLong |>
  filter(type == "p",
         test != 0) |>
  # Make methylation plot
  ggplot(mapping = aes(x = test)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Methylation Divergence", tag = "(A)") +
  labs(x = "KW P-Value",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Genetic divergence histogram
gen <- distribLong |>
  filter(type == "ZFst_mean",
         test != 0) |>
  ggplot(mapping = aes(x = test)) +
  geom_histogram() +
  theme_classic() +
  labs(title = "Genetic Divergence", tag = "(B)") +
  labs(x = "Mean of Z-Fst",
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# Combine the two plots at subplots
grid.arrange(methyl, gen, ncol = 2, 
             top=textGrob("Distribution of Divergence Tests",
                          gp=gpar(fontsize = 20),
                          just = "centre"))
