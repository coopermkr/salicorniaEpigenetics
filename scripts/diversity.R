#''''''''''''''''''''''''''''''''''''''''''''''''
#' Epigenetic Diversity Calculations
#' @date 2025-03-24
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

#### Setup
library(tidyverse)
library(data.table)

#### Diversity Calculation
# Methylation averaged over tiles to deal with runs of methylation
# Import chromosome sizes and divide by window size to get nwindows
# Then manually create two columns for start and end of window
# Double filter the bed file and sum the methylated reads/unmethylated to get percent
# Iterate over chromosomes and windows until done

# Start with the already calculate window densities:
meth <- read_csv("5.dma/methdens.csv")

# Calculate population-level standard deviation as a measure of methylation diversity for each window
methDev <- meth |>
  filter(npop == 4) |>
  group_by(chrom, window, pop) |>
  summarize(mean = mean(methper),
            sd = sd(methper)) |>
  filter(!is.na(sd),
         # Filter out nonvariable regions
         sd > 0)

# Graph deviations
methVio <- methDev |>
  # Filter out super high standard devs
  #filter(sd < 0.3) |>
  
  # Plot sd by population
  ggplot(mapping = aes(x = pop, y = sd, color = pop)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.05) +
  geom_violin(fill = NA, linewidth = 1) +
  
  # Put it on a log scale so you can actually see the lower end
  scale_y_continuous(trans = 'log10') +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Methylation Diversity",
       x = "Population",
       y = "Standard Deviation by 1000bp Window") +
  guides (size = "none") +
  theme_classic(base_size = 16)

png("methViolin.png", width = 600, height = 800)
methVio
dev.off()

# Make a boxplot
methBox <- methDev |>
  # Filter out super high standard devs
  #filter(sd < 0.3) |>
  
  # Plot sd by population
  ggplot(mapping = aes(x = pop, y = sd, color = pop)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.05) +
  geom_boxplot(fill = NA, linewidth = 1) +
  
  # Put it on a log scale so you can actually see the lower end
  scale_y_continuous(trans = 'log10') +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Methylation Diversity",
       x = "Population",
       y = "Standard Deviation by 1000bp Window") +
  guides (size = "none") +
  theme_classic(base_size = 16)

png("methBox.png", width = 600, height = 800)
methBox
dev.off()

# Test differences between the populations
library(performance)
library(emmeans)
library(glmmTMB)
library(DHARMa)

kruskal.test(sd ~ pop, data = methDev)
# The populations are technically different

# Make a model to do emmeans on
methlm <- lm(data = methDev,
              sd ~ pop)

r2(methlm)
check_model(methlm)
# This clearly needs a glm

# Make a gamma model
methglm <- glm(data = methDev,
               sd ~ pop,
               family = Gamma)

#check_model(methglm) # Looks lovely
check_posterior_predictions(methglm)

# Contrast populations with emmeans
methEM <- emmeans(methglm, specs = ~pop) |>
  contrast(method = "pairwise") |>
  confint()

methEM |> plot() +
  geom_vline(xintercept = 0, color = 'red', lty = 2)
# Basically The Creeks are very different from all the other populations
# I tested this with and without zeros
# But I'm not convinced this difference is actually relevant