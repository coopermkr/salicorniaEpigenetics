#''''''''''''''''''''''''''''''''''''''''''
#' Genomic Diversity
#' @date 2025-03-24
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''

# Tajima's D and Pi calculated at the population level using vcftools on server
# Output files pulled directly from server to process and graph here

# Load libraries
library(tidyverse)

#### Process output files
divProc <- function(pop){
  print(pop)
  dpath <- paste("6.diversity/", pop, ".pop.div.Tajima.D", sep = "")
  pipath <- paste("6.diversity/", pop, ".pop.div.windowed.pi", sep = "")
  
  d <- read_tsv(dpath) |>
    mutate(BIN_START = BIN_START + 1)
  pi <- read_tsv(pipath) |>
    rename(N_SNPS = N_VARIANTS)
  
  merge(d, pi) |>
    mutate(population = pop)
}

div <- map(.x = c("es", "ef", "et", "ew"), .f = divProc) |> list_rbind()

# Change NAs to 0
div <- div |>
  mutate(TajimaD = replace_na(TajimaD, 0),
         PI = replace_na(PI, 0))

# Sum number of populations represented in each window
divPops <- div |>
  group_by(CHROM, BIN_START) |>
  summarize(npops = n_distinct(population)) |>
  filter(npops == 4)

# Merge back and only keep windows with all four pops
divfilt <- merge(divPops, div, by.x = c("CHROM", "BIN_START")) |>
  arrange(CHROM, BIN_START, population)

# Make sure every window has 4 populations
unique(divfilt$npops)

write_tsv(divfilt, file = "6.diversity/snpDiversity.tsv")
divfilt <- read_tsv("6.diversity/snpDiversity.tsv")

##### Test differences between populations
library(emmeans)
library(performance)
library(glmmTMB)
library(DHARMa)
library(broom)

# How many windows?
divfilt |>
  select(CHROM, BIN_START, PI, population) |>
  pivot_wider(names_from = population,
              values_from = PI) |>
  nrow()

# Chart distribution
divfilt |>
  filter(PI < 0.0001) |>
  ggplot(mapping = aes(x = PI)) +
  geom_histogram() +
  facet_wrap(vars(population))

# Summarize PI Diversity
divfilt |>
  group_by(population) |>
  summarize(m = mean(PI),
            max = max(PI),
            min = min(PI),
            sd = sd(PI),
            dens = mean(N_SNPS/1000))

# Build gamma model- no need for zero inflation, pretty normal amount of 0s
piglm <- glm(data = divfilt,
             PI ~ population,
      family = Gamma)

#check_model(methglm) # Looks lovely
check_predictions(piglm)
tidy(piglm)

# Contrast populations with emmeans
piEM <- emmeans(piglm, specs = ~population) |>
  contrast(method = "pairwise") |>
  confint()

piEM |> plot() +
  geom_vline(xintercept = 0, color = 'red', lty = 2)



# Tajima's D
divfilt |>
  ggplot(mapping = aes(x = TajimaD)) +
  geom_histogram() +
  facet_wrap(vars(population))

# Fit normal distribution lm
dlm <- lm(data = divfilt,
           TajimaD ~ population)

check_model(dlm)

# Contrast populations
dEM <- emmeans(dlm, specs = ~population) |>
  contrast(method = "pairwise") |>
  confint()

dEM |> plot() +
  geom_vline(xintercept = 0, color = 'red', lty = 2)


#### Plotting
# Plot Pi diversity- same way we did methylation
piVio <- divfilt |>
  ggplot(mapping = aes(x = population, y = PI, color = population)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.01) +
  geom_violin(fill = NA, linewidth = 1) +
  
  # Put it on a log scale so you can actually see the lower end
  scale_y_continuous(trans = 'log10') +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Pi Diversity in 10kb Windows",
       x = "Population",
       y = "Pi over 10kb window") +
  guides (size = "none") +
  theme_classic(base_size = 16)

# Boxplot Pi diversity
piBox <- divfilt |>
  ggplot(mapping = aes(x = population, y = PI, color = population)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.01) +
  geom_boxplot(fill = NA, linewidth = 1) +
  
  # Put it on a log scale so you can actually see the lower end
  scale_y_continuous(trans = 'log10') +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Pi Diversity over 10kb Windows",
       x = "Population",
       y = "Pi over 10kb window") +
  guides (size = "none") +
  theme_classic(base_size = 16)

png("6.diversity/piBox.png", width = 600, height = 800)
piBox
dev.off()


# Plot Tajima violin plot
dvio <- divfilt |>
  ggplot(mapping = aes(x = population, y = TajimaD, color = population)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.01) +
  geom_violin(fill = NA, linewidth = 1) +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Tajima's D in 10kb Windows",
       x = "Population",
       y = "D over 10kb window") +
  guides (size = "none") +
  theme_classic(base_size = 16)


# Plot Tajima boxplot
dbox <- divfilt |>
  ggplot(mapping = aes(x = population, y = TajimaD, color = population)) +
  geom_point(position = position_jitter(width = .3, seed = 10), size = 1, alpha = 0.01) +
  geom_boxplot(fill = NA, linewidth = 1) +
  
  # Make it pretty
  scale_color_manual(labels = c("Folger's Marsh", "Savin Hill Cove",
                                "The Creeks Preserve", "Waquoit Bay"),
                     values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  labs(title = "Tajima's D in 10kb Windows",
       x = "Population",
       y = "D over 10kb window") +
  guides (size = "none") +
  theme_classic(base_size = 16)

png("5.diversity/tajdBox.png", width = 600, height = 800)
dbox
dev.off()

