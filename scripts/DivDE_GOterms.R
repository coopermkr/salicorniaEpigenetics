#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' S. depressa DE and Divergent Biological Processes
#' @date 2025-12-01
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Option 1: an upset plot
library(tidyverse)
library(UpSetR)

# Load in genetic and methylation divergent processes
marshSet <- read_delim(file = "8.go/marshDiverGO.txt", skip = 7, col_names = FALSE) |>
  rename("BiologicalProcess" = X1,
         "ArabRef" = X2,
         "marshGenes" = X3,
         "Expected" = X4,
         "Enrichment" = X5,
         "Reg" = X6,
         "pValue" = X7) |>
  mutate(MarshDiv = case_when(marshGenes > 0 ~ 1)) |>
  select(BiologicalProcess, MarshDiv, marshGenes)

# Load in N experiment processes
NSet <- read_delim(file = "8.go/allTreatDEGO.txt", skip = 7, col_names = FALSE) |>
  rename("BiologicalProcess" = X1,
         "ArabRef" = X2,
         "NGenes" = X3,
         "Expected" = X4,
         "Enrichment" = X5,
         "Reg" = X6,
         "pValue" = X7) |>
  mutate(ExpGenes = case_when(NGenes > 0 ~ 1)) |>
  select(BiologicalProcess, ExpGenes, NGenes)

# Merge process sets
goSets <- merge(marshSet, NSet, all = TRUE)
goSets[is.na(goSets)] <- 0 # Replace NAs with 0s

upset(goSets, sets = c("ExpGenes", "MarshDiv"))

#### Option 2: Venn diagram ####
library(ggVennDiagram)
goVenn <- goSets |>
  mutate(ExpGenes = ExpGenes*row_number(),
         marshDiv = MarshDiv*row_number())

# make a numbered list of genes for each treatment
n <- goVenn |> filter(ExpGenes != 0) |> select(ExpGenes)
m <- goVenn |> filter(marshDiv != 0) |> select(marshDiv)

x <- list(Experimental = n$ExpGenes, Marsh = m$marshDiv)

# Make the Venn diagram
goVennPlot <- ggVennDiagram(x) +
  theme(legend.position = "none") +
  labs(title = "Overrepresented Biological Processes") +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  coord_flip()

jpeg(filename = "8.go/processOverlap", width = 600, height = 600, quality = 100)
goVennPlot
dev.off()

# Is the overlap significant?
fe <- matrix(c(7, 87, 65, 24994),
             nrow = 2,
             dimnames = list(Marsh = c("Sig", "Insig"),
                             Exp = c("Sig", "Insig")))
fisher.test(fe)
