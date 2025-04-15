#'''''''''''''''''''''''''''''''''''''''''''
#' GO Analysis
#' @date 2025-04-15
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''''''''''''

# Load Libraries
library(tidyverse)

# Load in background and significant gene lists
back <- read_csv("8.go/background.ps.csv", col_names = FALSE) |>
  rename(gene = X1,
         prot = X2)

goi <- read_csv("8.go/goi.ps.csv", col_names = FALSE) |>
  rename(gene = X1,
         prot = X2)

# Make summary tables for goi and background
backsum <- back |>
  group_by(prot) |>
  summarize(background = n()) |>
  arrange(desc(background))

goisum <- goi |>
  group_by(prot) |>
  summarize(divergent = n()) |>
  arrange(desc(divergent))

# merge and keep rows with divergent hits
both <- full_join(y = backsum, x = goisum, by = "prot") |>
  replace_na(list(background = 0, divergent = 0)) |>
  filter(divergent > 0)

summarize(both, totdiv = sum(divergent),
          totback = sum(background))

# Write a function to do Fisher's Exact Test on each category
fetest <- function(prot, divergent, background) {
  # NOTE: These matrix values are for all genes we have prot-scriber hits for
  # But the table only needs to keep the pathways we're interested in testing
  ft <- matrix(c(divergent, background, 70-divergent, 24367-background),
            nrow = 2,
            dimnames = list(type = c("divergent", "background"),
                            pathway = c("same", "other"))) |>
    fisher.test(alternative = "greater")
  
  data.frame(id = prot, p = ft$p.value, orat = ft$estimate)
}

fres <- pmap(.l = both, .f = fetest) |> list_rbind() |>
  remove_rownames() |>
  arrange(p) |>
  mutate(padj = p.adjust(p, method = "BH"))

write_tsv(fres, "8.go/FEtest.tsv")

