#### R script for submission in a SLURM array

library(tidyverse)

setwd("/work/pi_brook_moyers_umb_edu/salicornia/epigenetics/")

# Get chromosome from task ID in the array
task_id <- commandArgs(trailingOnly=TRUE)[1] |> as.numeric()
#loop <- "Sdep1A"
#loop <- paste("Sdep", task_id, "A", sep = "")
loop <- paste("Sdep", task_id, "B", sep = "")
paste(loop)

# Load annotation
feature <- "promoter"

# anno <- read_tsv("data/Sdep_NCBI_annotation.gff3", col_names = FALSE, skip = 1) |>
#   rename("chrom" = X1,
#          "type" = X3,
#          "start" = X4,
#          "end" = X5,
#          "strand" = X7,
#          "name" = X9) |>
#   filter(chrom == loop,
#          type == feature) |>
#   select(chrom, start, end)

# OR Load Nitrogen differential genes
anno <- read_tsv("data/detectedTrans.tsv", col_names = TRUE) |>
    rename("chrom" = Chr,
           "start" = Start,
           "end" = End,
           "strand" = Strand,
           "name" = Geneid) |>
  mutate(chrom = substr(chrom, 1, 6),
         start = as.numeric(str_split_i(start, ";", 1)),
         end = as.numeric(str_extract(end, "[^;]+$"))) |>
# New line to calculate promoter start and end based on strand
  mutate(end = start,
	 start = ifelse(strand == "+", start - 3000, start + 3000)) |>
    filter(chrom == loop) |>
    select(chrom, start, end)

# Load in data and filter by task
methBed <- read_tsv("4.filter/filtered.2m.bed") |>
  select(chrom, id, start, end, Nmod, cov) |>
  filter(chrom == loop) |>
  mutate(pop = substr(id, 1,2))

# Calculate methylation densities for each feature
# We filtered by chromosome above, so we just need to filter by feature start and end position
methyldens <- function(featstart, featend) {
  methBed |>
    filter(start >= featstart,
           end <= featend) |>
    group_by(chrom, id) |>
    summarize(feat = featstart,
              end = featend,
              winmeth = sum(Nmod),
              wincov = sum(cov),
              winsites = n())
}

featdens <- map2(.x = anno$start, .y = anno$end, .f = methyldens) |> list_rbind()

# Save methylation density info
write_csv(x = featdens, file = paste("6.features/featdens_", loop, "_", feature, ".csv", sep = ""))

# Calculate number of populations with data per window
winpops <- function(featstart, featend) {
  methBed |>
    filter(start >= featstart,
           end <= featend) |>
    summarize(name = featstart,
              end = featend,
              npop = n_distinct(pop))
}

popdens <- map2(.x = anno$start, .y = anno$end, .f = winpops) |> list_rbind()

# Save population info
write_csv(x = popdens, file = paste("6.features/popdens_", loop, "_", feature, ".csv", sep = ""))

meth <- popdens |>
  rename(feat = name) |>
  merge(featdens) |>
  filter(npop > 2) |>
  mutate(pop = substr(id, 1, 2),
         feat = as.numeric(feat),
         methper = winmeth/wincov,
         # and transform it so we don't end up with 0s or 1s
         methtrans = (methper*(wincov-1)+0.5)/wincov,
	 grp = paste(chrom, feat, sep = "_"))

# Test for the effect of population on methylation density
kwtest <- function(start) {
  tmp <- meth |>
    filter(feat == start) %>%
    kruskal.test(methtrans ~ pop, data = .)
  
  data.frame(chrom = loop,
             feat = start,
             p = tmp$p.value)
}

kw <- map(.x = unique(meth$feat), .f = kwtest) |> list_rbind()

kw <- meth |> select(chrom, feat, end) |> distinct() |> merge(kw)

write_csv(x = kw, file = paste("6.features/kw_", loop, "_", feature, ".csv", sep = ""))

print("done")
