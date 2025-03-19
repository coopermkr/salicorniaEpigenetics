#''''''''''''''''''''''''''''''''''''''''''
#' Epigenetic DMA
#' @date 2025-03-05
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(data.table)

#### One Chromosome DMA
# load methylation and annotation data
meth <- fread("data/filtered.all.bed") |>
  mutate(pop = substr(id, 1,2))

genes <- fread("../Genome/Report/filtered.long.bed") |>
  filter(stringr::str_starts(V1, "^group1_")) |>
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

# Calculate windows
windows <- data.frame(winname = genes$name,
                      winstrand = genes$strand,
                      winstart = genes$start - 5000,
                      winend = genes$end + 5000,
                      winchrom = str_split_i(string = genes$chrom, pattern = "_", i = 1))

# Establish dataframe by running first window calculation
filt <- meth |>
  filter(chrom == windows[1,5],
         start > windows[1,3],
         end < windows[1,4])

winsum <- filt |>
  group_by(id) |>
  summarize(name = windows[1,1],
            winmeth = sum(Nmod),
            wincov = sum(cov),
            winsites = n())

winpop <- filt |>
  group_by(pop) |>
  summarize(name = windows[1,1],
            npop = n_distinct(pop))

# Loop to calculate methylation in all windows, saving to the dataframe above
for (i in 2:nrow(windows)) {
  print(i)
  
  filt <- meth |>
    
    # Filter to the proper gene location
    # We may not actually care about strand here--hemimethylation still matters
    filter(#strand == windows[i,2],
           start > windows[i,3],
           end < windows[i,4],
           chrom == windows[i,5])
    
    # summarize number of methylated and total reads and sites
    winsum <<- filt |>
          group_by(id) |>
          summarize(name = windows[i,1],
              winmeth = sum(Nmod),
              wincov = sum(cov),
              winsites = n()) |>
    
    # Merge with last loop iteration
    merge(winsum, all = TRUE)
  
    # And count up distinct populations for each gene
  winpop <<- filt |>
    summarize(name = windows[i,1],
              npop = n_distinct(pop)) |>
    merge(winpop, all = TRUE)
  
  # Increase loop by one
  i <- i + 1
}

# Note, you only return windows that you have data for, otherwise the filter returns
# a zero row frame for summarize to work with. We're not interested in those anyway
library(betareg)
library(lmtest)

# Merge together the ndistinct population info into the winsum table
winper <- winpop |> 
  select(!pop) |>
  merge(winsum) |>
  
  # filter for genes that have info from all four pops
  filter(npop > 3) |>
  mutate(pop = substr(id, 1, 2)) |>
  
  # Calculate methylation percentage across the window for each sample
  mutate(methper = winmeth/wincov,
         # and transform it so we don't end up with 0s or 1s
         methtrans = (methper*(wincov-1)+0.5)/wincov)

# Model the effect of population on methylation as a betareg
winper |>
  filter(name == "anno1.g13784.t1") |>
  betareg(formula = methtrans ~ pop,
                  weights = wincov,
                  link = "logit") |>
  waldtest()





