#''''''''''''''''''''''''''''''''''''''''''''''''
#' Methylation Window Density
#' @date 2025-03-05
#' @author Cooper Kimball-Rhines
#''''''''''''''''''''''''''''''''''''''''''''''''

library(tidyverse)
library(furrr)
library(data.table)

meth <- read_tsv("data/filtered.all.bed") |>
  mutate(pop = substr(id, 1,2))

#### Calculate window methylation density
# Calculate windows around genes
# genes <- fread("filtered.long.bed") |>
#   filter(stringr::str_starts(V1, "^group")) |>
#   rename(chrom = V1,
#          start = V2,
#          end = V3,
#          name = V4,
#          score = V5,
#          strand = V6,
#          thicks = V7,
#          thicke = V8,
#          rgb = V9,
#          bcount = V10,
#          bsizes = V11,
#          bstart = V12)

#windows <- data.frame(winname = genes$name,
#                      winstrand = genes$strand,
#                      winstart = genes$start - 5000,
#                      winend = genes$end + 5000,
#                      winchrom = str_split_i(string = genes$chrom, pattern = "_", i = 1))


# OR Calculate standard size windows
sizes <- read.delim("data/sizes.tetra.scaff.18.txt", header = FALSE, sep = " ") |>
  filter(stringr::str_starts(V1, "^group")) |>
  rename(chrom = V1, length = V2)

wincalc <- function(chr, len) {
  print(chr)

  sizes |>
    reframe(winchrom = chr,
            winstart = as.numeric(seq(from = 1, to = len, by = 10000)))
}

windows <- map2(.x = sizes$chrom, .y = sizes$length, .f = wincalc) |> list_rbind()

# Calculate methylation densities in each window
windens <- function(chr, winstart) {
  print(chr)
  winend = winstart + 10000
  
  meth |>
    filter(chrom == !!chr,
           start >= winstart,
           end < winend) |>
    group_by(chrom, id) |>
    summarize(winchrom = chr,
              window = winstart,
              winmeth = sum(Nmod),
              wincov = sum(cov),
              winsites = n())
}

methTest <- meth |>
  select(chrom, id, start, end, Nmod, cov) |>
  filter(chrom == "group1")

# Run without parallelization
system.time(dens1unpar <- map2(.x = "group1", .y = windows$winstart[1:50], .f = windens) |> list_rbind())
# 50 calculations takes 64.52   14.53   79.52 seconds

# Fix data size error- only necessary on Windows
options(future.globals.maxSize = 1e11)

# Run the function in parallel to calculate methylation densities
plan(multisession, workers = 4)
system.time(dens1 <- future_map2(.x = "group1", .y = windows$winstart[1:50], 
                                 .f = windens, .progress = TRUE) |> list_rbind())

# 50 calculations with 4 cores takes 75.69  222.25  529.05
# check that the parallel and unparallel output are the same
identical(dens1, dens1unpar) #TRUE
# The starting dataset is so large that the time it takes to copy it to workers
# outweighs the time saved from faster calculations
# stick with the unparallelized version and loop it over many jobs or figure out a way
# to chunk the input so it's faster to calculate. Maybe one core per chromosome?

# Calculate number of populations with data per window
winpops <- function(chr, winstart) {
  print(chr)
  winend = winstart + 10000
  
  meth |>
    filter(chrom == !!chr,
           start >= winstart,
           end < winend) |>
    summarize(name = windows[i,1],
              npop = n_distinct(pop))
}

denpop <- map2(.x = windows$chr, .y = windows$winstart, .f = windpops) |> list_rbind()

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
      chrom == windows[i,5],
      start > windows[i,3],
      end < windows[i,4])
  
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

write_tsv(winsum, "5.dma/winsum.tsv")
write_tsv(winpop, "5.dma/winpop.tsv")

