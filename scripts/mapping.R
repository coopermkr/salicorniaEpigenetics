#''''''''''''''''''''''''''''''''''''''''''''''''
#' 
#' Mapping Sample Sites
#' @date 2025-06-18
#' @author Cooper Kimball-Rhines
#' 
#''''''''''''''''''''''''''''''''''''''''''''''''


# Load libraries
library(maptools)
library(tidyverse)
library(sf)
library(ggrepel)

# Load MA shapefile
mass <- st_read("mapping/outline25k/OUTLINE25K_ARC.shp")

# Convert coords to lat long
mass <- st_transform(mass, 4326)

massMap <- ggplot(data = mass) +
  geom_sf()
massMap

# Load in site coords
sites <- read_csv("mapping/site_coordinates.csv") |>
  filter(site %in% c("The Creeks Preserve", "Folger's Marsh", "Waquoit Bay", "Savin Hill Cove"))

# Make sample site map
siteMap <- massMap +
  geom_point(data=sites, aes(x=long, y=lat), color=c("#004D40", "#FFC107", "#D81B60", "#1E88E5"), size=2.5, show.legend = FALSE) +
  geom_label_repel(data=sites, aes(x=long, y=lat, label=site), nudge_y = 0, nudge_x = -0.03) +
  theme_classic(base_size = 16) +
  xlab("Longitude") + ylab("Latitude")
siteMap

png(filename = "report/sampleMap.png", width = 600, height = 400)
siteMap
dev.off()