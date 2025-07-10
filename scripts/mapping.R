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


# Convert coords to lat long and filter out western MA
mass <- st_transform(mass, 4326)

massMap <- ggplot(data = mass) +
  geom_sf(color = "darkgrey")
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

#### Create inset metals plot
metals <- read_csv("data/xrf_soil_data_intensity_concentration.csv")

# Filter and rename sites
metFilt <- metals |>
  filter(element == "Pb",
         site %in% c("FM", "WB", "SHC", "TC")) |>
  mutate(site = str_replace_all(site,
                                      c("FM" = "Folger's\nMarsh",
                                        "SHC" = "Savin Hill\nCove",
                                        "TC" = "The Creeks\nPreserve",
                                        "WB" = "Waquoit\nBay")))

# Order factor to change boxplot series
metFilt$site <- factor(metFilt$site, levels = c("Folger's\nMarsh", "The Creeks\nPreserve",
                                                "Waquoit\nBay", "Savin Hill\nCove"))
leadBox <- metFilt |>
  ggplot(mapping = aes(x = site, y = concentration, color = site)) +
  geom_boxplot() +
  scale_color_manual(labels = c("Folger's Marsh", "The Creeks Preserve",
                                "Waquoit Bay", "Savin Hill Cove"),
                     values = c("#D81B60", "#FFC107", "#004D40", "#1E88E5")) +
  labs(title = "Lead ppm in Study Marshes") +
  guides (size = "none", color = "none") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

png(filename = "report/lead.png", width = 400, height = 400)
leadBox
dev.off()
