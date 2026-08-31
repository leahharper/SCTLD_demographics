library(tidyverse)
library(rgeos)
library(ggmap)
library(leaflet)
library(webshot2)
library(ggnewscale)
library(maps)
library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(dplyr)

# Change working directory
setwd("C:/Users/harperl/OneDrive - Smithsonian Institution/Documents/GitHub/SCTLD_demographics/")


# Import the survey data
coords <- read.csv("metadata.csv")

coords$Tagged.Corals.Y.N <- recode(coords$Tagged.Corals.Y.N,
                                   "Y" = TRUE,
                                   "N" = FALSE)

coords <- coords %>%
  mutate(marker_label = ifelse(Tagged.Corals.Y.N, "*", ""))

##############
#Map Sites
##############
SiteColor <- colorFactor(c("#882255", "#117733", "#DDCC77", "#332288", "#44AA99", "#999933", "#CC6677"),
                        domain = c("CBC Central", "CBC Lagoon", "CBC30N",
                                   "Curlew Patch", "House Reef", "South Reef Central", "SR30N"))
sitecolors = c('CBC Central'='#882255','CBC Lagoon'='#117733','CBC30N'='#DDCC77',
               'Curlew Patch' = '#332288', 'House Reef' ='#44AA99', 
               'South Reef Central' = '#999933', 'SR30N' = '#CC6677')

# A map of all CBC sites
CBCsites <- leaflet(coords) %>%
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  setView(lng = -88.09, lat = 16.79, zoom = 12.2) %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude, weight = 1, 
                   radius = 7, popup = ~Site,  fillColor = ~SiteColor(Site),
                   fillOpacity = 0.8, 
    label = ~marker_label,
    labelOptions = labelOptions(
      noHide = TRUE, # Always show the label
      direction = 'center', # Center the label over the marker
      textOnly = TRUE, # Ensure only the text is displayed
      style = list("font-size" = "24px", "font-weight" = "italic", "color" = "black"))) %>%
  addLegend("bottomright", SiteColor, values = ~Site,
            title = "Study Sites",
            opacity = 1) %>%
  addScaleBar("topright", options = scaleBarOptions(imperial = FALSE))

CBCsites



# 1. Save as an HTML widget file
htmlwidgets::saveWidget(CBCsites, "temp_map.html", selfcontained = FALSE)

# 2. Convert the HTML file to a static PNG image
webshot2::webshot("temp_map.html", "Figures/current/inset.png", delay = 0.2)

# Wider map


# high-resolution coastline/land polygons (10m scale - resolves small islands/cays)
land <- ne_countries(scale = 10, returnclass = "sf")

# optional: minor islands layer, which rnaturalearth stores separately
# and often includes cays/small islands missed by the main coastline layer
minor_islands <- ne_download(scale = 10, type = "minor_islands",
                             category = "physical", returnclass = "sf")

BZArea <- c(min(-85, na.rm = TRUE),
            min(20,  na.rm = TRUE),
            max(-85, na.rm = TRUE),
            max(20,  na.rm = TRUE))

# add a small buffer so points near the edge aren't clipped
buffer <- 11
xlim <- c(BZArea[1] - buffer, BZArea[3] + buffer)
ylim <- c(BZArea[2] - buffer, BZArea[4] + buffer)

tiff("Figures/current/Fig 1_wide.png", width = 3, height = 3, units = "in", res = 300)

ggplot() +
  geom_sf(data = land, fill = "grey85", color = "grey40", linewidth = 0.2) +
  geom_sf(data = minor_islands, fill = "grey85", color = "grey40", linewidth = 0.2) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  labs(x = NULL, y = NULL)

dev.off()

# 1. Build a minimal ggplot just to generate the legend you want
legend_plot <- ggplot(coords, aes(x = Longitude, y = Latitude, color = Site)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = sitecolors, name = "Study Sites") +
  theme_void() +  # strips all plot elements, keeps only the legend
  theme(
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white")
  )

# 2. Extract just the legend grob using cowplot
legend_grob <- cowplot::get_legend(legend_plot)

# 3. Save the legend alone as its own image
legend_only_plot <- ggdraw(legend_grob)
ggsave("Figures/current/legend_only.png", legend_only_plot, width = 2, height = 2, dpi = 300, bg = "white")

# 4. Read it back in as a magick image
legend_img <- image_read("Figures/current/legend_only.png") %>%
  image_border("black", "2x2")

library(cowplot)
library(magick)
wide_img <- image_read("Figures/current/Fig 1_wide.png")
inset_img <- image_read("Figures/current/inset.png") %>%
  image_resize("880x660")  # scale inset to desired size

# optional: add a border/frame around the inset for visual separation

inset_cropped <- image_crop(inset_img, "400x300+220+170") %>%
  image_border("black", "1x1")

inset_cropped
# composite: position inset in a corner (adjust offset as needed)


final_map <- ggdraw() +
  draw_image(wide_img) +
  draw_image(inset_cropped, x = 0.47, y = 0.1, width = 0.5, height = 0.5) +
  draw_line(
    x = c(0.4, 0.47),   # start proportion, end proportion (left point -> inset)
    y = c(0.36, 0.35),  # start proportion, end proportion
    color = "black", size = 1, linetype = "solid"
  )  +
  draw_image(legend_img, 
             x = 0.65, y = 0.585,      # position: proportion from left, bottom
             width = 0.4, height = 0.3)

final_map

ggsave("Figures/current/Fig 1.png", final_map, width = 10, height = 8, dpi = 300, bg = "white")
