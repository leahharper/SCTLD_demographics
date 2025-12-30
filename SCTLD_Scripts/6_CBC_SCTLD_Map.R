library(tidyverse)
library(rgeos)
library(ggmap)
library(leaflet)
library(webshot2)

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
SiteColor <- colorFactor(c("red3", "darkorchid", "gold", "blue", "seagreen4", "orange", "pink"),
                        domain = c("CBC Central", "CBC Lagoon", "CBC30N",
                                   "Curlew Patch", "House Reef", "South Reef Central", "SR30N"))


# A map of all CBC sites
CBCsites <- leaflet(coords) %>%
  addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
  setView(lng = -88.09, lat = 16.79, zoom = 12.2) %>%
  addCircleMarkers(lng = ~Longitude, lat = ~Latitude, weight = 1, 
                   radius = 7, popup = ~Site,  fillColor = ~SiteColor(Site),
                   fillOpacity = 0.6, 
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
webshot2::webshot("temp_map.html", "Figures/current/Fig 1.png", delay = 0.2)

