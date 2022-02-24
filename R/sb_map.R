library(raster)
library(tidyverse)
library(sf)
library(magrittr)

# load VCR DEM polygon
load(here::here("data/vcr_land_polygon.rdata"))

## Load SB polygon (bounding box)
sb_bbox <- st_read(here::here("data/SouthBay.kml")) %>% 
  ## The .kml needs to be reprojected to match CRSs with the VCR DEM polygon
  sf::st_transform(crs = st_crs(rast_poly)) %>% 
  st_zm()

# intersect the two polygons
sb_poly <- rast_poly %>% 
  st_intersection(.,sb_bbox) %>% 
  st_as_sf() %>% 
  dplyr::rename(land = x)
  

ggplot() + 
  geom_sf(data = sb_poly,
          color = "transparent",
          fill = "#BDB76B") +
  coord_sf(xlim = c(-75.885, -75.795),
           ylim = c(37.195, 37.29))

