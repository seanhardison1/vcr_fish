library(tidyverse)
library(stars)
library(sf)
library(rgdal)
library(raster)
library(landscapetools)
library(landscapemetrics)
library(magrittr)

# import data
load(here::here("data/fish_processed_spatial.rdata"))
load(here::here("data/processed_phys_vars.rdata"))

sample_locs <- fish_sf %>% 
  dplyr::select(year, season, insideSeagrass, site) %>%
  distinct()

# reproject wind fetch raster so that
# resolutions match the DEM (50 m x 50 m)
fetch2015 <- raster::projectRaster(fetch2015a,r2)
fetch2014 <- raster::projectRaster(fetch2014a,r2)

# residence time is already reprojected
dem <- raster::projectRaster(dem,r2)

# average the two wind fetch rasters
fetch <- mosaic(fetch2014, fetch2015, fun = mean, na.rm = T)

# generate 150 m buffer around sampling sites
d <- 150
loc_buffs <- 
  sample_locs %>% 
  dplyr::select(site) %>% 
  distinct() %>% 
  st_buffer(., dist = d) %>% 
  as_Spatial()

# intersect buffers with sampling sites
nlist_ft <- extract(fetch, loc_buffs)
names(nlist_ft) <- unique(sample_locs$site)

nlist_rt <- extract(res_time, loc_buffs)
names(nlist_rt) <- unique(sample_locs$site)

nlist_depth <- extract(dem, loc_buffs)
names(nlist_depth) <- unique(sample_locs$site)

# combine all in data.frame
phys_vars <- NULL
for (i in 1:length(unique(sample_locs$site))){
  assign("phys_vars",
         rbind(phys_vars,
               tibble(depth = mean(nlist_depth[[unique(sample_locs$site)[i]]]),
                      ft = mean(nlist_ft[[unique(sample_locs$site)[i]]]),
                      rt = mean(nlist_rt[[unique(sample_locs$site)[i]]]),
                      site = unique(sample_locs$site)[i])))
}

save(phys_vars, file = here::here("data/phys_vars.rdata"))