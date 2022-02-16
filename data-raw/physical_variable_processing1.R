library(tidyverse)
library(raster)
library(mgcv)
library(nngeo)
library(magrittr)
library(glmmTMB)
library(sf)
library(readxl)
process_raw <- F

#raw data
if (process_raw){
  vcr <- st_read(here::here("data/VCR_poly3.kml")) %>% 
    sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
    st_zm() %>% 
    st_as_sf()
  
  dem_raw <- raster(here::here("data/V2_topobathy_tif.tif"))
  sb_raw <- raster(here::here("data/SouthBay_Bathymetry_Shareable/SB_Bathy_krig/SB1_Krig1.tif"))
  
  # residence time is a large file and need to be aggregated before pushed to github
  res_time <- raster(here::here("data/dxv_blended/w001001.adf"))
  
  # template raster for reprojecting data
  r2 <- raster(crs = crs(vcr),
               ext = extent(vcr),
               resolution = 50)
  res_time <- raster::projectRaster(res_time,r2)
  
  fetch2015a <- raster(here::here("data/Summers2014_2015_Fetch/combined15a.tif"))
  fetch2014a <- raster(here::here("data/Summers2014_2015_Fetch/combined14a.tif"))
  
  save(r2, vcr, dem_raw, sb_raw, res_time, fetch2014a, fetch2015a,
       file = here::here("data/raw_physical_variables.rdata"))
} else {
  load(here::here("data/raw_physical_variables.rdata"))
}

# VCR Polygon for clipping raster
vcr_poly <- st_read(here::here("data/VCR_poly3.kml")) %>% 
  sf::st_transform(crs = st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>% 
  st_zm() %>% 
  as(., "Spatial")

# Cropped DEM
dem <-
    mask(
      crop(
        dem_raw,y = extent(vcr_poly)
      ), 
      vcr_poly)

# Cropped residence time raster
res_time <- 
    mask(
      crop(
        res_time,y = extent(vcr_poly)
      ), 
      vcr_poly)

#here I remove all values in the DEM that are > 0 m elevation
dem[dem > 0] <- NA

# the DEM is a composite of the DEM found here: http://www.vcrlter.virginia.edu/cgi-bin/showDataset.cgi?docid=knb-lter-vcr.210
# and a second DEM created specifically for South Bay by Tyler Barnes (UVA). 
# The following code reprojects the older DEM and merges it with the new data. The merge step
# simply replaces the overlapping pixels within the two DEMs with the newer data.
dem <- raster::projectRaster(dem, crs = crs(sb_raw),
                             res = res(sb_raw))
origin(dem) <- origin(sb_raw)
dem <- merge(sb_raw, dem) # 
dem <- raster::projectRaster(dem, crs = crs(dem_raw))
save(dem, res_time, r2,
     fetch2014a, fetch2015a, 
     file = here::here("data/processed_phys_vars.rdata"))

