library(raster)
library(tidyverse)
library(sf)
library(magrittr)
library(rnaturalearth)

process_dem <- F
# Load packages
load(here::here("data/sav_cover.rdata"))

# sampling sites in 2012
load(here::here("data/fish_processed_spatial.rdata"))
sample_locs <- 
  fish_sf %>% 
  filter(year == 2012) %>% 
  dplyr::select(geometry, site, insideSeagrass) %>% 
  distinct()

# same CRS for all operations
ncrs <- "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs"

# coastline
shoreline <- ne_countries(scale='large',returnclass = 'sf')

## Load VCR polygon to trim DEM
vcr_poly <- st_read(here::here("data/vcr_poly3.kml")) %>% 
  ## The shapefile needs to be reprojected to match CRSs with the DEM
  sf::st_transform(crs = st_crs(ncrs)) %>% 
  st_zm() %>% 
  as(., "Spatial")

# process sample sites
sample_locs.df <- sample_locs %>% 
  st_transform(st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>% 
  dream::sfc_as_cols(.,names = c("Longitude","Latitude")) %>% 
  st_set_geometry(NULL) %>% 
  distinct() %>% 
  mutate(`Sampling Habitat` = ifelse(!str_detect(site, "Bare"),
                                     "Seagrass",
                                     "Unvegetated"),
         fill = ifelse(`Sampling Habitat` == "Seagrass",
                       "#1F77B4FF",
                       "#FF7F0EFF"))

# boundinng box for DEM
ymin <- 37.175
xmin <- -75.95
ymax <- 37.5
xmax <- -75.6
VCR_bb <- extent(c(xmin = xmin, xmax = xmax,
                   ymin = ymin, ymax = ymax))

dem_bbox <- st_bbox(VCR_bb,
                    crs = 4326) %>% 
  st_as_sfc() %>% 
  st_transform(.,crs = st_crs(ncrs)) %>% 
  as_Spatial()

if (process_dem){
  dem_raw <- raster("/Users/seanhardison/Documents/git/seagrass-landscape-metrics/data/V2_topobathy_tif.tif")
  # process DEM----
  reduced_dem <- aggregate(mask(crop( dem_raw,y = dem_bbox), dem_bbox),
                           fact = 5)
  plot(reduced_dem)
  land_dem <- reduced_dem
  # Remove all values in the DEM that are > 0 m elevation
  land_dem[land_dem < -0] <- NA
  plot(land_dem)
  
  # convert DEM to polygon
  sf_use_s2(FALSE)
  rast_poly <- sf::st_as_sf(stars::st_as_stars(land_dem),
                            as_points = FALSE, merge = TRUE) %>%
    st_union()
  rast_poly %<>% st_transform(st_crs("+proj=longlat +datum=WGS84 +no_defs"))
  save(rast_poly, file = here::here("data/vcr_land_polygon.rdata"))
} else {
  load(here::here("data/vcr_land_polygon.rdata"))
}


#convert raster to df for plotting ocean
sea_rast <- dem_reproj
sea_rast[sea_rast>0] <- NA

xmin2 <- -78
xmax2 <- -73
ymin2 <- 33
ymax2 <- 41

EC_bb <- extent(c(xmin = xmin2, xmax = xmax2,
                  ymin = ymin2, ymax = ymax2))

EC_bbox <- st_bbox(EC_bb,
                   crs = 4326) #%>% 
  # st_as_sfc() %>% 
  # st_transform(.,crs = st_crs(ncrs))

shoreline_vis <- shoreline %>% 
  st_crop(.,EC_bbox) %>% 
  summarise() %>% 
  st_transform(st_crs("+proj=longlat +datum=WGS84 +no_defs"))

# visualize -------------------------------------------------------------------------
land <- "#BDB76B"

sav_vis <- sav %>% 
  st_crop(.,xmin = xmin, ymin = ymin, ymax = ymax, xmax = xmax) %>%
  dplyr::filter(DENSITY != 0, year %in% c(2012, 2018)) %>%
  group_by(year) %>% 
  summarise() %>% 
  mutate(year = factor(year, levels = c(2012, 2018))) %>% 
  st_transform(st_crs("+proj=longlat +datum=WGS84 +no_defs"))

ec_zoom <- ggplot() +
  geom_sf(data = shoreline_vis,
          size = 0.5,
          fill = "#BDB76B",
          color = "transparent") +
  geom_sf(data = dem_bbox %>%  
            as(.,"sf") %>% 
            st_transform(st_crs("+proj=longlat +datum=WGS84 +no_defs")),
          fill = "transparent",
          color = "black",
          size = 0.5) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent", linetype = 1),
        panel.background = element_rect(fill = "#E0FFFF", color = NA),
        panel.ontop = F,
        panel.border = element_rect(size = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  coord_sf(xlim = c(xmin2, xmax2),
           ylim = c(ymin2 + 1, ymax2 - 1)) +
  scale_x_continuous(expand = c(0.001,0.001)) +
  scale_y_continuous(expand = c(0.001,0.001))

site_map <- ggplot() + 
  geom_sf(data = rast_poly,
          color = "transparent",
          fill = land,
          inherit.aes = F) +
  
  geom_sf(data = sav_vis,
          aes(fill = year),
          alpha = 0.5,
          color = "transparent") + 
  geom_sf(data = sample_locs.df, aes(color = `Sampling Habitat`),
             size = 2.2) +
  geom_sf(data = sample_locs.df, aes(shape = `Sampling Habitat`),
             color = "black",
             size = 2.25, alpha = 0.5) +
    coord_sf(xlim = c(-75.9, -75.69),
             ylim = c(37.2, 37.45)) +
  
  scale_x_continuous(expand = c(0.001,0.001)) +
  scale_y_continuous(expand = c(0.001,0.001)) +
  scale_fill_manual(values = c("purple","grey30")) +
  scale_shape_manual(name = "Initial Sampling\n Habitat",
                     labels = c("Seagrass",
                                "Unvegetated"),
                     values = c(21,21)) +
  scale_color_manual(values = c("Seagrass" = "#1F77B4FF",
                                "Unvegetated" = "#FF7F0EFF")) +
  
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent", linetype = 1),
        panel.background = element_rect(fill = "#E0FFFF", color = NA),
        panel.ontop = F,
        panel.border = element_rect(size = 0.5),
        legend.key=element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.position = c(1.3, 0.3),
        axis.title = element_blank()) +
  
  labs(fill = "Seagrass cover",
       color = "Initial Sampling\n Habitat") +
  ylab("Latitude") +
  xlab("Longitude")
# site_map

sb_inset <- 
  ggplot() + 
  geom_sf(data = rast_poly,
          color = "transparent",
          fill = land,
          inherit.aes = F) +
  geom_sf(data = sav_vis,
          aes(fill = year),
          alpha = 0.5,
          color = "transparent") + 
  geom_point(data = sample_locs.df, aes(x = Longitude, y = Latitude,
                                        color = `Sampling Habitat`),
             size = 2) +
  coord_sf(xlim = c(-75.87,-75.845),
           ylim = c(37.235, 37.25)) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_line(color = "transparent", linetype = 1),
        panel.background = element_rect(fill = "#E0FFFF", color = NA),
        panel.ontop = F,
        panel.border = element_rect(size = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0.001,0.001)) +
  scale_y_continuous(expand = c(0.001,0.001)) +
  scale_color_manual(values = c("Seagrass" = "#1F77B4FF",
                                "Unvegetated" = "#FF7F0EFF")) +
  scale_fill_manual(values = c("purple","grey30")) +
  guides(shape = "none", color = "none", fill = "none")

inset_map <- ggplot() +
  coord_equal(xlim = c(0, 32), ylim = c(0, 28), expand = FALSE) +
  cowplot::draw_plot(ec_zoom, x = 25.75, y = 14.95, 
                     width = 10, height = 10) +
  cowplot::draw_plot(sb_inset, x = -4, y = 2.5, 
                     width = 10, height = 10) +
  annotation_custom(ggplotGrob(site_map), xmin = 0, xmax = 32, ymin = 0,
                    ymax = 28) +
  geom_segment(aes(x = 18.6, y = 6,
                   xend = 15, yend = 5),
               arrow = arrow(length = unit(0.15, "inches")), lineend = "butt",
               linejoin = "mitre") +
  geom_text(data = landmarks,
            aes(x = x,
                y = y ,
                label = name),
            color = "grey20",
            size = 5.5)+
  theme_void()

inset_map

ggsave(inset_map, 
       width = 10.5,
       height = 7.36,
       filename = here::here("figures/VCR_map2.png"),
       dpi = 300)
