library(tidyverse)
library(magrittr)
library(here)
library(GGally)
library(rfishbase)
library(sf)
library(raster)
library(lubridate)
library(readxl)

###### GENERAL DATA CLEANING --------------------------------------------------------------------------
#Load raw data--------------------------------------------------------
data.dir <- here::here("data")
helper.dir <- here::here("R/helper")

process <- F

#Load fish seine survey database
fish_in <- readr::read_csv(file.path(data.dir,"FishQueryTable.csv"))

#Add species names to data while retaining common names
fish_classed <- 
  fish_in %>% 
  mutate(speciesName = tolower(speciesName)) %>% 
  mutate(sciName = case_when(speciesName %in% c("silversides",
                                                "atlantic silverside") ~ "Menidia menidia",
                             speciesName == "oyster toad" ~ "Opsanus tau",
                             speciesName == "pipefish" ~ "Syngnathus spp.",
                             speciesName == "seahorse"~"Hippocampus erectus",
                             speciesName == "mojarra"~"Eucinostomus argenteus",
                             speciesName == "goby"~"Gobiidae",
                             speciesName == "sand mullet"~"Mugilidae",
                             speciesName == "squid"~"Lolliguncula brevis",
                             speciesName == "unknown sciaenid"~"Sciaenidae",
                             speciesName == "filefish"~"Stephanolepis hispidus",
                             speciesName == "pufferfish"~"Sphoeroides maculatus",
                             speciesName %in% c("american anchovy",
                                                "bay anchovy",
                                                "unknown anchovy",
                                                "unknown bay anchovy")~"Anchoa spp.",
                             speciesName %in% c("croaker",
                                                "atlantic croaker") ~ "Micropogonias undulatus",
                             speciesName == "speckled trout"~"Cynoscion nebulosus",
                             speciesName == "striped blenny"~"Chasmodes bosquianus",
                             speciesName == "black sea bass"~"Centropristis striata",
                             speciesName %in% c("halfbeak",
                                                "american halfbeak")~ "Hyporhamphus meeki",
                             speciesName == "mummichog" ~ "Fundulus heteroclitus",
                             speciesName == "pinfish" ~ "Lagodon rhomboides",
                             speciesName %in% c("perch",
                                                "silver perch") ~ "Bairdiella chrysoura",
                             speciesName == "speckled seatrout" ~ "Cynoscion nebulosus",
                             speciesName == "striped burrfish" ~ "Chilomycterus schoepfii",
                             speciesName == "tautog" ~ "Tautoga onitis",
                             speciesName == "spot" ~ "Leiostomus xanthurus",
                             speciesName == "northern sennet"~"Sphyraena borealis",
                             speciesName == "pigfish" ~ "Orthopristis chrysoptera",
                             speciesName == "scup"~"Stenotomus chrysops",
                             speciesName == "bluespotted cornetfish" ~ "Fistularia commersonii",
                             speciesName == "sheepshead" ~ "Archosargus probatocephalus",
                             speciesName == "butterfish" ~ "Peprilus triacanthus",
                             speciesName == "summer flounder" ~ "Paralichthys dentatus",
                             speciesName == "mangrove snapper" ~ "Lutjanus griseus",
                             speciesName == "bluefish" ~ "Pomatomus saltatrix",
                             speciesName == "atlantic needlefish" ~ "Strongylura marina",
                             speciesName == "menhaden"~"Brevoortia tyrannus",
                             speciesName == "skilletfish" ~ "Gobiesox strumosus")) 

#Add lat/lon coordinates----------------------------------------------------
#There are two sources of sampling locations that overlap
# fish_sites1 <- readxl::read_excel(file.path(data.dir, "2012_18_fish_sites.xlsx"))  %>% 
#   dplyr::rename(longitude = long, latitude = lat) %>% 
#   st_as_sf(.,coords = c("longitude","latitude"), crs = 4326) %>% 
#   st_transform("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs") %>% 
#   dream::sfc_as_cols(names = c("longitude","latitude")) %>% 
#   st_set_geometry(NULL) %>% 
#   filter(site %in% unique(fish_in$site))
# 
# fish_sites2 <- st_read(file.path(data.dir, "Syn_Site_Points_2018")) %>% 
#   st_transform("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs") %>% 
#   dplyr::select(site = SiteName) %>% 
#   mutate(site = str_replace_all(site, " |-", "_")) %>% 
#   dream::sfc_as_cols(names = c("longitude","latitude")) %>% 
#   sf::st_set_geometry(NULL) %>% 
#   filter(site %in% unique(fish_in$site),
#          !site %in% unique(fish_sites1$site))

fish_sites <- 
  read_csv("/users/seanhardison/downloads/Historical_Fish_Sampling_Coordinates.csv") %>% 
  st_as_sf(.,coords = c("X","Y"), crs = 4326) %>%
  st_transform(st_crs("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")) %>%
  dream::sfc_as_cols(names = c("longitude","latitude")) %>% 
  sf::st_set_geometry(NULL) %>% 
  dplyr::select(site = Name, longitude , latitude ) %>% 
  mutate(site = str_replace_all(site, " ", "_")) %>% 
  dplyr::filter(!site %in% c("SB_105","SB_152",
                      "SB_Bare_4", "SB_Bare_5"))
# 
# t <- 
#   fish_in %>%
# left_join(
# fish_sites %>% 
#   mutate(site = str_replace_all(site, " ", "_")) 
# ) %>% 
#   dplyr::select(site, longitude, latitude) %>% 
#   distinct() 
# 
# 
# fish_sites <- bind_rows(fish_sites1, fish_sites2)
# all(unique(t$site) %in% unique(fish_in$site))

#join classified fish data with site coordinate info
fish1 <- 
  fish_classed %>% 
  left_join(.,fish_sites, by = "site") %>% 
  as.data.frame()

zeros <- fish1 %>% 
  filter(Count == 0)

fish1 %<>% filter(Count != 0) 

#expand rows where Count > 1--------------------------------------------------
fish <- NULL
for (i in 1:nrow(fish1)){
  
  df <- fish1 %>% dplyr::slice(i)
  if (df$Count > 1 & !is.na(df$Count)){
    ndf <- df[rep(row.names(df), df$Count), ]
  } else if (is.na(df$Count)){
    ndf <- df
  } else {
    ndf <- df
  }
  assign('fish', rbind(fish, ndf))
}

#Add seasonal categories----------------------------------------------------------
fish %<>% 
  dplyr::select(-Count) %>% 
  mutate(sampleDate = lubridate::as_date(sampleDate, format = "%m/%e/%y"),
         season = case_when(month(sampleDate) %in% c(5,6) ~ "summer",
                            month(sampleDate) %in% c(9,10) ~ "fall")) 

#Fill with zeros for species not caught-------------------------------------------
#retaining lengths
f_lengths <- 
  fish %>% 
  dplyr::select(sampleDate, site, sampleTime, sciName, Length, latitude, longitude) 

fish %<>% 
  dplyr::select(-Length) %>% 
  group_by(season, sampleDate, sampleTime, site, 
           latitude, longitude, 
            dissolvedOxygen, Temperature,
           salinity, conductivity, sciName) %>% 
  dplyr::summarise(nFish = n()) %>% 
  
  #spread and gather data to fill in missing species
  tidyr::spread(key = sciName, value = nFish, fill = 0) %>% 
  tidyr::gather(key = "sciName", value = "nFish", -c(season:conductivity)) %>% 

  #join lengths back in. rows with 0 fish caught have NA values for length
    left_join(.,f_lengths, by = c("sampleDate","sampleTime","site",
                                          "latitude","longitude","sciName")) %>% 
    group_by(season, sampleDate, sampleTime, site, 
             latitude, longitude, 
             dissolvedOxygen, Temperature,
             salinity, conductivity, sciName) %>%
  {. ->> fish_length_df} %>% 
  #nFish was the same for each row in sciName grouping, so mean(nFish) is just 
  #collapsing the row
  dplyr::summarise(nFish = mean(nFish),
                   meanLength = mean(Length),
                   medianLength = median(Length)) %>% 
  as.data.frame() 


# Join in zero counts for each species in the list
z1 <- zeros %>% 
  mutate(sampleDate = lubridate::as_date(sampleDate, format = "%m/%e/%y"),
         season = case_when(month(sampleDate) %in% c(5,6) ~ "summer",
                            month(sampleDate) %in% c(9,10) ~ "fall")) %>% 
  group_by(sampleDate, site) %>% 
  left_join(.,tibble(speciesName = unique(fish1$sciName),
                     Count = 0), by = "Count") %>%
  dplyr::select(sampleDate:conductivity,
                Length:season, 
                sciName = speciesName.y,
                -sciName) %>% 
  mutate(nFish = 0,
         meanLength = 0,
         medianLength = 0)

fish %<>% 
  bind_rows(.,z1)

#Remove outliers in water quality data----
fish %<>% 
  mutate(dissolvedOxygen = ifelse(dissolvedOxygen > 200|dissolvedOxygen == 0,
                                  NA, dissolvedOxygen),
         Temperature = ifelse(Temperature < 1, NA, Temperature),
         salinity = ifelse(salinity > 100|salinity < 5, NA, salinity),
         conductivity = ifelse(conductivity > 70|conductivity == 0, NA, conductivity))


save(fish, file = file.path(data.dir,"fish_processed.rdata"))

# Merge data from VIMS aerial surveys and UVA synoptic surveys-------------------

#load fish and SAV data
load(file.path(data.dir, "sav_cover.rdata"))
sav %<>% st_transform("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")

# This data set is drawn from fish seining data sheets. It refers to 
# notes determining whether or not there was eelgrass present during seining.
qual_class <- read_xlsx(file.path(data.dir, "seagrass_comments.xlsx"))

# Convert fish data.frame to sf object
fish_sf1 <- fish %>% 
  st_as_sf(., coords = c("longitude","latitude"), 
           crs = "+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs") %>% 
  mutate(year = year(sampleDate))

# Intersect seagrass and fish sampling geometries by year. 
# This merges VIMS seagrass classes and fish data
fish_sf <- NULL
for (i in c(2012:2018)){
  sav_filt <- sav %>% filter(year == i) 
  fish_filt <- fish_sf1 %>% filter(year == i)
  
  fish_sav_int <- st_join(fish_filt,  st_make_valid(sav_filt))
  assign('fish_sf', rbind(fish_sav_int, fish_sf))
}

fish_sf %<>% 
  dplyr::select(-year.y, -BEDLABEL, -ZONE_, -Surveyed) %>% 
  dplyr::rename(year = year.x, vimsSeagrassDens = DENSITY) %>% 
  mutate(vimsSeagrassDens = ifelse(str_detect(site, "HI") & year == 2013,
                                   NA, vimsSeagrassDens))

# Merge UVA Synoptic survey data --------------------------------------------------
load(file.path(data.dir, "syn_seagrass.rdata"))
synoptic_seagrass %<>% 
  dplyr::select(-Latitude, -Longitude)

# Join the data frames for further analysis
fish_sf2 <- fish_sf %>% 
  dream::sfc_as_cols(., names = c("Longitude","Latitude")) %>% 
  st_set_geometry(NULL) %>%
  left_join(.,synoptic_seagrass, by = c("year","season","site")) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"))  

# There are some instances where the VIMS data indicates no seagrass present, 
# but the UVA data include shoot densities. In these instances,
# the UVA data is defaulted to (i.e. seagrass is considered present)
sg_key <- fish_sf2 %>% 
  dplyr::select(year, season, site, vimsSeagrassDens, uvaShootDens, uvaPlotAge,
                canopy, biomass, ag_bio, bg_bio) %>% 
  st_set_geometry(NULL) %>% 
  distinct() %>% 
  left_join(.,qual_class) %>% 
  # dplyr::select(year, season, site, vimsSeagrassDens, uvaShootDens, qual_class, comment) %>% 
  dplyr::select(-comment) %>% 
  dplyr::rename(insideSeagrass = qual_class) %>% 
  group_by(site, season) %>% 
  mutate(plot_age = ifelse(str_detect(site, "SB") & !str_detect(site, "Bare"),
                           year - 2001,
                           year - 2006),
         plot_age = ifelse(str_detect(site, "Bare"), NA, plot_age))

#join in seagrass AND distance from sample sites to deep areas of the VCR
fish_sf %<>% 
  left_join(.,sg_key, by = c("year","season","site"))

st_crs(fish_sf) <- st_crs(sav)

# convert sample dates and time to datetime, and merge with Kiptopeke water levels and nearshore
# heatwave day data
fish_sf %<>% 
  mutate(sampleTime = format(strptime(substr(as.POSIXct(sprintf("%04.0f", sampleTime), 
                                                        format="%H%M"), 12, 16), '%H:%M'), '%I:%M %p'),
         sampleDate = as_date(sampleDate),
         datetime =  round_date(ymd_hm(paste(sampleDate, sampleTime)), "6 minutes")) %>%
  mutate(month = month(sampleDate)) 

save(fish_sf, fish_length_df,
     file = file.path(data.dir, "fish_processed_spatial.rdata"))

