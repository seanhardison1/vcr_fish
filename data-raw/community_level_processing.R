library(vegan)
library(nlme)
library(tidyverse)
library(sf)
library(patchwork)
library(MeanRarity)
library(lubridate)
library(magrittr)
library(DHARMa)
library(ggeffects)
library(broom)
library(broom.mixed)
library(ggeffects)
library(MASS)
library(glmmTMB)
library(geosphere)
library(lme4)
library(dream)
library(GGally)
library(landscapemetrics)

helper.dir <- here::here("R/helper")
r.dir <- here::here("R")
data.dir <- here::here("data")

load(file.path(data.dir, "fish_processed_spatial.rdata"))
load(file.path(data.dir, "biodiversity_table.rdata"))
div_summ %<>%  dplyr::select(-speciesName_ital, -nFish)
load(file.path(data.dir, "phys_vars.rdata"))

#remove geometry for ease of use
fish <- fish_sf %>% 
  sfc_as_cols(names = c("longitude","latitude")) %>% 
  st_set_geometry(NULL) %>% 
  mutate( meadow = str_extract(site, "(?:(?!_).)*")) 

# Disaggregated data for use with GLMMs-------------------------------------------------------------
specs_ran <- 
  fish %>% 
  group_by(year, plot_age, season, meadow, 
           insideSeagrass, 
           sampleDate, sampleTime, 
           site, 
           uvaShootDens,
           latitude, longitude) %>% 
  dplyr::summarise(richness = sum(nFish != 0),
                   hill_shannon = rarity(nFish, 0),
                   hill_simpson = rarity(nFish, -1),
                   abundance = sum(nFish),
                   pipefish_catch = sum(nFish[sciName == "Syngnathus spp."]),
                   shannon = diversity(nFish),
                   evenness = shannon/log(richness),
                   salinity = mean(salinity, na.rm = T),
                   DO = mean(dissolvedOxygen, na.rm = T),
                   conductivity = mean(conductivity, na.rm = T),
                   shoot_dens = mean(uvaShootDens, na.rm = T)) %>% 
  left_join(.,phys_vars, by = "site") %>% 
  ungroup() %>% 
  mutate_at(vars(season, meadow, 
                 insideSeagrass, site), as.factor) %>%  
  unite("yearseason",c("year","season"), remove = F) %>% 
  mutate(insideSeagrass = factor(ifelse(`insideSeagrass` == 0,
                                        "Unvegetated",
                                        "Seagrass")),
         doy = yday(sampleDate),
         month = factor(month(sampleDate)),
         month2 = factor(case_when(month == 5 ~ "May",
                           month == 6 ~ "June",
                           month == 9 ~ "Sept.",
                           month == 10 ~ "Oct.")),
         year = factor(year),
         season = factor(str_to_title(season), levels = c("Summer","Fall")),
         yearseas = factor(paste(year, season)),
         bimonthly = factor(ifelse(month == 5 & day(sampleDate) <= 15, 5,
                                   ifelse(month == 5 & day(sampleDate) > 15, 5.5,
                                          ifelse(month == 6 & day(sampleDate) <= 15, 6,
                                                 ifelse(month == 6 & day(sampleDate) > 15, 6.5,
                                                        ifelse(month == 9 & day(sampleDate) <= 15, 9,
                                                               ifelse(month == 9 & day(sampleDate) > 15, 9.5,
                                                                      ifelse(month == 10 & day(sampleDate) <= 15, 10,
                                                                             ifelse(month == 10 & day(sampleDate) > 15, 10.5,NA))))))))),
         site_month = factor(paste(month,
                                   ifelse(str_detect(site, "SB_Bare"),
                                          "SB-bare", 
                                          ifelse(str_detect(site, "HI_Bare"),
                                                 "HI-Bare",
                                                 as.character(meadow)))))) %>% 
  as.data.frame() 

# wide data
specs_wide <-
  fish %>%
  left_join(.,div_summ) %>% 
  dplyr::select(nFish, site, uvaShootDens,
                year, meadow, season, insideSeagrass, 
                sciName, common, sampleDate, sampleTime) %>% 
  group_by(year, season, meadow, 
           insideSeagrass, 
           sampleDate, sampleTime, uvaShootDens,
           site, common) %>% 
  dplyr::summarise(totFish = sum(nFish)) %>%
  ungroup() %>% 
  pivot_wider(names_from = common, values_from = totFish) %>% 
  mutate(n = rowSums(dplyr::select(.,`American halfbeak`:Tautog)),
         month = factor(month(sampleDate)),
         site_month = factor(paste(month,
                                   ifelse(str_detect(site, "SB_Bare"),
                                          "SB-bare", 
                                          ifelse(str_detect(site, "HI_Bare"),
                                                 "HI-Bare",
                                                 as.character(meadow))))),
         site2 = factor(ifelse(str_detect(site, "SB_Bare"),
                               "SB-bare", 
                               ifelse(str_detect(site, "HI_Bare"),
                                      "HI-Bare",
                                      as.character(meadow))))) %>% 
  {. ->> specs_wide2} %>% 
  ungroup() %>% 
  dplyr::select(-n)


specs_ran2 <- 
  fish %>% 
  dplyr::rename(vims_sg_dens = vimsSeagrassDens.x) %>% 
  group_by(year, plot_age, season, meadow, 
           insideSeagrass, 
           sampleDate, sampleTime, 
           site, 
           vims_sg_dens,
           uvaShootDens,
           ag_bio,
           canopy,
           latitude, longitude) %>% 
  dplyr::summarise(richness = sum(nFish != 0),
                   hill_shannon = rarity(nFish, 0),
                   hill_simpson = rarity(nFish, -1),
                   abundance = sum(nFish),
                   pipefish = sum(nFish[sciName == "Syngnathus spp."]),
                   silverside = sum(nFish[sciName == "Menidia menidia"]),
                   pinfish = sum(nFish[sciName == "Lagodon rhomboides"]),
                   perch = sum(nFish[sciName == "Bairdiella chrysoura"]),
                   spot = sum(nFish[sciName == "Leiostomus xanthurus"]),
                   shannon = diversity(nFish),
                   evenness = shannon/log(richness),
                   salinity = mean(salinity, na.rm = T),
                   DO = mean(dissolvedOxygen, na.rm = T),
                   conductivity = mean(conductivity, na.rm = T),
                   shoot_dens = mean(uvaShootDens, na.rm = T)) %>% 
  ungroup() %>% 
  mutate_at(vars(season, meadow, 
                 insideSeagrass, site), as.factor) %>%  
  unite("yearseason",c("year","season"), remove = F) %>% 
  mutate(insideSeagrass = factor(ifelse(`insideSeagrass` == 0,
                                        "Unvegetated",
                                        "Seagrass")),
         doy = yday(sampleDate),
         month = factor(month(sampleDate)),
         month2 = factor(case_when(month == 5 ~ "May",
                                   month == 6 ~ "June",
                                   month == 9 ~ "Sept.",
                                   month == 10 ~ "Oct.")),
         yearseas = factor(paste(year, season)),
         bimonthly = factor(ifelse(month == 5 & day(sampleDate) <= 15, 5,
                                   ifelse(month == 5 & day(sampleDate) > 15, 5.5,
                                          ifelse(month == 6 & day(sampleDate) <= 15, 6,
                                                 ifelse(month == 6 & day(sampleDate) > 15, 6.5,
                                                        ifelse(month == 9 & day(sampleDate) <= 15, 9,
                                                               ifelse(month == 9 & day(sampleDate) > 15, 9.5,
                                                                      ifelse(month == 10 & day(sampleDate) <= 15, 10,
                                                                             ifelse(month == 10 & day(sampleDate) > 15, 10.5,NA))))))))),
         site_month = factor(paste(month,
                                   ifelse(str_detect(site, "SB_Bare"),
                                          "SB-bare", 
                                          ifelse(str_detect(site, "HI_Bare"),
                                                 "HI-Bare",
                                                 as.character(meadow)))))) %>% 
  left_join(.,phys_vars)


save(specs_wide, specs_wide2, specs_ran2, specs_ran,
     file = here::here("data/fish_processed_aggregated.rdata"))
