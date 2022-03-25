library(readxl)
library(tidyverse)
library(lubridate)
library(heatwaveR)
library(tsibble)
library(magrittr)
library(sf)

# Oyster water temps
load(here::here("data/oyster_water_temps.rdata"))

# Get Bay heatwave data
wch_hw <- read_csv(here::here("data/tmsWACH.csv")) %>% 
  dplyr::select(-1)

# compare Wachapreague and Oyster data
wch_hw %>% 
  dplyr::rename(wch_temp = temp) %>% 
  right_join(.,may_2016, c("t" = "day")) %>% 
  as_tsibble() %>% 
  fill_gaps() %>% 
  ggplot() +
    geom_line(aes(y = temp, x = t))+
    geom_line(aes(y = wch_temp, x = t), color = "red")

comp_df <- 
  wch_hw %>% 
  dplyr::rename(wch_temp = temp) %>% 
  right_join(.,may_2016, c("t" = "day")) %>% 
  na.exclude()

# r = 0.99
cor(comp_df$wch_temp, comp_df$temp)

# bind in missing data from available data at Oyster
wch_hw2 <- 
  wch_hw %>% 
  dplyr::rename(wch_temp = temp) %>% 
  left_join(.,may_2016, by = c("t" = "day")) %>%
  mutate(wch_temp = ifelse(is.na(wch_temp) &
                             month(t) == 5 &
                             year(t) == 2016,
                           temp, wch_temp)) %>% 
  dplyr::select(t, temp = wch_temp, seas, thresh, X5Ev)

# Get fish data
load(here::here("data/fish_processed_spatial.rdata"))

sample_dates <- 
  fish_sf %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(sampleDate) %>% 
  distinct() %>% 
  pull()

sst_anom <- NULL
for (i in 1:length(sample_dates)){
  
  #DD index
  prev_month_anom <-
    wch_hw2 %>% 
      mutate(anom = temp - seas) %>% 
      filter(between(t, sample_dates[i] %m-% months(1), 
                     sample_dates[i])) %>% 
      pull(anom) %>% 
      mean()
  
  assign('sst_anom', rbind(sst_anom, 
                                tibble(prev_month_anom,
                                       sampleDate = sample_dates[i])))
}

save(sst_anom, file = here::here("data/sst_anom.rdata"))
