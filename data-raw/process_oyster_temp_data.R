library(dplyr)
library(lubridate)
library(tsibble)
library(stringr)
library(magrittr)
library(xts)
library(tidyr)

# Porter, J., D. Krovetz, J. Spitler, J. Spitler, T. Williams and K. Overman. 2019. Tide Data for Hog Island (1991-), 
# Redbank (1992-), Oyster (2007-) . Virginia Coast Reserve Long-Term Ecological Research Project Data Publication
# knb-lter-vcr.61.33 (http://www.vcrlter.virginia.educgi-bin/showDataset.cgi?docid=knb-lter-vcr.61.33).

print(paste("Time of data pull is", Sys.time()))

# Read in data from VCR database
fname <- "http://www.vcrlter.virginia.edu/data/metdata/metgraphs/tidedata/VCRTide.csv"
# fname <- "http://www.vcrlter.virginia.edu/data/metdata/metgraphs/csv/hourly/todayTide.csv"
infile1 <- data.table::fread(fname, col.names = c("station",
                                                  "date",
                                                  "time",
                                                  "relative_tide_level",
                                                  "water_temperature"),
                             quote = '"',
                             select = c(1:5)) %>% 
  tibble::tibble() %>% 
  mutate(time = as.numeric(time),
         relative_tide_level = as.numeric(relative_tide_level))

# infile2 <- readr::read_csv(fname, col_names = c("station",
#                                      "date",
#                                      "time",
#                                      "relative_tide_level",
#                                      "water_temperature",
#                                      "barometric_pressure"),
#                 quote = '"')

# Process for use in package format
tides_new_df <-
  infile1 %>% 
    
    # select data from the Oyster station only
    dplyr::filter(stringr::str_detect(station, "OYST")) %>% 
    dplyr::select(-station) %>% 
  
    # convert missing values to NA
    mutate_all(function(x)ifelse(x == ".", NA, x)) %>% 
  
    # process the datetimes
    dplyr::mutate(date = as.Date(date, format = "%d%b%Y"),
           time = format(strptime(substr(as.POSIXct(sprintf("%04.0f", as.numeric(time)),
                                                    format="%H%M", tz = "America/New_York"), 12, 16), 
                                  '%H:%M'), '%I:%M %p'),
           datetime = as.POSIXct(paste(date, time), 
                                 format="%Y-%m-%d %I:%M %p", tz = "America/New_York"),
           water_temperature = as.numeric(water_temperature)) %>% 
    dplyr::filter(date <= Sys.Date()) %>% 
    dplyr::select(-date, -time) %>%
    dplyr::filter(!is.na(datetime)) %>%
  
    # Convert to tsibble to fill missing datetimes
    dplyr::filter(!duplicated(datetime)) %>% 
    tsibble(index = datetime) %>% 
    fill_gaps(.full  = TRUE) %>% 
    as_tibble() %>%
  
    # get the hourly mean
    group_by(y = year(datetime),
             m = month(datetime),
             d = day(datetime),
             h = hour(datetime)) %>% 
    dplyr::summarise(relative_tide_level = mean(relative_tide_level, na.rm = T),
                     water_temperature = mean(water_temperature, na.rm = T)) %>% 
  
    # convert back to datetime
    dplyr::mutate(datetime = lubridate::ymd_h(paste(y, m, d, h, sep = "-"), tz = "America/New_York")) %>%
    
    # select variables of interest
    dplyr::ungroup() %>% 
    dplyr::select(datetime, relative_tide_level, water_temperature)
 
may_2016 <- 
  tides_new_df %>% 
  filter(year(datetime) == 2016,
         month(datetime) %in% 4:8) %>% 
  group_by(day = date(datetime)) %>% 
  dplyr::summarise(temp = mean(water_temperature))
save(may_2016, file = here::here("data/oyster_water_temps.rdata"))
