library(tidyverse)
library(sf)
library(mgcv)
library(lubridate)
data.dir <- here::here("data")

load(file.path(data.dir, "fish_processed_aggregated.rdata"))
load(file.path(data.dir, "fish_processed_spatial.rdata"))

fish_length_df %>% 
  filter(!is.na(Length)) %>% 
  mutate(month = factor(month(sampleDate), levels = c(5,6,9,10)),
         year = factor(year(sampleDate))) %>% 
  ggplot() +
  geom_boxplot(aes(x = month, y = Length, color = year))

t <- 
  fish_length_df %>% 
  filter(!is.na(Length),
         sciName %in% c("Anchoa spp.",
    "Bairdiella chrysoura",
    "Leiostomus xanthurus",
    "Menidia menidia",
    "Lagodon rhomboides",
    "Syngnathus spp.")) %>% 
    mutate(common = plyr::mapvalues(sciName,
                                    from = unique(.$sciName),
                                  to = c("Anchovies",
                                         "Silver perch",
                                         "Pinfish",
                                         "Spot",
                                         "Atl. Silverside",
                                         "Pipefish")))

lam <- 
  tibble(sciName = c("Anchoa spp.",
                     "Bairdiella chrysoura",
                     "Lagodon rhomboides",
                     "Leiostomus xanthurus",
                     "Menidia menidia",
                     "Syngnathus spp."),
         common = c("Anchovies",
                    "Silver perch",
                    "Pinfish",
                    "Spot",
                    "Atl. Silverside",
                    "Pipefish"),
         Length_min = c(40, 91, 110, 135, 91, 91),
         Length_max = c(45 ,95, NA, 158, 98, 125),
         Measure = c("FL", "SL", "SL", "SL", "TL","SL"))


# silverside -> fay
# spot -> johnson
# pinfish -> ohs
# newburger -> anchovies
# pipefish -> ripley and foran

ggplot(t) +
  geom_histogram(aes(x = Length, color = month(sampleDate))) +
  geom_rect(data = lam,
            aes(xmin=Length_min, 
                xmax=Length_max, ymin=0, ymax=Inf),
            alpha = 0.3) +
  geom_vline(data = lam %>% filter(sciName == "Lagodon rhomboides"),
             aes(xintercept = Length_min)) +
  facet_wrap(~common) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.001, 0.001)) +
  labs(x = "Total Length (mm)",
       y = "Frequency") +
  ecodata::theme_facet()


mod_df <- 
  fish_length_df %>% 
  filter(!is.na(Length)) %>% 
  left_join(.,specs_ran %>% 
              filter(!is.na(sg_dist)) %>% 
              dplyr::select(-sampleTime),
            by = c("site","latitude","longitude","sampleDate")) %>% 
  # filter(sciName %in% c("Anchoa spp.",
  #                       "Bairdiella chrysoura",
  #                       "Leiostomus xanthurus",
  #                       "Menidia menidia",
  #                       "Lagodon rhomboides",
  #                       "Syngnathus spp.")) %>%
  mutate(season = factor(ifelse(month %in% c(5,6), "Summer","Fall")),
         site2 = factor(ifelse(str_detect(site, "SB_Bare"),
                                "SB-bare", 
                                ifelse(str_detect(site, "HI_Bare"),
                                       "HI-Bare",
                                       as.character(meadow)))),
          site = factor(site),
          year = factor(year),
          sciName = factor(sciName),
          site_time = factor(paste(sampleDate, site)),
         ys = factor(paste(year, season))) %>% 
  group_by(year, season, ys, site_time,site, 
           longitude, latitude,sciName, mn_sg_dens) %>%
  dplyr::summarise(med_len = mean(Length, na.rm = T))

mod <- 
  gam(med_len ~ mn_sg_dens +
            # season +
            # s(mn_sg_dens,sciName, bs = "fs") + 
            # s(sciName, bs = "re") + 
            te(longitude, latitude, by = ys, bs = "gp"),
            # s(site2, bs = "re"),
            # (1 + mn_sg_dens|sciName) +
            # (1|sciName),
            # (1|year) + (1|site),
          family = Gamma(link = "log"),
        data = mod_df)

summary(mod)
s <- DHARMa::simulateResiduals(mod, n = 1000)
gam.check(mod)
plot(s)
gratia::appraise(mod)
dream::validate(list(mod))
gratia::draw(mod)
 
ggplot(mod_df) +
  geom_point(aes(x = mn_sg_dens, y = Length)) +
  facet_wrap(~sciName)

mod_df %>% 
  ggplot() +
    geom_point(aes(x = sg_cover_400, y = Length, group = site, color = site2)) +
    facet_wrap(~season) +
    geom_smooth(aes(x = sg_cover_400, y = Length),
                method = "glm")


fish_length_df %>% 
  filter(!is.na(Length)) %>% 
  mutate(month = factor(month(sampleDate), levels = c(5,6,9,10)),
         year = factor(year(sampleDate))) %>% 
  ungroup() %>% 
  group_by(year, month, sampleDate) %>% 
  summarise(m = mean(Length)) %>% 
  ggplot() +
  geom_point(aes(x = sampleDate, y = m, color = month))

fish_length_df %>% 
  filter(!is.na(Length)) %>% 
  mutate(month = factor(month(sampleDate), levels = c(5,6,9,10)),
         year = factor(year(sampleDate))) %>% 
  group_by(sciName) %>% 
  dplyr::summarise(mlen = mean(Length),
                   range_low = range(Length)[1],
                   range_high = range(Length)[2]) %>% View

syng_lens <- 
  fish_length_df %>%  
  filter(!is.na(Length)) %>% 
  mutate(month = factor(month(sampleDate), levels = c(5,6,9,10)),
         year = factor(year(sampleDate))) %>% 
  filter(sciName == "Syngnathus spp.", month %in% c(5,6)) %>% 
  group_by(year = year(sampleDate), month, site, sampleDate) %>% 
  dplyr::summarise(length = median(Length, na.rm = T))

ggplot(syng_lens) +
  geom_histogram(aes(length))

