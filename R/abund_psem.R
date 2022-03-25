library(tidyverse)
library(lavaan)
library(DHARMa)
library(MASS)
library(sf)
library(semEff)
library(dream)
library(tsibble)
library(ggeffects)
library(forecast)
library(glmmTMB)
library(piecewiseSEM)
library(lme4)
library(lubridate)

plot_pred <- function(cov = NULL, resp = NULL,
                      mod = NULL, mod_df, binom = T,
                      lab_x = NULL, lab_y = NULL,
                      scale = 100,
                      b = 20){
  if (is.null(all(c(lab_x,lab_y)))) {
    lab_x <- cov
    lab_y <- resp
  }
  if (binom){
    p1 <- ggpredict(mod, 
                    terms = paste0(cov,"[all]")) %>% 
      as_tibble() 
    ggplot() +
      geom_ribbon(data = p1, aes(x = x,
                                 ymin = conf.low,
                                 ymax = conf.high),
                  alpha = 0.5) +
      geom_histogram(data = mod_df %>% 
                     filter(get(resp) == 1),
                   aes(x = get(cov), 
                       y = -1*stat(count/scale)), 
                   bins = b, na.rm = TRUE, 
                   position = position_nudge(y = 1)) + 
      geom_histogram(data = mod_df %>% 
                       filter(get(resp) == 0),
                     aes(x = get(cov), y = stat(count/scale)),
                     bins = b, na.rm = TRUE) +
      geom_line(data = p1, aes(x = x, y = predicted, group = group)) +
      labs(x = lab_x,
           y = lab_y) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      theme(axis.title.x = element_text(size = 7))
  } else {
    p1 <- ggpredict(mod, 
                    terms = paste0(cov,"[all]")) %>% 
      as_tibble() 
    ggplot() +
      geom_ribbon(data = p1, aes(x = x,
                                 ymin = conf.low,
                                 ymax = conf.high),
                  alpha = 0.5) +
      geom_point(data = mod_df, aes(y = get(resp), 
                                    x = get(cov))) +
      geom_line(data = p1, aes(x = x, y = predicted, group = group)) +
      labs(x = lab_x,
           y = lab_y) 
  }
  
}

# load data----
load(here::here("data/fish_processed_spatial.rdata"))
load(here::here("data/fish_processed_aggregated.rdata"))
source(here::here("R/multigroup.R"))
mod_df <- specs_ran2 %>% 
  mutate(ft = as.numeric(scale(ft)),
         rt = as.numeric(scale(rt)),
         depth = as.numeric(scale(depth)),
         canopy = log(canopy),
         meadow = factor(ifelse(str_detect(site,
                                           "HI"),"HI","SB")),
         tsr = ifelse(str_detect(site, "SB"),
                      year - 2001,
                      year - 2006),
         sg_pres = ifelse(insideSeagrass == "Seagrass",1,0),
         ag_bio = log(ag_bio * uvaShootDens),
         season = factor(season),
         year = factor(year),
         meadow = factor(meadow),
         site2 = factor(ifelse(str_detect(site, "SB_Bare"),
                               "SB-bare", 
                               ifelse(str_detect(site, "HI_Bare"),
                                      "HI-bare",
                                      as.character(meadow)))),
         ys = factor(paste(year, site2)),
         site = factor(site)) %>% 
  as.data.frame() 

# ggplot(mod_df %>% filter(meadow == "HI")) +
#   geom_point(aes(x = longitude, y = latitude, color = site)) +
#   ggsci::scale_color_d3()
# 
# ggplot(mod_df %>% filter(meadow == "HI")) +
#   geom_point(aes(x = sampleDate, y = insideSeagrass, color = year))  +
#   facet_wrap(~site) +
#   ggsci::scale_color_d3()
# 
# ggplot() +
#   geom_point(data = mod_df %>% filter(meadow == "HI",
#                                       !str_detect(site, "Bare|72")),
#              aes(x = sampleDate, y = abundance, group = site, color = insideSeagrass))  +
#   geom_line(data = mod_df %>% filter(meadow == "HI",
#                                       !str_detect(site, "Bare|72")),
#              aes(x = sampleDate, y = abundance, group = site, color = insideSeagrass)) 
# 
# ggplot(mod_df) +
#   geom_point(aes(x = depth, y  = rt))
# 
# ggplot(mod_df) +
#   geom_point(aes(x = sg_dist, y = sg_pres, color = meadow),
#              position = position_jitter(height = 0, width = 0.1))
# 
# ggplot(mod_df) +
#   geom_point(aes(y = sg_pres, x = ft)) +
#   facet_wrap(~year)

# submodel 1----
# m_lm_meadow <- 
#   glmer(sg_pres ~ (tsr + ft + rt + depth)*meadow + 
#           (1|site) + (1|year),
#         family = "binomial",
#         data = mod_df,
#         control=glmerControl(optimizer="bobyqa", 
#                              optCtrl=list(maxfun=1000000)))
# s <- simulateResiduals(m_lm_meadow, n = 1000);plot(s)

# m_lm_season <- 
#   glmer(sg_pres ~ (tsr + ft + rt + depth)*season + 
#           (1|site) + (1|year),
#         family = "binomial",
#         data = mod_df,
#         control=glmerControl(optimizer="bobyqa", 
#                              optCtrl=list(maxfun=1000000)))
# s <- simulateResiduals(m_lm_season, n = 1000);plot(s)

m_lm <- 
  glmer(sg_pres ~ ft + rt + depth + prev_month_anom +
          (1|site) + (1|year),
        family = "binomial",
        data = mod_df,
        control=glmerControl(optimizer="bobyqa", 
                             optCtrl=list(maxfun=1000000)))

# submodel 2----
# s <- simulateResiduals(m_glmm_seas, n = 1000);plot(s)
# 
# m_glmm_meadow <- glmer.nb(abundance ~ (sg_pres + rt)*meadow + 
#                           (1|site) + (1|year),
#                         data = mod_df)
# s <- simulateResiduals(m_glmm_meadow, n = 1000);plot(s)
m_glmm <- glmer.nb(abundance ~ sg_pres + rt + prev_month_anom +
                            (1|year),
                          data = mod_df)
summary(m_glmm)
m_psem_abund <- psem(
  m_lm,
  m_glmm,
  data = mod_df
)
summary(m_psem_abund)
multigroup(m_psem_abund, group = "season", model_sim = F)

# Visualizing submodels----
(depth_fig <- 
  plot_pred(cov = "depth", resp = "sg_pres", mod = m_lm,
            lab_y = "Probability of\neelgrass presence",
            lab_x = "Bathymetry (normalized NAVD88 m)",
            mod_df,
            b = 30) +
  theme_bw())
(rt_fig <- 
  plot_pred(cov = "rt", resp = "sg_pres", mod = m_lm,
            lab_y = "Probability of\neelgrass presence",
            lab_x = "Residence time (normalized hours)",
            mod_df,
            b = 30) +
  theme_bw())

eelgrass_sem_fig <- 
  depth_fig + rt_fig +
  plot_annotation(tag_levels = "A")

ggplot(mod_df) +
  geom_point(aes(y = abundance, x = prev_month_anom, color = month))+
  geom_smooth(aes(y = abundance, x = prev_month_anom, color = month),
              method = "glm.nb") +
  facet_wrap(~month) +
  labs(x = "SST anomaly in previous month (°C)") +
  theme_bw()

ggplot(mod_df) +
  geom_point(aes(y = abundance, x = prev_month_anom, color = month))+
  geom_smooth(aes(y = abundance, x = prev_month_anom, color = month),
              method = "glm.nb") +
  facet_wrap(~season) +
  labs(x = "SST anomaly in previous month (°C)") +
  theme_bw()

save(eelgrass_sem_fig,
       file = here::here("figures/eelgrass_sem_fig.rdata"))

 