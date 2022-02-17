library(tidyverse)
library(mgcv)
library(lavaan)
library(DHARMa)
library(MASS)
library(sf)
library(dream)
library(tsibble)
library(ggeffects)
library(forecast)
library(glmmTMB)
library(piecewiseSEM)
library(lme4)
library(lubridate)

plot_pred <- function(cov = NULL, resp = NULL,
                      mod = NULL, mod_df){
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
    labs(x = cov,
         y = resp)
}

# load data----
load(here::here("data/fish_processed_spatial.rdata"))
load(here::here("data/fish_processed_aggregated.rdata"))
source(here::here("R/multigroup.R"))
mod_df <- specs_ran2 %>% 
  mutate(ft = as.numeric(scale(ft)),
         rt2 = rt,
         rt = as.numeric(scale(rt)),
         d2 = depth,
         depth = as.numeric(scale(depth)),
         canopy = log(canopy),
         meadow = factor(ifelse(str_detect(site,
                                           "HI"),"HI","SB")),
         tsr = ifelse(str_detect(site, "SB"),
                      year - 2001,
                      year - 2006),
         sg_pres = ifelse(insideSeagrass == "Seagrass",1,0),
         season = factor(season),
         year = factor(year),
         meadow = factor(meadow),
         site2 = factor(ifelse(str_detect(site, "SB_Bare"),
                               "SB-bare", 
                               ifelse(str_detect(site, "HI_Bare"),
                                      "HI",
                                      as.character(meadow)))),
         site = factor(site),
         insideSeagrass = factor(insideSeagrass)) %>% 
  as.data.frame() 

# submodel 1----
# m_lm_meadow <-
#   glmer(sg_pres ~ (tsr + ft + rt + depth)*meadow +
#           (1|site) + (1|year),
#         family = "binomial",
#         data = mod_df,
#         control=glmerControl(optimizer="bobyqa",
#                              optCtrl=list(maxfun=1000000)))
# s <- simulateResiduals(m_lm_meadow, n = 1000);plot(s)
# 
# m_lm_season <-
#   glmer(sg_pres ~ (tsr + ft + rt + depth)*season +
#           (1|site) + (1|year),
#         family = "binomial",
#         data = mod_df,
#         control=glmerControl(optimizer="bobyqa",
#                              optCtrl=list(maxfun=1000000)))
# summary(m_lm_season)
# s <- simulateResiduals(m_lm_season, n = 1000);plot(s)

# submodel 1----
m_lm <- 
  glmer(sg_pres ~ ft + rt + depth + 
          (1|site) + (1|year),
        family = "binomial",
        data = mod_df,
        control=glmerControl(optimizer="bobyqa", 
                             optCtrl=list(maxfun=1000000)))

# m_lm2 <- glmmTMB(sg_pres ~ ft + rt + depth + 
#                    (1|site) + (1|year),
#                  family = "binomial",
#                  data = mod_df)

# s <- simulateResiduals(m_lm, n = 1000);plot(s)

# submodel 2----
# m_glmm_season <- glmer(richness ~ (sg_pres + rt)*season + 
#                          (1|year) + (1|site),
#                 family = "poisson",
#                 data = mod_df)
# s <- simulateResiduals(m_glmm_season, n = 1000);plot(s)
# 
# m_glmm_meadow <- glmer(richness ~ (sg_pres + rt)*meadow +
#                          (1|year) + (1|site),
#                 family = "poisson",
#                 data = mod_df)
# s <- simulateResiduals(m_glmm_meadow, n = 1000);plot(s)

m_glmm <- glmer(richness ~ sg_pres + rt + depth + (1|year),
                family = "poisson",
               data = mod_df)
performance::check_collinearity(m_glmm)
# s <- simulateResiduals(m_glmm, n = 1000);plot(s)

# m_psem_season <- psem(
#   m_lm,
#   m_glmm,
#   data = mod_df
# )
# multigroup(m_psem_season, group = "season")

# summary(m_glmm2)
# P = 0, not good fitting model
m_psem_rich <- psem(
  m_lm,
  m_glmm,
  data = mod_df
)
summary(m_psem_rich)
multigroup(m_psem_rich, group = "meadow", model_sim = F)
multigroup(m_psem_rich, group = "insideSeagrass", model_sim = F)
