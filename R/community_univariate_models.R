library(tidyverse)
library(mgcv)
library(patchwork)
library(gratia)
library(DHARMa)
library(magrittr)
library(colorspace)
library(glmmTMB)
library(ggeffects)

load(here::here("data/fish_processed_aggregated.rdata"))

# process data to exclude years without VIMS aerial imagery
mod_df <- specs_ran %>% 
  mutate(meadow_month = factor(paste(month, meadow)),
         month = factor(case_when(month == 5 ~ "May",
                           month == 6 ~ "June",
                           month == 9 ~ "Sept.",
                           month == 10 ~ "Oct."), levels = c("May","June","Sept.","Oct.")),
         site2 = factor(ifelse(str_detect(site, "SB_Bare"),
                               "SB-bare", 
                               ifelse(str_detect(site, "HI_Bare"),
                                      "HI-Bare",
                                      as.character(meadow)))),
         hill_shannon = hill_shannon - 1,
         season = factor(season))

# ANOVA models----
# What's the difference in CPUE between inside and outside the meadow?----
sg_mod_abund <- 
  glmmTMB(abundance ~  insideSeagrass + 
            (1|site) + (1|year),
          family="nbinom2",
          data = mod_df)
# summary(sg_mod_abund)
# plot(ggpredict(sg_mod_abund))
# ggpredict(sg_mod_abund) %>% as_tibble()
# car::Anova(sg_mod_abund)
# dream::validate(list(sg_mod_abund))
# s <- DHARMa::simulateResiduals(sg_mod_abund, n = 1000);plot(s)

# What's the difference in Hill-Shannon diversity between inside and outside the meadow?----
# Note: Site random intercept dropped from model because of near-zero variance explained.
sg_mod_div <- 
  glmmTMB(hill_shannon ~ insideSeagrass +
            (1|site) + (1|year),
          ziformula = ~.,
          family = ziGamma(link = "log"),
          data = mod_df)

# ggpredict(sg_mod_div, type = "fe.zi") %>% as_tibble()
# t <- ggpredict(sg_mod_div, type = "fe.zi") %>% as_tibble()
# ggplot(mod_df) +
#   geom_point(aes(x = insideSeagrass, y = hill_shannon)) +
#   geom_point(data = t$insideSeagrass,
#              aes(x = x,
#                  y = predicted), color = "orange")
# 
# ggplot(mod_df) +
#   geom_histogram(aes(fill = insideSeagrass, x = hill_shannon)) +
#   geom_vline(data = t$insideSeagrass,
#              aes(xintercept = predicted, color = x))
# car::Anova(sg_mod_div)
# dream::validate(list(sg_mod_div))
# s <- DHARMa::simulateResiduals(sg_mod_div, n = 1000);plot(s)

# What's the difference in sample richness between inside and outside of the meadow?----
sg_mod_rich <- 
  glmmTMB(richness ~ insideSeagrass +
            (1|site) + (1|year),
          family = "poisson",
          data = mod_df)
# summary(sg_mod_rich)
# ggpredict(sg_mod_rich) %>% as_tibble()
# plot(ggpredict(sg_mod_rich))
# car::Anova(sg_mod_rich)
# dream::validate(list(sg_mod_rich))
# s <- DHARMa::simulateResiduals(sg_mod_rich, n = 1000);plot(s)


# Does CPUE vary seasonally?
season_mod_abund <- 
  glmmTMB(abundance ~ season + 
            (1|site) + (1|year),
          family="nbinom2",
          data = mod_df)
# summary(season_mod_abund)
# plot(ggpredict(season_mod_abund))
# ggpredict(season_mod_abund) %>% as_tibble()
# car::Anova(season_mod_abund)
# em_mon_mod <- emmeans::emmeans(season_mod_abund, "season")
# pairs(em_mon_mod)
# dream::validate(list(season_mod_abund))
# s <- DHARMa::simulateResiduals(sg_mod_abund, n = 1000);plot(s)


season_mod_div <- 
  glmmTMB(hill_shannon ~ season +
            (1|site) + (1|year),
          ziformula = ~.,
          family =  ziGamma(link = "log"),
          data = mod_df)
# summary(season_mod_div)
# plot(ggpredict(season_mod_div))
# car::Anova(season_mod_div)
# em_div_mod <- emmeans::emmeans(season_mod_div, "season")
# pairs(em_div_mod)
# dream::validate(list(season_mod_div))

season_mod_rich <- 
  glmmTMB(richness ~ season +
            (1|site) + (1|year),
          family = "poisson",
          data = mod_df)
# summary(season_mod_rich)
# plot(ggpredict(season_mod_rich))
# car::Anova(season_mod_rich)
# dream::validate(list(season_mod_rich))

meadow_mod_rich <- 
  glmmTMB(richness ~ meadow +
            (1|site) + (1|year),
          family = "poisson",
          data = mod_df)
# summary(meadow_mod_rich)
# plot(ggpredict(meadow_mod_rich))
# car::Anova(meadow_mod_rich)
# dream::validate(list(meadow_mod_rich))

meadow_mod_abund <- 
  glmmTMB(abundance ~ meadow +
            (1|site) + (1|year),
          family = "nbinom2",
          data = mod_df)
summary(meadow_mod_abund)

sg_mod_dens <- 
  glmmTMB(shoot_dens ~ meadow +
            (1|site) + (1|year),
          family = Gamma(link = "log"),
          data = mod_df)
summary(sg_mod_dens)
# plot(ggpredict(sg_mod_dens))
car::Anova(sg_mod_dens)
dream::validate(list(sg_mod_dens))
s <- DHARMa::simulateResiduals(sg_mod_dens); plot(s)

summary(meadow_mod_abund)

# ANOVA visualizations----
ilink <- family(sg_mod_abund)$linkinv
sg_mod_abund_pred <- 
  cbind(mod_df, 
        predict(sg_mod_abund,
          mod_df %>% dplyr::select(abundance, year, site, month, insideSeagrass),
          se.fit=TRUE,
          type = "link")) %>%
  mutate(sg_abund_upr = ilink(fit + (2 * se.fit)),
         sg_abund_lwr = ilink(fit - (2 * se.fit)),
         sg_abund_fit = ilink(fit))

ilink <- family(sg_mod_div)$linkinv
sg_mod_div_pred <- 
  cbind(mod_df, 
        predict(sg_mod_div,
                mod_df,
                se.fit=TRUE,
                type = "link")) %>%
  mutate(sg_div_upr = ilink(fit + (2 * se.fit)),
         sg_div_lwr = ilink(fit - (2 * se.fit)),
         sg_div_fit = ilink(fit))

ilink <- family(sg_mod_rich)$linkinv
sg_mod_rich_pred <- 
  cbind(mod_df, 
        predict(sg_mod_rich,
                mod_df,
                se.fit=TRUE,
                type = "link")) %>%
  mutate(sg_rich_upr = ilink(fit + (2 * se.fit)),
         sg_rich_lwr = ilink(fit - (2 * se.fit)),
         sg_rich_fit = ilink(fit))


ilink <- family(season_mod_abund)$linkinv
season_mod_abund_pred <- 
  cbind(mod_df, 
        predict(season_mod_abund,
                mod_df,
                se.fit=TRUE,
                type = "link")) %>%
  mutate(season_abund_upr = ilink(fit + (2 * se.fit)),
         season_abund_lwr = ilink(fit - (2 * se.fit)),
         season_abund_fit = ilink(fit))

ilink <- family(season_mod_div)$linkinv
season_mod_div_pred <- 
  cbind(mod_df, 
        predict(season_mod_div,
                mod_df,
                se.fit=TRUE,
                type = "link")) %>%
  mutate(season_div_upr = ilink(fit + (2 * se.fit)),
         season_div_lwr = ilink(fit - (2 * se.fit)),
         season_div_fit = ilink(fit))

ilink <- family(season_mod_rich)$linkinv
season_mod_rich_pred <- 
  cbind(mod_df, 
        predict(season_mod_rich,
                mod_df,
                se.fit=TRUE,
                type = "link")) %>%
  mutate(season_rich_upr = ilink(fit + (2 * se.fit)),
         season_rich_lwr = ilink(fit - (2 * se.fit)),
         season_rich_fit = ilink(fit))


# Inside vs outside CPUE 
# sig_df1 <- tibble(insideSeagrass = "Seagrass",
#                  abundance = 105,
#                  label = "**")
(sg_abund_figure <- 
  ggplot(data = sg_mod_abund_pred) +
    geom_boxplot(aes(x = insideSeagrass, y = abundance,
                    fill = insideSeagrass)) +
  labs(y = "Catch per seine tow") +
  guides(color = "none",
         fill = "none") +
  ggsci::scale_color_d3() +
  ggsci::scale_fill_d3() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.title.y = element_text(size = 11)))

# Inside vs outside diversity
# sig_df2 <- tibble(insideSeagrass = "Seagrass",
#                  hill_shannon = 5.75 ,
#                  label = "***")
(sg_div_figure <- 
    ggplot(sg_mod_div_pred) +
    geom_violin(aes(x = insideSeagrass, y = hill_shannon + 1),
                size = 0.75) +
    geom_boxplot(aes(x = insideSeagrass, y = hill_shannon + 1,
                     color = insideSeagrass, fill = insideSeagrass),
                 width=0.1, alpha = 0.25) +
    labs(y = "Hill-Shannon Diversity") +
    guides(color = "none",
           fill = "none") +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 11)))

# Inside vs outside diversity
# sig_df2 <- tibble(insideSeagrass = "Seagrass",
#                  hill_shannon = 5.75 ,
#                  label = "***")
(sg_rich_figure <- 
    ggplot(sg_mod_rich_pred) +
    geom_boxplot(aes(x = insideSeagrass, y = richness,
                    fill = insideSeagrass),
                width = 0.4) +
    labs(y = "Sample richness") +
    guides(color = "none",
           fill = "none") +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 11)))

# seasonal CPUE 
(mon_abund_figure <- 
    ggplot(season_mod_abund_pred) +
    geom_boxplot(aes(x = month, fill = insideSeagrass, 
                     color = insideSeagrass, 
                     y = abundance),
                alpha = 0.25) +
    labs(y = "CPUE") +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    guides(color = "none",
           fill = "none") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 11)))

# Inside vs outside diversity
(mon_div_figure <- 
    ggplot(season_mod_div_pred) +
    geom_boxplot(aes(x = month, fill = insideSeagrass, 
                     color = insideSeagrass, 
                     y = hill_shannon + 1),
                 alpha = 0.25) +
    labs(y = "Hill-Shannon Diversity") +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    labs(color = "Sampling habitat",
           fill = "Sampling habitat") +
    theme_bw() +
    ylim(0, 6) +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 11),
          legend.position = "bottom"))

(mon_rich_figure <- 
    ggplot(season_mod_rich_pred) +
    geom_boxplot(aes(x = month, fill = insideSeagrass, 
                     color = insideSeagrass, 
                     y = richness),
                 alpha = 0.25) +
    labs(y = "Sample Richness") +
    ggsci::scale_color_d3() +
    ggsci::scale_fill_d3() +
    labs(color = "Sampling habitat",
         fill = "Sampling habitat") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = 10),
          axis.title.y = element_text(size = 11),
          legend.position = "bottom"))

(anova_fig1 <- 
  # sg_abund_figure + sg_div_figure + 
    mon_abund_figure + mon_div_figure + plot_annotation(tag_levels = "A") + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom"))

(anova_fig2 <- 
    sg_abund_figure + sg_rich_figure + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom"))

ggsave(anova_fig2, filename = here::here("figures/inside_outside2.png"),
       width = 8, height= 4)

save(anova_fig1,
     file = here::here("data/sg_cover_mods_figs.rdata"))
