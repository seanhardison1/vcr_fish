library(tidyverse)
library(glmmTMB)
library(ggtext)
library(sf)
library(dream)
library(ggeffects)
library(patchwork)
library(magrittr)

# load data----
load(here::here("data/fish_processed_spatial.rdata"))
load(here::here("data/fish_processed_aggregated.rdata"))

# process----
mod_df <- 
  specs_wide2 %>% 
  dplyr::select(`American halfbeak`:Tautog, insideSeagrass,
                sampleDate, site, month, season, year) %>% 
  gather(sciName, nFish, -site, -season, -year,
         -month, -sampleDate,-insideSeagrass) %>% 
  mutate(month = factor(case_when(month == 5 ~ "May",
                                  month == 6 ~ "June",
                                  month == 9 ~ "Sept.",
                                  month == 10 ~ "Oct."),
                        levels = c("May","June","Sept.","Oct.")),
         season = factor(str_to_title(season), levels = c("Summer","Fall"))) %>% 
  filter(year != 2011)

# model function----
fit_glmm <- function(mod_df, sciName, cov = "insideSeagrass"){
  if (cov == "insideSeagrass"){
    mod <- glmmTMB::glmmTMB(nFish ~ insideSeagrass + (1|site) + (1|year),
                            family = "nbinom2", 
                            data = mod_df)
    pdf <- ggpredict(mod) %>% 
      as_tibble() 
    pdf <- pdf$insideSeagrass
    outside <- pdf %>% filter(x == 0) %>% pull(predicted)
    inside <- pdf %>% filter(x == 1) %>% pull(predicted)
    print(paste0(sciName," ", round(inside/outside,2), "X more common inside meadow"))
  } else if (cov == "season"){
    mod <- glmmTMB::glmmTMB(nFish ~ season + (1|site) + (1|year),
                            family = "nbinom2", 
                            data = mod_df)
  }
  
  print(dream::validate(list(mod)))
  # s <- simulateResiduals(mod, n = 500); plot(s)
  
  return(mod)
}



# prediction function----
pred_func <- function(mod, spec, cov = "insideSeagrass"){
  
  pred <- 
    ggeffects::ggpredict(mod, terms = cov,
                         type = "fe") %>% 
    as_tibble()
  
  
  if (cov == "insideSeagrass"){
    pred %<>% 
      dplyr::rename(insideSeagrass = x, 
                    fit = predicted) %>% 
      mutate(sciName = spec,
             insideSeagrass = factor(ifelse(insideSeagrass == 1, "Seagrass", 
                                            "Unvegetated"),
                                     levels = c("Seagrass", "Unvegetated"))) 
  } else if (cov == "season"){
    pred %<>% 
      dplyr::rename(season = x, 
                    fit = predicted) %>% 
      mutate(sciName = spec)
  }
  
  return(pred)
}

# extract estimates----
fit_stat <- function(mod){
  broom.mixed::tidy(mod)
}

most_common <- c("Pipefish", 
                 "Atlantic silverside",
                 "Pinfish",
                 "Silver perch",
                 "Anchovies",
                 "Spot")

cov <- c("insideSeagrass","season")

pred_df <- list()
pval_df <- list()
for (i in cov){
  # fit models----
  mods <- 
    mod_df %>%
    filter(sciName %in% most_common) %>%
    {. ->> mod_df2} %>% 
    group_by(sciName) %>%
    nest() %>%
    mutate(model = map(data, fit_glmm, sciName, cov = i),
           fit_stat = map(model, fit_stat)) %>% 
    unnest(fit_stat) %>% 
    group_by(sciName) 
  
  # get predictions----
  predictions <- 
    mods %>% 
    group_by(sciName, model) %>% 
    dplyr::select(model, sciName) %>% 
    distinct() %>% 
    # nest() %>% 
    group_by(sciName) %>% 
    mutate(pred = map(model, .f = pred_func, spec = sciName, cov = i)) %>% 
    ungroup() %>% 
    dplyr::select(-sciName, -model) %>% 
    unnest(pred) %>% 
    mutate(cov = i)
  
  # get species with perfect separation
  dnc <- 
    predictions %>% 
    filter(!is.finite(conf.high)) %>% 
    pull(sciName)
  
  #extract p values and adjust for multiple comparison
  pvals <- mods %>% 
    filter(str_detect(term, i)) %>% 
    dplyr::select(sciName, p.value) 
  pvals$p.value <- p.adjust(pvals$p.value, method = "BH")
  pvals %<>% mutate(label = ifelse(p.value < 0.05 & p.value > 0.01, "*",
                                      ifelse(p.value < 0.01, "**",
                                             NA)),
                       nFish = ifelse(i == "season", 4.5, 8.25),
                    cov = i)
  pval_df[[i]] <- pvals
  pred_df[[i]] <- predictions
}


(cpue_sg_figure <- 
    ggplot(pred_df$insideSeagrass %>% 
             mutate(sciName = ifelse(sciName == "Atlantic silverside",
                                     "Atl. silverside",sciName))) +
    geom_point(aes(x = sciName, y = fit, color = insideSeagrass),
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = sciName, ymin = conf.low, 
                      ymax = conf.high,
                      color = insideSeagrass),
                  position = position_dodge(width = 0.5),
                  width = 0) +
    geom_text(data = pval_df$insideSeagrass%>% 
                mutate(sciName = ifelse(sciName == "Atlantic silverside",
                                        "Atl. silverside",sciName)), 
              aes(x = sciName, y = nFish, label = label)) +
    ggsci::scale_color_d3() +
    labs(y = "Mean CPUE",
         color = "Sampling habitat") +

    scale_y_continuous(expand = c(0.01, 0.01),
                       limits = c(0, 9),
                       breaks = seq(0,8,2)) +
    # ylim(0, 9) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 0.7),
          legend.position="top",
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          axis.title.y = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)))

(cpue_seas_figure <- 
    ggplot(pred_df$season%>% 
             mutate(sciName = ifelse(sciName == "Atlantic silverside",
                                     "Atl. silverside",sciName))) +
    geom_point(aes(x = sciName, y = fit, color = season),
               position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = sciName, ymin = conf.low, 
                      ymax = conf.high,
                      color = season),
                  position = position_dodge(width = 0.5),
                  width = 0) +
    geom_text(data = pval_df$season%>% 
                mutate(sciName = ifelse(sciName == "Atlantic silverside",
                                        "Atl. silverside",sciName)), 
              aes(x = sciName, y = nFish, label = label)) +
    scale_color_manual(values = c("Summer" = "#925E9FFF",
                                  "Fall" = "#AD002AFF"))+
    labs(y = "Mean CPUE",
         color = "Season") +
    scale_y_continuous(expand = c(0.01, 0.01),
                       limits = c(0, 5)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, 
                                     vjust = 0.7),
          legend.position="top",
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8)))

spec_sg_fig <- 
  cpue_seas_figure + cpue_sg_figure +
  plot_annotation(tag_levels = "A") 

save(spec_sg_fig, file = here::here("data/spec_sg_figure.rdata"))

