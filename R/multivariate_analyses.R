library(tidyverse)
library(vegan)
library(ggvegan)
library(GGally)
library(ggrepel)
library(ggtext)
library(mgcv)
library(ggordiplots)

load(here::here("data/fish_processed_aggregated.rdata"))
load(here::here("data/biodiversity_table.rdata"))

specs_wide2 %<>% 
  mutate(tot = rowSums(dplyr::select(.,`American halfbeak`:Tautog))) %>% #this step highlights where there are rows with zero catches
  filter(tot > 0) %>% 
  dplyr::select(-tot)

get_ellipsoids <- function(nmds, input, var){
  input[[paste(var)]] <- as.factor(input[[paste(var)]])
  plot.new()
  ord <- ordiellipse(nmds, input[[paste(var)]], 
                     display = "sites", kind ="sd", 
                     conf = 0.95, label = T)
  dev.off()
  
  df_ell <- tibble()
  for(g in levels(input[[paste(var)]])){
    if(g != "" && (g %in% names(ord))){
      df_ell <- rbind(df_ell, 
                      cbind(as.data.frame(
                        with(
                          input[input[[paste(var)]] == g,],
                          vegan:::veganCovEllipse(ord[[g]]$cov,
                                          ord[[g]]$center,
                                          ord[[g]]$scale)
                        )
                      ),
                      treatment = g)
      )
    }
  }
  return(df_ell)
}


# extract community matrix and transform (log(x) + 1 for x > 0; zeros left alone)
spp <- specs_wide2 %>% 
  dplyr::select(`American halfbeak`:Tautog) %>% 
  as.data.frame() %>% 
  decostand(.,method = "log")

# extract matrix of environmental covariates
env <- specs_wide2 %>% 
  dplyr::select(-`American halfbeak`:-Tautog) %>% 
  mutate(meadow = factor(meadow),
         season = factor(str_to_title(season),
                         levels = c("Summer", "Fall")),
         month = factor(case_when(month == 5 ~ "May",
                                  month == 6 ~ "June",
                                  month == 9 ~ "Sept.",
                                  month == 10 ~ "Oct."),
                        levels = c("May","June","Sept.","Oct.")),
         site2 = ifelse(site2 == "HI", "Hog Island - Veg.",
                        ifelse(site2 == "HI-Bare", "Hog Island - Unveg.",
                               ifelse(site2 == "SB", "South Bay - Veg.",
                                      ifelse(site2 == "SB-bare", "South Bay - Unveg.", NA)))),
         insideSeagrass = factor(ifelse(insideSeagrass == 1, "Seagrass","Unvegetated"),
                                 levels = c("Seagrass","Unvegetated"))) %>% 
  dplyr::select(year, meadow, month, site2, site,
                site_month, season, insideSeagrass)

# PERMANOVA----
mod <- 
  adonis2(spp ~ insideSeagrass + season + year, 
          permutations = how(blocks = env$meadow,
                             plots = Plots(env$site),
                             within = Within("free"),
                             nperm = 5000),
        data = env)
mod
# NDMS----
mds_mod <- metaMDS(spp,  trymax = 1000, noshare = 0.01, k = 2)

# Plotting----
bcd_stand <- bind_cols(spp, env)

# Extract NMDS ellipse points
nmds_ellipsoids <- NULL
for (i in c("meadow","season", "insideSeagrass")){
  intermediate <- get_ellipsoids(nmds = mds_mod, 
                                 input = bcd_stand,
                                 var  = i) %>% 
    mutate(Var = i) 
  
  assign('nmds_ellipsoids', rbind(nmds_ellipsoids, intermediate))
}

# Calculate centroids from NMDS plot and extract points
nmds_df <- mds_mod$points %>%
  as_tibble %>% 
  bind_cols(., bcd_stand %>% dplyr::select(season,year, meadow, insideSeagrass)) %>%
  mutate_at(vars(season, meadow, insideSeagrass), as.factor) 
  
nmds_seas <- nmds_df %>% 
  group_by(season) %>% 
  dplyr::summarise(MDS2 = mean(MDS2),
                   MDS1 = mean(MDS1))

nmds_sg <- nmds_df %>% 
  group_by(insideSeagrass) %>% 
  dplyr::summarise(MDS2 = mean(MDS2),
                   MDS1 = mean(MDS1))

# Habitat plot----
sg_mds <- 
  ggplot() + 
  geom_point(data = nmds_df, 
             aes(x = MDS1, y = MDS2, color = insideSeagrass), shape = 1) +
  geom_point(data = nmds_sg,
             aes(x = MDS1, y = MDS2, color = insideSeagrass),
             size = 2) +
  # scale_shape_manual(values = c(3,2)) +
  geom_polygon(data = nmds_ellipsoids %>% 
                 filter(Var == "insideSeagrass"),
               aes(x = NMDS1, y = NMDS2, 
                   color = treatment, 
                   group = treatment), 
               fill = "transparent") +
  ggsci::scale_color_d3() +
  dream::theme_fade() +
  labs(color = "Sample habitat") +
  theme(panel.border = element_rect(colour = "black", fill=NA,
                                    size=1),
        # legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.title =  element_text(size = 8),
        axis.text =  element_text(size = 8))
# season plot----
seas_mds <- 
  ggplot() + 
  geom_point(data = nmds_df, 
             aes(x = MDS1, y = MDS2, color = season), shape = 1) +
  geom_point(data = nmds_seas,
             aes(x = MDS1, y = MDS2, color = season),
             size = 2) +
  # scale_shape_manual(values = c(3,2)) +
  geom_polygon(data = nmds_ellipsoids %>% 
                 filter(Var == "season"),
               aes(x = NMDS1, y = NMDS2, 
                   color = treatment, 
                   group = treatment), 
               fill = "transparent") +
  labs(color = "Season") +
  ggsci::scale_color_jama() +
  dream::theme_fade() +
  theme(panel.border = element_rect(colour = "black", fill=NA,
                                    size=1),
        # legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        axis.title =  element_text(size = 8),
        axis.text =  element_text(size = 8))

# year plot----
(yr_mds <- 
  ggplot() + 
  geom_point(data = nmds_df, 
             aes(x = MDS1, y = MDS2, color = factor(year,
                                                    levels = 2012:2018)), size = 2) +
  ggsci::scale_color_lancet() +
  dream::theme_fade() +
  labs(color = "Year",
       alpha = "Year") +
  theme(panel.border = element_rect(colour = "black", fill=NA,
                                    size=1),
        # legend.position = "bottom",
        legend.spacing.x = unit(0.05, "cm"),
        legend.title = element_text(size = 8),
        # legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.title =  element_text(size = 8),
        axis.text =  element_text(size = 8)))
 
all_nmds_figs <- 
  sg_mds + seas_mds + yr_mds +
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

ggsave(all_nmds_figs, 
       filename = here::here("figures/nmds_figs.png"),
       device = "png",
       dpi = 300,
       width = 8, height = 5.5)
