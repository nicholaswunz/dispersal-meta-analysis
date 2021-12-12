# Load packages
library(ggplot2)
library(cowplot) 
library(dplyr)
library(sp)
library(tidyverse)
library(ggeffects)

## FIG 1 - STUDY LOCATIONS & DATA SUMMARY ##----------------------------------------------------------------------------------------
# Raw data from disp_analysis.R
world    <- rgdal::readOGR("C:/Users/chwu6540/Dropbox (Personal)/Turtle project/Data/ne_50m_land/ne_50m_land.shp")
world_df <- sp::SpatialPolygonsDataFrame(world, world@data) #turn the data into a spatial data frame

str(ind_data_clean)
ind_map <- ggplot() +
  geom_polygon(data = world_df, aes(x = long, y = lat, group = group), fill = "#d6d6d6", size = 0.2) +
  geom_point(data = ind_data_clean, aes(x = lon, y = lat, colour = thermal_strategy, shape = thermal_strategy), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("#023FA5", "#8E063B")) +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  ylab(NULL) + xlab(NULL) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_fixed() +
  mytheme()

ind_data_clean %>% 
  dplyr::group_by(origin) %>% 
  dplyr::summarise(study_n = length(unique(study_ID)))

pop_map <- ggplot() +
  geom_polygon(data = world_df, aes(x = long, y = lat, group = group), fill = "#d6d6d6", size = 0.2) +
  geom_point(data = pop_data_clean, aes(x = lon_front, y = lat_front, colour = taxa), size = 2, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = pop_data_clean, aes(x = lon_core, y = lat_core, fill = taxa), size = 1.5, shape = 24) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(-60, 90), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  ylab(NULL) + xlab(NULL) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  coord_fixed() +
  mytheme()

# Percentage of taxa represented
taxa_ind_piechart <- ind_data_clean %>%
  dplyr::group_by(taxa) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(per = `n`/sum(`n`),
                label = scales::percent(per)) %>% 
  dplyr::arrange(desc(taxa)) %>% 
  ggplot(aes(x = "", y = per, fill = fct_reorder(taxa, per))) +
  geom_bar( stat="identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Dark2") +
  theme_void() +
  geom_text(aes(label = scales::percent(round(per, 3))), position = position_stack(vjust = 0.5), size = 3)

taxa_pop_piechart <- pop_data_clean %>%
  dplyr::group_by(taxa) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(per = `n`/sum(`n`),
                label = scales::percent(per)) %>% 
  dplyr::arrange(desc(taxa)) %>%
  ggplot(aes(x = "", y = per, fill = fct_reorder(taxa, per))) +
  geom_bar( stat="identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(palette = "Dark2") +
  theme_void() +
  geom_text(aes(label = scales::percent(round(per, 3))), position = position_stack(vjust = 0.5), size = 3)

taxa_piechart <- cowplot::plot_grid(taxa_ind_piechart + ggtitle("Individual movement"), 
                                    taxa_pop_piechart + ggtitle("Population range expansion"), 
                                    ncol = 2, align = 'h', axis = 'tb')

# Combine all plots
cowplot::plot_grid(ind_map + ggtitle("Studies on individual variation in dispersal"), 
                   pop_map + ggtitle("Studies on population range expansion"),
                   taxa_piechart,
                   ncol = 1, align = 'v', axis = 'r', rel_heights = c(1, 1, 0.7, 0.5))