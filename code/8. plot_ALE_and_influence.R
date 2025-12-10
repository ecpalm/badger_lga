
################################################################################
#
# Article title: Roads, soil, snow, and topography influence genetic 
#                connectivity: A machine learning approach for a peripheral
#                American badger population 
#
# Article authors: Eric Palm, Erin Landguth, Karina Lamy, Jamieson Gorrell,
#                  Richard Weir, Emma Richadson, Krystyn Forbes, Helen Davis,
#                  Joanna Burgar
#
# Script description: Create plots of accumulated local effects (ALE) and
#                     variable relative influence from top model 
#
# Script author: Eric Palm
#
################################################################################

# This script requires "ale" package version 0.3.1
# Subsequent versions changed the syntax and the ale::ale function below will not work

# Load required packages
sapply(c("ale", "caret", "dplyr", "tidyr", "gbm", "stringr"),
       require, character.only = T) 

# Import saved model
model <- readRDS("models/model_straight.rds")

# Generate accumulated local effects, using parallel argument to speed things up
ale_m <- ale::ale(data = model$trainingData, model = model$finalModel, y_col = ".outcome", pred_type = "response",
                  parallel = 50, boot_it = 100,
                  seed = 1234)

# Calculate relative influence of variables for plot ordering
imp <-
  caret::varImp(model) %>% 
  .$importance %>% 
  tibble::rownames_to_column(., "covariate") %>% 
  dplyr::mutate(covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Coll|Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads",
                                             T ~ covariate)) %>% 
  dplyr::mutate(covariate = forcats::fct_reorder(covariate, Overall))

# Preparing ALE data for plotting
for_plot <- 
  ale_m$data %>% 
  dplyr::bind_rows(., .id = "covariate") %>% 
  dplyr::mutate(ale_x = dplyr::if_else(covariate == "euc_geog", ale_x/1000, ale_x),
                covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Coll|Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads",
                                             T ~ covariate)) %>% 
  dplyr::mutate(covariate = factor(covariate, levels = levels(imp$covariate)))

# Calculate confidence intervals around median predicted value
ci <-
  tibble::enframe(ale_m$y_summary) %>% 
  dplyr::distinct() %>% 
  tidyr::pivot_wider(names_from = name)

# Calculate rugs to add to plot
rugs <- 
  model$trainingData %>% 
  dplyr::select(-.outcome) %>% 
  tidyr::pivot_longer(., cols = dplyr::everything(), names_to = "covariate", values_to = "ale_x") %>% 
  dplyr::mutate(ale_y = 3.19) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Coll|Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads",
                                             T ~ covariate),
                ale_x = dplyr::if_else(covariate == "Geographic distance (km)", ale_x/1000, ale_x),
                covariate = factor(covariate, levels = levels(imp$covariate)))

# Create ALE plot
for_plot %>% 
  ggplot(., aes(x = ale_x, y = ale_y)) +
  geom_rect(data = ci, aes(xmin = -Inf, xmax = Inf, ymin = med_lo, ymax = med_hi), 
            fill = "#0072B2", color = NA, alpha = .2, inherit.aes = F) +
  geom_line(linewidth = 1, color = "#0072B2") +
  geom_hline(data = ci, aes(yintercept = `50%`), 
             linetype = "dotted",  linewidth = .8, color = "#0072B2") +
  geom_rug(data = rugs, sides = "b", alpha = .05, length = unit(0.05, "npc"), color = "#0072B2",
           show.legend = F) +
  theme_classic(base_size = 20) +
  facet_wrap(~forcats::fct_rev(covariate), scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text = element_text(color = "black"),
        legend.position.inside = c(.85, .20),
        text = element_text(family = "Roboto Condensed"),
        strip.text = element_text(vjust = 1.5),
        plot.margin = margin(r = 15, l = 5, t = 5),
        panel.spacing.y = unit(-.3, "cm"),
        legend.key.width = unit(1, "cm")) +
  labs(color = NULL, 
       fill = NULL, 
       x = NULL, 
       y = "Predicted genetic distance")

# Save plot to file
ggsave("figures/ALE_straight.tiff", height = 8, width = 9, compression = "lzw")


# Plot variable relative influence
imp %>% 
  ggplot(., aes(y = covariate, x = Overall)) +
  geom_col(aes(fill = effect)) +
  theme_classic(base_size = 20) +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        legend.background = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        strip.placement = "outside",
        strip.text = element_text(size = 15),
        text = element_text(family = "Roboto Condensed")) +
  labs(y = NULL, fill = NULL,
       x = "Relative influence on\npredicted genetic distance") +
  scale_x_continuous(expand = expansion(mult = c(0, .02))) +
  ggokabeito::scale_fill_okabe_ito(order = c(3,6))

# Save plot to file
ggsave("figures/relative_influence_straight.tiff", height = 7, width = 7, dpi = 600, compression = "lzw")

