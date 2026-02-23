
################################################################################
#
# Article title: Roads, soil, snow, and topography influence genetic 
#                connectivity: A machine learning approach for a peripheral
#                American badger population 
#
# Article authors: Eric Palm, Erin Landguth, Karina Lamy, Jamieson Gorrell,
#                  Richard Weir, Emma Richardson, Krystyn Forbes, Helen Davis,
#                  Joanna Burgar
#
# Script description: Create plots of accumulated local effects (ALE) and
#                     variable relative influence from top model 
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("ale", "caret", "dplyr", "tidyr", "gbm", "stringr", "tibble", "purrr"),
       require, character.only = T) 

# Import saved model
model <- readRDS("models/model_straight.rds")

m_data <- model$trainingData %>% rename(y = .outcome)

# Get the final hyperparameters to use in the model calls below
model$bestTune

# Model call for use in model bootstrapping; the "random variable" is a purely random
# variable that is used to estimate the ALER band (described below)
m_call <- "gbm::gbm(y ~ euc_geog + slope + snow + Organic + Fluvial + Colluvial + 
                    roads_major + lacustrine + random_variable,
                    data = rand_data, 
                    n.trees = 1710, interaction.depth = 2, shrinkage = .005, n.minobsinnode = 25)"

# Model call for full model bootstrapping
m_call_boot <- "gbm::gbm(y ~ euc_geog + slope + snow + Organic + Fluvial + Colluvial + 
                         roads_major + lacustrine,
                         data = boot_data, 
                         n.trees = 1710, interaction.depth = 2, shrinkage = .005, n.minobsinnode = 25)"

# This function allows us to get ALER band around the median predicted genetic distance.
# 95% of random variables used in model bootstrapping have ALE values that fully lay within this band.
# If ALE values for a model variable full completely within this ALER band,
# then it has no greater effect than 95% of purely random variables.
m_alepdist <- ale::ALEpDist(model = model$finalModel, data = m_data, 
                            y_col = "y", pred_type = "response",
                            random_model_call_string = m_call,
                            rand_it = 500, 
                            parallel = parallel::detectCores(logical = FALSE) - 2)

# Conduct full model bootstrapping 500 times to generate confidence intervals around predictions
m_boot <- ale::ModelBoot(model = model$finalModel, data = m_data, 
                         y_col = "y", pred_type = "response",
                         model_call_string = m_call_boot,
                         ale_options = list(max_num_bins = 25),
                         model_packages = c("caret", "gbm"),
                         ale_p = m_alepdist,
                         boot_centre = "median",  
                         output_model_stats = FALSE,
                         output_model_coefs = FALSE,
                         boot_it = 500, 
                         parallel = parallel::detectCores(logical = FALSE) - 2)


# Calculate relative influence of variables for plot ordering
imp <-
  caret::varImp(model) %>% 
  .$importance %>% 
  tibble::rownames_to_column(., "covariate") %>% 
  dplyr::mutate(covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material\n(binary)"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads/n(binary)",
                                             T ~ covariate)) %>% 
  dplyr::mutate(covariate = forcats::fct_reorder(covariate, Overall))

# Preparing ALE data for plotting
for_plot <- 
  get(m_boot, type = "boot") %>%
  purrr::map(., ~dplyr::rename(., !!"ale_x" := names(.)[1])) %>%
  purrr::map_dfr(., ~.x, .id = "covariate") %>%
  dplyr::rename(ale_y = .y_median) %>% 
  dplyr::mutate(ale_x = dplyr::if_else(covariate == "euc_geog", ale_x/1000, ale_x),
                covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material\n(binary)"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads\n(binary)",
                                             T ~ covariate)) %>% 
  dplyr::mutate(covariate = factor(covariate, levels = levels(imp$covariate)))

# Extract ALER band around median predicted value.
ci <-
  m_boot@ale$single@params$y_summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "stat") %>%
  tidyr::pivot_wider(names_from = stat, values_from = y)

# Calculate rugs to add to plot
rugs <- 
  model$trainingData %>% 
  dplyr::select(-.outcome) %>% 
  tidyr::pivot_longer(., cols = dplyr::everything(), names_to = "covariate", values_to = "ale_x") %>% 
  dplyr::mutate(ale_y = 3.19) %>% 
  dplyr::as_tibble() %>% 
  dplyr::mutate(covariate = dplyr::case_when(covariate == "euc_geog" ~ "Geographic distance (km)",
                                             stringr::str_detect(covariate, "Lacust|luvial|Organic") ~
                                               stringr::str_c(covariate, " parent material\n(binary)"),
                                             covariate == "slope" ~ "Slope (%)",
                                             covariate == "snow" ~ "Annual snowfall (mm)",
                                             covariate == "roads_major" ~ "Major roads\n(binary)",
                                             T ~ covariate),
                ale_x = dplyr::if_else(covariate == "Geographic distance (km)", ale_x/1000, ale_x),
                covariate = factor(covariate, levels = levels(imp$covariate)))

# Create ALE plot
for_plot %>%
  ggplot(., aes(x = ale_x, y = ale_y)) +
  geom_ribbon(aes(ymin = .y_lo, ymax = .y_hi), color = NA, fill = "#0072B2", alpha = .2) +
  geom_line(size = 1, color = "#0072B2") +
  geom_hline(data = ci, aes(yintercept = aler_lo), linetype = "dashed",  linewidth = .4, color = "gray40") +
  geom_hline(data = ci, aes(yintercept = aler_hi), linetype = "dashed",  linewidth = .4, color = "gray40") +
  geom_rug(data = rugs, sides = "b", alpha = .05, length = unit(0.05, "npc"), color = "#0072B2",
           show.legend = F) +
  theme_classic(base_size = 18) +
  facet_wrap(~covariate, scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text = element_text(color = "black"),
        legend.position.inside = c(.85, .20),
        text = element_text(family = "Roboto Condensed"),
        strip.text = element_text(vjust = 3),
        plot.margin = margin(r = 15, l = 5, t = 5),
        panel.spacing.y = unit(-.3, "cm"),
        legend.key.width = unit(1, "cm")) +
  labs(color = NULL, fill = NULL, x = NULL, y = "Predicted genetic distance") +
  scale_x_continuous(expand = expansion(mult = c(.02, 0)))

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

