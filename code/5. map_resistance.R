
################################################################################
#
# Article title: Roads, soil, snow, and topography influence genetic 
#                connectivity: A machine learning approach for a peripheral
#                American badger population 
#
# Article authors: Eric Palm, Erin Landguth, Karina Lamy, Jamieson Gorrell,
#                  Richard Weir, Emma Richadson, Krystyn Forbes, Helen Davis
#                  Joanna Burgar
#
# Script description: Create resistance surface by predicting from fitted 
#                     gradient boosting machine
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("dplyr", "stringr", "caret", "gbm", "terra"),
       require, character.only = T)

iteration <- "straight" # or "lcp_1", "lcp_2", etc.

gbm_final <- list.files("models", full.names = T, pattern = iteration)[13] %>% 
  readRDS()

# Use the median geographic distance as a constant value for spatial predictions
env <- list.files("spatial/raster", full.names = T) %>% 
  stringr::str_subset(., "water", negate = T) %>% 
  terra::rast() %>% 
  c(., r_soil) %>% 
  c(., .$slope %>% terra::setValues(., median(gbm_final$trainingData$euc_geog)) %>% 
      `names<-`("euc_geog"))

# Make spatial predictions of landscape resistance from the model using the covariate raster stack
pred_rast <- predict(env[[head(colnames(gbm_final$trainingData), -1)]], gbm_final) %>% 
  terra::mask(., terra::vect("spatial/vector/spatial_extent.shp"))

# Plot to make sure landscape resistance predictions look reasonable
terra::plot(pred_rast)

# Save raster to file
pred_rast %>% 
  terra::writeRaster(., stringr::str_c("maps/resistance_", iteration, ".tif"), 
                     datatype = "FLT4S", 
                     overwrite = T)
