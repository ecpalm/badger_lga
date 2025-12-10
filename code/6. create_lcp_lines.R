
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
# Script description: Create 500-m buffered least cost path lines from predicted  
#                     resistance surface for use in covariate extraction
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("dplyr", "raster", "sf", "terra", "parallel", "doParallel", "stringr",
         "sp", "gdistance", "foreach"),
       require, character.only = T)

# Set previous iteration (e.g., "straight", "lcp_1", "lcp_2", etc.)
iteration_previous <- "straight"

# lcp_iteration (e.g., 1, 2, 3)
iteration_next <- "1"

# Set minimum pairwise geographic distance
min_dist <- 3000

# Import extracted covariates from previous iteration
dat <- 
  readRDS(stringr::str_c("results/extracted_", iteration_previous, "_geomasked_shuffled.rds")) %>% 
  dplyr::filter(euc_geog >= min_dist) %>% 
  dplyr::mutate(index = dplyr::row_number())

# Import resistance raster from previous iteration
res <- terra::rast(stringr::str_c("maps/resistance_", iteration_previous, ".tif"))

# Add water mask and set water to maximum resistance 
res_masked <- res %>% 
  terra::mask(., terra::rast("spatial/raster/mask_water.tif") %>% 
                terra::crop(., res, mask = T), 
              maskvalues = 1, 
              updatevalue = max(terra::values(res, na.rm = T))) %>% 
  terra::app(., function(x)
    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)) + .01) %>%  
  raster::raster()

# Creating transition layer to calculate least cost paths with resistance surface
tr_predr <- gdistance::transition(res_masked, function(x) 1/mean(x), directions=8) %>% 
  gdistance::geoCorrection(., type = "c", multpl = F)

# Begin loop through each pair of locations and calculate buffered least cost path.
# Use as many cores as available to increase speed.

cls <- 100
cl <- parallel::makeCluster(cls)
doParallel::registerDoParallel(cl)

system.time(  
  lcp <- foreach(i = dat$index, .combine = rbind, .verbose = F, 
                 .packages = c("sp", "sf", "raster", "dplyr", "gdistance")) %dopar% 
    {

      y_1 <- dat$y_1[i]
      x_1 <- dat$x_1[i]
      y_2 <- dat$y_2[i]
      x_2 <- dat$x_2[i]
      
      lcp_i <- gdistance::shortestPath(tr_predr, sp::SpatialPoints(cbind(x_1, y_1)), 
                                       sp::SpatialPoints(cbind(x_2, y_2)), 
                                       output = "SpatialLines") %>%
        sf::st_as_sf(.) %>% 
        sf::st_set_crs(., 3005) %>%
        sf::st_buffer(., 500)
      
      return(lcp_i)
    }  
)
# Stop cluster
parallel::stopCluster(cl)

# Save to file.
# This saved file is not included in the GitHub repository due to its large size,
# so you will need to run this script to create it.
# Data frame with extracted covariates from least-cost paths from real dataset 
# is included for running LCP GBM models.
dat %>% 
  dplyr::select(index, euc_gen:euc_geog) %>% 
  dplyr::bind_cols(., lcp) %>% 
  saveRDS(stringr::str_c("results/lcp_", lcp_iteration, "_lines_geomasked_shuffled.rds"))

