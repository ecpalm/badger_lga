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
# Script description: Extract covariates within 500-m buffered least-cost path
#                     pairwise transects
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("dplyr", "raster", "sf", "terra", "parallel", "doParallel", "stringr",
        "foreach", "exactextractr"),
       require, character.only = T)

# Set LCP iteration
lcp_iteration <- 1

# Import saved buffered LCP lines.
# This file is not provided in the Github repository due to size, so you will
# need to create it and save it using the 'create_lcp_lines.R' script
lcps <- 
  readRDS(stringr::str_c("results/lcp_", lcp_iteration, "_lines_geomasked_shuffled.rds")) %>% 
  sf::st_sf() %>% 
  sf::st_make_valid() %>% 
  dplyr::mutate(index = dplyr::row_number())

# List of annual raster stacks, named by year
env <- list.files("spatial/raster/", full.names = T) %>% 
  stringr::str_subset(., "water", negate = T) %>% 
  raster::stack()

# Extract mean covariate values within buffered least cost paths
# set up 'foreach' loop here
noCLS <- raster::nlayers(env)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

system.time(  
  env_list <- foreach(i=seq_along(names(env)), .verbose = F, 
                      .packages = c("sf", "raster", "exactextractr", "terra", "dplyr")) %dopar% 
    {
      exactextractr::exact_extract(x = terra::rast(env)[[i]], y = lcps, append_cols = "index", fun = "mean") %>% 
        `names<-`(c("index", names(env[[i]])))
      
    }
)

# Stop the cluster
parallel::stopCluster(cl)

# Rename columns in the output data frame of covariate values
df_extracted <- Reduce(dplyr::inner_join, env_list) %>% 
  dplyr::inner_join(sf::st_drop_geometry(lcps), .)

saveRDS(df_extracted, stringr::str_c("results/extracted_lcp_", lcp_iteration, "_geomasked_shuffled.rds"))

# Return to '4. model_gbm' script and use this saved file as input data for the 
# next iteration of LCP models 
