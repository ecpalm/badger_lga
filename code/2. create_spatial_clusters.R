
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
# Script description: Create spatial clusters for use in cross-validation
#                     within gradient boosting machines
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("readr",  "dplyr", "CAST", "sf"),
       require, character.only = T)

# Import geomasked genetic sample locations
samples <- readr::read_csv("data/samples_geomasked_shuffled.csv")

# Import spatial extent
extent <- sf::read_sf("spatial/vector/spatial_extent.shp")

# Use K-fold nearest neighbor distance matching to create 6 spatial clusters
knn_folds <- CAST::knndm(tpoints = sf::st_as_sf(samples, coords = c("x", "y"), crs = 3005), 
                         k = 6, 
                         modeldomain = extent %>% 
                           sf::st_transform(3005))

# Add the cluster numbers to the samples and create a list of six tibbles, 
# each of which represents the badger IDs in a spatial cluster. This will be
# used in the model cross validation in the 'caret' package.
samples %>%
  dplyr::mutate(cluster = knn_folds$clusters) %>% 
  dplyr::select(id, cluster) %>% 
  dplyr::nest_by(cluster) %>% 
  dplyr::pull(data) %>% 
  saveRDS(., "results/spatial_clusters_knndm_index_6_geomasked.rds")
