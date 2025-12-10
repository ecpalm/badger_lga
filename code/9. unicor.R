
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
# Script description: Prepare resistance surface for resistant kernel 
#                     connectivity algorithm in UNICOR and save resistant
#                     kernel raster output from UNICOR
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("ale", "caret", "dplyr", "tidyr", "gbm", "stringr"),
       require, character.only = T) 

# Import predicted resistance surface
r <- terra::rast("maps/resistance_straight.tif") 

# Set water values to maximum resistance and normalize between 1 and 101
r_mask_norm <- terra::rast("maps/resistance_straight.tif")  %>% 
  terra::mask(., terra::rast("spatial/raster/mask_water.tif") %>% 
                terra::crop(., r, mask = T), maskvalues = 1, updatevalue = max(values(., na.rm = T))) %>% 
  terra::app(., function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))*100 + 1) 

# Write this normalized raster to an .asc file for UNICOR
terra::writeRaster(r_mask_norm, "unicor/resistance_straight.asc", datatype = "FLT4S", overwrite = T) 

# Pull out the different components of the the .asc text/header file
tx <- readLines("unicor/resistance_straight.asc")
header <- c(tolower(tx[(1:5)]), tx[6]) 
bdy <- (tx)[-(1:6)]

# Then rewrite the .asc file in the preferred format for UNICOR
# This saved file not included in github repository due to its large size,
# so you will need to run this script to create it.
writeLines(append(header, bdy, after = 6), "unicor/resistance_straight.asc")

# Sample 2000 start locations with a regular distribution from within the study area.
# spatial extent, but randomly jitter the locations a bit.
set.seed(1234)

sf::read_sf("spatial/vector/spatial_extent.shp") %>% 
  dplyr::mutate(area_km2 = as.numeric(units::set_units(sf::st_area(.), "km2"))) %>% 
  sf::st_buffer(-200) %>% 
  sf::st_sample(., 2000, type = "regular") %>% 
  sf::st_sf() %>% 
  dplyr::mutate(X = sf::st_coordinates(.)[,1] + runif(nrow(.), min = -1500, max = 1500),
                Y = sf::st_coordinates(.)[,2] + runif(nrow(.), min = -1500, max = 1500)) %>% 
  sf::st_drop_geometry() %>% 
  readr::write_csv(., "unicor/start_points.csv")

################################################################################

# Run resistant kernel connectivity algorithm in UNICOR, available here:
# https://github.com/ComputationalEcologyLab/UNICOR

# Install UNICOR in "unicor" directory
# use command line to navigate to "unicor" directory. 
# Then run the following line:
# python UNICOR.py resistant_kernel.rip

################################################################################

# Import the resistant kernel output from UNICOR, set the coordinate reference
# system, mask it by the study area extent, and write it to a file
terra::rast("unicor/resistance_straight.asc_start_points.csv.addedpaths.txt") %>% 
  `crs<-`("EPSG:3005") %>% 
  terra::mask(., terra::vect("spatial/vector/spatial_extent.shp")) %>%
  terra::writeRaster(., "maps/resistant_kernel_straight.tif", 
                     datatype = "INT2U", overwrite = T)
