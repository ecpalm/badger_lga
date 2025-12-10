
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
# Script description: Extract covariates within 500-m buffered straight-line
#                     pairwise transects
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("dplyr", "sf", "data.table", "sfheaders", "parallel", "doParallel",
         "raster", "terra", "exactextractr", "foreach"),
       require, character.only = T)

# Load geomasked and shuffled pairwise genetic data
gen_table <- readRDS("results/samples_pairwise_geomasked_shuffled.rds") %>% 
  dplyr::mutate(index = dplyr::row_number())

############# CREATE LINE SEGMENTS
# Use the 'data.table' and 'sfheaders' packages to quickly create straight-line transects
dt <- data.table::as.data.table(gen_table)

## To use `sfheaders` the data needs to be in long form
dt1 <- dt[, .(index, x = x_1, y = y_1)]
dt2 <- dt[, .(index, x = x_2, y = y_2)]

## Add a 'sequence' variable so we know which one comes first
dt1[, seq := 1L ]
dt2[, seq := 2L ]

# Bind the rows and make sure it's in the correct order
dt <- data.table::rbindlist(list(dt1, dt2), use.names = TRUE)
data.table::setorder(dt, index, seq)

# Create a sf linestring object and buffer by 500 m to create polygons
line_sf <- sfheaders::sf_linestring(
  obj = dt
  , x = "x"
  , y = "y"
  , linestring_id = "index"
) %>% 
  sf::st_set_crs(3005) %>% 
  sf::st_buffer(., 500)

# Load environmental covariate data and crop to transect polygon spatial extent
env <- list.files("spatial/raster", full.names = T) %>% 
  stringr::str_subset(., "water", negate = T) %>% 
  terra::rast() %>% 
  terra::crop(., line_sf) %>% 
  raster::stack()

# Extract mean covariate values within buffered straight lines
# set up 'foreach' loop here
noCLS <- raster::nlayers(env)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

system.time(  
  env_list <- foreach(i=seq_along(names(env)), .verbose = F, 
                       .packages = c("sf", "raster", "exactextractr", "terra", "dplyr")) %dopar% 
    {
      exactextractr::exact_extract(x = terra::rast(env)[[i]], y = line_sf, append_cols = "index", fun = "mean") %>% 
        `names<-`(c("index", names(env[[i]])))
      
    }
)

# Stop the cluster
parallel::stopCluster(cl)

# Rename columns in the output data frame of covariate values
df_extracted <- Reduce(dplyr::inner_join, env_list) %>% 
  dplyr::inner_join(gen_table, .)

# Save the file
saveRDS(df_extracted, "results/extracted_straight_geomasked_shuffled.rds")
