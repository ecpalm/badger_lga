
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
# Script description: Run gradient boosting machines in 'caret' incorporating
#                     spatial cross-validation and variable selection
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("doParallel", "caret", "raster", "dplyr", "CAST", "gbm", "stringr"),
       require, character.only = T) 

# Set iteration (e.g., "straight", "lcp_1", "lcp_2", etc.)
iteration <- "straight"

# Set minimum pairwise geographic distance
min_dist <- 3000

dat <- 
  readRDS(stringr::str_c("results/extracted_", iteration, ".rds")) %>% 
  dplyr::filter(euc_geog >= min_dist) %>% 
  dplyr::mutate(index = dplyr::row_number())

# Select all covariates from Table 1
traindat <- dat %>% 
  dplyr::select(euc_geog, slope, snow, canopy_cover, # continuous
                Shrubland, Grassland, Wetland, Colluvial, Fluvial, Glaciofluvial, Lacustrine, Organic, # binary soil
                roads_major) # binary

# Load spatial CV clusters
cluster_list <- readRDS("results/spatial_clusters_knndm_index_6.rds")

# Empty list for storing results of for loop
indices_cv <- list()

# For each spatial cluster, this identifies the indices in the input data that 
# do NOT contain a genetic sample from that cluster 
for (i in 1:length(cluster_list)){
  indices_cv[[i]] <- dat %>%
    dplyr::filter(!(id_1 %in% cluster_list[[i]]$id |
                      id_2 %in% cluster_list[[i]]$id)) %>%
    dplyr::pull(index)
}

# This is the format that the 'caret' and 'CAST' packages use
indices_cv <- Filter(function(x) length(x) < nrow(dat), indices_cv)

# specify model cross-validation parameters here
ctrl_fit <- caret::trainControl(index = indices_cv,
                                allowParallel = T,
                                number = length(indices_cv),
                                method = "cv",
                                classProbs = FALSE,
                                savePredictions = T,
                                seeds = NULL)

# Choose how many clusters you want to run for parallelization
cls <- 60
cl <- parallel::makeCluster(cls)
doParallel::registerDoParallel(cl)

# Run forward feature selection in 'caret' using the 'CAST' package 
# Set tuning hyperparameters for 'caret::train'
# This example uses a narrow search grid to facilitate faster run time 
n_trees <- seq(50, 2000, 50)
learning_rate <- 0.005
int_depth <- c(1, 2)
n_min_obs <- c(20)

system.time(
  gbm_ffs <- CAST::ffs(predictors = traindat, 
                       response = dat$euc_gen,
                       method = "gbm",
                       trControl = ctrl_fit,
                       metric = "RMSE",
                       minVar = 2,
                       seed = 1234,
                       tuneGrid=expand.grid(interaction.depth = int_depth,
                                            n.trees = n_trees,
                                            shrinkage = learning_rate,
                                            n.minobsinnode = n_min_obs))
)

# Now that you have the feature-selected variables, run a model that only
# uses those variables, and conduct a slightly more thorough grid search

# Run a wider search grid for n.trees on final model
n_trees <- seq(50, 3000, 10)
learning_rate <- .005
int_depth <- c(1,2)
n_min_obs <- c(5, 15, 25, 50)

# Set seed for reproducibility
set.seed(1234)

system.time(
  gbm_model <- caret::train(x = traindat[, gbm_ffs$selectedvars], 
                            y = dat$euc_gen,
                            method = "gbm",
                            trControl = ctrl_fit,
                            metric = "RMSE",
                            tuneGrid=expand.grid(interaction.depth = int_depth,
                                                 n.trees = n_trees,
                                                 shrinkage = learning_rate,
                                                 n.minobsinnode = n_min_obs))
)

# Check RMSE_test value of the final model
caret::getTrainPerf(gbm_model)

# Get best set of hyperparameters from the final model
best <- gbm_model$bestTune

# Run a model with only the best tune to reduce file size for github repository
ctrl_fit$savePredictions <- "none"

set.seed(1234)
system.time(
  gbm_final <- caret::train(x = traindat[, gbm_ffs$selectedvars], 
                            y = dat$euc_gen,
                            method = "gbm",
                            trControl = ctrl_fit,
                            metric = "RMSE",
                            tuneGrid = best)
)

# Stop cluster
parallel::stopCluster(cl)

# Save the model output
saveRDS(gbm_final, stringr::str_c("models/model_", iteration, ".rds"))
