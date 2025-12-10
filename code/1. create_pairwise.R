
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
# Script description: Create data frame of pairwise connections between genetic
#                     sample locations
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("readr",  "dplyr", "gstudio", "sf"),
       require, character.only = T)

# Import genetic samples with geomasked locations and shuffled genotypes
samples <- readr::read_csv("data/samples_geomasked_shuffled.csv") %>% 
  dplyr::mutate(id = as.character(id))

# Function to convert a distance matrix to a tibble where each row is a sample pair
lower_df <-
  function (x) {
    mat <- as.matrix(x)
    ind <- which(lower.tri(mat, diag = F), arr.ind = TRUE)
    nn <- dimnames(mat)
    tibble(id_1 = nn[[1]][ind[, 1]],
           id_2 = nn[[2]][ind[, 2]],
           d = mat[ind])
  }

# Convert columns to "locus" class for use in 'gstudio' package
samples_locus <-
  samples %>% 
  dplyr::select(-c(subregion:y)) %>% 
  dplyr::mutate(across(-id, ~gstudio::locus(., type = "separated"))) %>% 
  as.data.frame() 

# Calculate pairwise genetic distances: euclidean, proportion of shared alleles,
# and Queller-Goodnight genetic relatedness
euc <- gstudio::genetic_distance(samples_locus, mode = "Euclidean", stratum = "id") %>% 
  lower_df() %>% 
  dplyr::rename(euc_gen = d)

Dps <- gstudio::genetic_distance(samples_locus, mode = "Dps", stratum = "id") %>% 
  lower_df() %>% 
  dplyr::rename(Dps = d)

QG <- gstudio::rel_queller(samples_locus) %>% 
  `colnames<-`(sort(samples_locus$id)) %>% 
  `rownames<-`(sort(samples_locus$id)) %>% 
  lower_df() %>% 
  dplyr::rename(QG = d)

# Join the three genetic distance tibbles together, reorder and rename columns.
# Adjust Dps to 1 minus Dps and flip sign of QG so both for both metrics,
# larger values indicate greater genetic distance. 
# Drop pairwise connections between the same animal.
# Add pairwise euclidean geographic distance ("euc_geog")
gen_dists_filtered <-
  purrr::reduce(list(euc, Dps, QG), dplyr::inner_join) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(p1 = pmin(id_1, id_2),
                p2 = pmax(id_1, id_2)) %>% 
  dplyr::select(-c(id_1, id_2)) %>% 
  dplyr::rename(id_1 = p1, id_2 = p2) %>% 
  dplyr::ungroup() %>% 
  dplyr::inner_join(., samples %>% dplyr::select(id_1 = id, x_1 = x, y_1 = y,
                                                    subregion_1 = subregion)) %>% 
  dplyr::inner_join(., samples %>% dplyr::select(id_2 = id, x_2 = x, y_2 = y,
                                                    subregion_2 = subregion)) %>% 
  dplyr::mutate(across(Dps, ~ 1 - .)) %>% 
  dplyr::mutate(across(QG, ~ -1 * .)) %>% 
  dplyr::filter(!id_1 == id_2) %>% 
  dplyr::mutate(euc_geog = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2))

# Save file
saveRDS(gen_dists_filtered, "results/samples_pairwise_geomasked_shuffled.rds")
