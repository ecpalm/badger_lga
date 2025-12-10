# Roads, soil, snow, and topography influence genetic connectivity: A machine learning approach for a peripheral American badger population  
This repository includes data and code for landscape genetics analysis of American badgers in British Columbia, Canada. 

# NOTE 
- Because badgers in BC are an endangered species, we do not include data with actual GPS locations from badger genetic samples. 
- Instead, we include data with randomly generated GPS locations for use in scripts that require GPS locations.
- For running gradient boosting machines and creating spatial predictions of landscape resistance and connectivity (which don't require GPS locations), we include the actual badger data.
- Running connectivity algorithm (resistant kernels) will require installation of `UNICOR`, available for download here: https://github.com/ComputationalEcologyLab/UNICOR

# Details
DOI: https://doi.org/10.5281/zenodo.17885809


Please open the `badger_lga.Rproj` file to start RStudio before opening individual scripts to ensure that relative file paths in the code work correctly.


Here is a list of script files with descriptions: 

`1. create_pairwise.R` – Create data frame of pairwise connections between genetic sample locations.

`2. create_spatial_clusters.R` – Create spatial clusters for use in cross-validation within gradient boosting machines.

`3. extract_straight_line.R` – Extract covariates within 500-m buffered straight-line pairwise transects.

`4. model_gbm.R` – Run gradient boosting machines in `caret` incorporating spatial cross-validation and variable selection.

`5. map_resistance.R` – Create resistance surface by predicting from fitted gradient boosting machine.

`6. create_lcp_lines.R` – Create buffered least-cost path lines from predicted resistance surface for use in covariate extraction

`7. extract_lcp.R` – Extract covariates within 500-m buffered least-cost path pairwise transects.

`8. plot_ALE_and_influence.R` – Create plots of accumulated local effects (ALE) and variable relative influence from top model.

`9. unicor.R` – Prepare resistance surface for resistant kernel connectivity algorithm in UNICOR and save resistant kernel raster output from UNICOR.


