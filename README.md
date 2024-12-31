
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ensembleDownscaleR

<!-- badges: start -->

[![R-CMD-check](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WyattGMadden/ensembleDownscaleR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of ensembleDownscaleR is to …

## Installation

You can install the development version of ensembleDownscaleR like so:

``` r
remotes::install_github("WyattGMadden/ensembleDownscaleR")
```

## Example

``` r
library(ensembleDownscaleR)
```

### Stage 1

``` r

?cmaq_aqs_matched

cmaq_fit <- grm(
    Y = cmaq_aqs_matched$pm25,
    X = cmaq_aqs_matched$ctm,
    L = cmaq_aqs_matched[, c("elevation", "forestcover",
                             "hwy_length", "lim_hwy_length", 
                             "local_rd_length", "point_emi_any")],
    M = cmaq_aqs_matched[, c("tmp", "wind")],
    n.iter = 50,
    burn = 10,
    thin = 4,
    covariance = "matern",
    matern.nu = 0.5,
    coords = cmaq_aqs_matched[, c("x", "y")],
    space.id = cmaq_aqs_matched$space_id,
    time.id = cmaq_aqs_matched$time_id,
    spacetime.id = cmaq_aqs_matched$spacetime_id,
    verbose.iter = 10
)

?modis_aqs_matched

modis_fit <- grm(
    Y = modis_aqs_matched$pm25,
    X = modis_aqs_matched$aod,
    L = modis_aqs_matched[, c("elevation", "forestcover",
                              "hwy_length", "lim_hwy_length", 
                              "local_rd_length", "point_emi_any")],
    M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                              "windaod", "elevationaod")],
    n.iter = 50,
    num_neighbors = 6,
    burn = 10,
    thin = 4,
    coords = modis_aqs_matched[, c("x", "y")],
    space.id = modis_aqs_matched$space_id,
    time.id = modis_aqs_matched$time_id,
    spacetime.id = modis_aqs_matched$spacetime_id
)
```

### Stage 2

``` r

?cmaq_full
 
cmaq_pred <- grm_pred(
    grm.fit = cmaq_fit,
    X = cmaq_full$ctm,
    L = cmaq_full[, c("elevation", "forestcover",
                      "hwy_length", "lim_hwy_length", 
                      "local_rd_length", "point_emi_any")],
    M = cmaq_full[, c("tmp", "wind")],
    coords = cmaq_full[, c("x", "y")],
    time.id = cmaq_full$time_id,
    space.id = cmaq_full$space_id,
    spacetime.id = cmaq_full$spacetime_id,
    n.iter = 10,
    verbose = T
)

?modis_full

modis_pred <- grm_pred(
    grm.fit = modis_fit,
    X = modis_full$aod, 
    L = modis_full[, c("elevation", "forestcover",
                       "hwy_length", "lim_hwy_length", 
                       "local_rd_length", "point_emi_any")],
    M = modis_full[, c("tmp", "wind", "cmaq", "tempaod", 
                       "windaod", "elevationaod")],
    coords = modis_full[, c("x", "y")],
    time.id = modis_full$time_id, 
    space.id = modis_full$space_id, 
    spacetime.id = modis_full$spacetime_id,
    n.iter = 10,
    verbose = T)
```

### Stage 3

``` r

cv_id_cmaq_ord <- create_cv(
    space.id = cmaq_aqs_matched$space_id,
    time.id = cmaq_aqs_matched$time_id, 
    type = "ordinary"
)

cmaq_fit_cv <- grm_cv(
    Y = cmaq_aqs_matched$pm25,
    X = cmaq_aqs_matched$ctm,
    cv.object = cv_id_cmaq_ord,
    L = cmaq_aqs_matched[, c("elevation", "forestcover",
                             "hwy_length", "lim_hwy_length", 
                             "local_rd_length", "point_emi_any")],
    M = cmaq_aqs_matched[, c("tmp", "wind")],
    n.iter = 50,
    burn = 10,
    thin = 4,
    coords = cmaq_aqs_matched[, c("x", "y")],
    space.id = cmaq_aqs_matched$space_id,
    time.id = cmaq_aqs_matched$time_id,
    spacetime.id = cmaq_aqs_matched$spacetime_id
)

cv_id_modis_ord <- create_cv(
    create.from = cv_id_cmaq_ord,
    space.id = modis_aqs_matched$space_id,
    time.id = modis_aqs_matched$time_id,
)

modis_fit_cv <- grm_cv(
    Y = modis_aqs_matched$pm25,
    X = modis_aqs_matched$aod,
    cv.object = cv_id_modis_ord,
    L = modis_aqs_matched[, c("elevation", "forestcover",
                              "hwy_length", "lim_hwy_length", 
                              "local_rd_length", "point_emi_any")],
    M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                              "windaod", "elevationaod")],
    n.iter = 50,
    burn = 10,
    thin = 4,
    coords = modis_aqs_matched[, c("x", "y")],
    space.id = modis_aqs_matched$space_id,
    time.id = modis_aqs_matched$time_id,
    spacetime.id = modis_aqs_matched$spacetime_id
)
```

### Stage 4

``` r

ensemble_fit <- ensemble_spatial(
    grm.fit.cv.1 = cmaq_fit_cv,
    grm.fit.cv.2 = modis_fit_cv,
    n.iter = 50, 
    burn = 10, 
    thin = 4,
    tau.a = 0.001, 
    tau.b = 0.001, 
    theta.tune = 0.2, 
    theta.a = 5, 
    theta.b = 0.05
)
```

### Stage 5

``` r

weight_preds <- weight_pred(
    ensemble.fit = ensemble_fit,
    coords = cmaq_full[, c("x", "y")],
    space.id = cmaq_full$space_id
)
```

### Stage 6

``` r

results <- gap_fill(
    grm.pred.1 = cmaq_pred,
    grm.pred.2 = modis_pred,
    weights = weight_preds
)
```
