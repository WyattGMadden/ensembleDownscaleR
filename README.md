
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

### Stage 1

``` r
library(ensembleDownscaleR)

?cmaq_aqs_matched

ctm_fit <- grm(Y = cmaq_aqs_matched$pm25,
               X = cmaq_aqs_matched$ctm,
               L = cmaq_aqs_matched[, c("elevation", "forestcover",
                                        "hwy_length", "lim_hwy_length", 
                                        "local_rd_length", "point_emi_any")],
               M = cmaq_aqs_matched[, c("tmp", "wind")],
               n.iter = 100,
               burn = 20,
               thin = 1,
               covariance = "matern",
               matern.nu = 0.5,
               coords = cmaq_aqs_matched[, c("x", "y")],
               space.id = cmaq_aqs_matched$space_id,
               time.id = cmaq_aqs_matched$time_id,
               spacetime.id = cmaq_aqs_matched$spacetime_id,
               verbose.iter = 10)




cv_id_ctm_ord <- create_cv(space.id = cmaq_aqs_matched$space_id,
                           time.id = cmaq_aqs_matched$time_id, 
                           type = "ordinary")



ctm_fit_cv <- grm_cv(Y = cmaq_aqs_matched$pm25,
                     X = cmaq_aqs_matched$ctm,
                     cv.object = cv_id_ctm_ord,
                     L = cmaq_aqs_matched[, c("elevation", "forestcover",
                                      "hwy_length", "lim_hwy_length", 
                                      "local_rd_length", "point_emi_any")],
                     M = cmaq_aqs_matched[, c("tmp", "wind")],
                     n.iter = 500,
                     burn = 100,
                     thin = 4,
                     coords = cmaq_aqs_matched[, c("x", "y")],
                     space.id = cmaq_aqs_matched$space_id,
                     time.id = cmaq_aqs_matched$time_id,
                     spacetime.id = cmaq_aqs_matched$spacetime_id)


?modis_aqs_matched

cv_id_maia_ord <- create_cv(space.id = modis_aqs_matched$space_id,
                            time.id = modis_aqs_matched$time_id,
                            type = "ordinary")

maia_fit <- grm(Y = modis_aqs_matched$pm25,
                X = modis_aqs_matched$aod,
                L = modis_aqs_matched[, c("elevation", "forestcover",
                                  "hwy_length", "lim_hwy_length", 
                                  "local_rd_length", "point_emi_any")],
                M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                                  "windaod", "elevationaod")],
                n.iter = 500,
                num_neighbors = 6,
                burn = 100,
                thin = 4,
                coords = modis_aqs_matched[, c("x", "y")],
                space.id = modis_aqs_matched$space_id,
                time.id = modis_aqs_matched$time_id,
                spacetime.id = modis_aqs_matched$spacetime_id)

maia_fit_cv <- grm_cv(Y = modis_aqs_matched$pm25,
                      X = modis_aqs_matched$aod,
                      cv.object = cv_id_maia_ord,
                      L = modis_aqs_matched[, c("elevation", "forestcover",
                                        "hwy_length", "lim_hwy_length", 
                                        "local_rd_length", "point_emi_any")],
                      M = modis_aqs_matched[, c("tmp", "wind", "cmaq", "tempaod", 
                                        "windaod", "elevationaod")],
                      n.iter = 500,
                      burn = 100,
                      thin = 4,
                      coords = modis_aqs_matched[, c("x", "y")],
                      space.id = modis_aqs_matched$space_id,
                      time.id = modis_aqs_matched$time_id,
                      spacetime.id = modis_aqs_matched$spacetime_id)
```

### Stage 2

``` r
?cmaq_full

ctm_pred <- grm_pred(grm.fit = ctm_fit,
                     X.pred = cmaq_full$ctm,
                     L.pred = cmaq_full[, c("elevation", "forestcover",
                                            "hwy_length", "lim_hwy_length", 
                                            "local_rd_length", "point_emi_any")],
                     M.pred = cmaq_full[, c("tmp", "wind")],
                     coords.Y = cmaq_aqs_matched[, c("x", "y")],
                     space.id.Y = cmaq_aqs_matched$space_id,
                     coords.pred = cmaq_full[, c("x", "y")],
                     space.id = cmaq_full$space_id,
                     time.id = cmaq_full$time_id,
                     spacetime.id = cmaq_full$spacetime_id,
                     include.additive.annual.resid = T,
                     include.multiplicative.annual.resid = T,
                     n.iter = 20,
                     verbose = T)



?modis_full

maia_pred <- grm_pred(grm.fit = maia_fit,
                      X.pred = modis_full$aod, 
                      L.pred = modis_full[, c("elevation", "forestcover",
                                              "hwy_length", "lim_hwy_length", 
                                              "local_rd_length", "point_emi_any")],
                      M.pred = modis_full[, c("tmp", "wind", "cmaq", "tempaod", 
                                              "windaod", "elevationaod")],
                      coords.Y = modis_aqs_matched[, c("x", "y")],
                      space.id.Y = modis_aqs_matched$space_id,
                      coords.pred = modis_full[, c("x", "y")],
                      space.id = modis_full$space_id, 
                      time.id = modis_full$time_id, 
                      spacetime.id = modis_full$spacetime_id,
                      include.additive.annual.resid = T,
                      include.multiplicative.annual.resid = T,
                      n.iter = 100,
                      verbose = T)
```

### Stage 3

``` r

ensemble_fit <- ensemble_spatial(grm.fit.cv.1 = ctm_fit_cv,
                                 grm.fit.cv.2 = maia_fit_cv,
                                 date.Y.1 = cmaq_aqs_matched$date,
                                 date.Y.2 = modis_aqs_matched$date,
                                 coords.Y.1 = cmaq_aqs_matched[, c("x", "y")],
                                 space.id.Y.1 = cmaq_aqs_matched$space_id,
                                 n.iter = 5000, 
                                 burn = 1000, 
                                 thin = 4,
                                 tau.a = 0.001, 
                                 tau.b = 0.001, 
                                 theta.tune = 0.2, 
                                 theta.a = 5, 
                                 theta.b = 0.05)
```

### Stage 4

``` r

weight_preds <- weight_pred(ensemble.fit = ensemble_fit,
                            coords.Y.1 = cmaq_aqs_matched[, c("x", "y")], 
                            space.id.Y.1 = cmaq_aqs_matched$space_id, 
                            coords.pred.1 = cmaq_full[, c("x", "y")],
                            space.id.pred.1 = cmaq_full$space_id)


results <- gap_fill(grm.pred.1 = ctm_pred,
                    grm.pred.2 = maia_pred,
                    date.pred.1 = cmaq_full$date,
                    date.pred.2 = modis_full$date, 
                    space.id = ctm_pred$space.id, 
                    weights = weight_preds)
```
