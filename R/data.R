#' Community Multiscale Air Quality (CMAQ) Chemical Transport Model (CMT) and USEPA Air Quality System (AQS) PM2.5 Data For Atlanta Metro Area, 2003-2005
#'
#' Daily AQS PM2.5 readings matched with CMAQ CMT output and covariates of interest in Atlanta Metro Area from 2003-2005. 
#' Each row of the data is a unique AQS monitor location/day. 
#' Data included only if CMAQ CMT output is available for that day and monitoring location.
#' (ie. if the monitor is within a 12km x 12km grid square for which CMT data is available).
#'
#' @format ## `cmaq_aqs_matched`
#' A data frame with 11,404 rows and 16 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{pm25}{PM2.5 reading}
#'   \item{ctm}{Chemical transport model value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"cmaq_aqs_matched"


#' Moderate Resolution Imaging Spectro-radiometer (MODIS) Aerosol Optical Depth (AOD) and USEPA Air Quality System (AQS) PM2.5 Data For Atlanta Metro Area, 2003-2005
#'
#' Daily AQS PM2.5 readings matched with MODIS AOD and covariates of interest in Atlanta Metro Area from 2003-2005. 
#' Each row of the data is a unique AQS monitor location/day. 
#' Data included only if AOD output is available for that day and monitoring location.
#' (ie. if the monitor is within a 12km x 12km grid square for which AOD data is available).
#'
#' @format ## `modis_aqs_matched`
#' A data frame with 6531 rows and 20 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{pm25}{PM2.5 reading}
#'   \item{aod}{Aerosol optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{cmaq}{CMT value spatio-temporal covariate (USEPA Models-3/Community Multiscale Air Quality)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://pubmed.ncbi.nlm.nih.gov/31465992/>
"modis_aqs_matched"

#' Community Multiscale Air Quality (CMAQ) Chemical Transport Model (CMT) Data For Atlanta Metro Area, June 2004
#'
#' Daily CMAQ CMT output, and covariates of interest in Atlanta Metro Area for June 2004. 
#' Data included for every 12km x 12km grid square in study area. 
#' Each row of the data is a unique location/day. 
#'
#' @format ## `cmaq_full`
#' A data frame with 72,000 rows and 15 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{ctm}{Chemical transport model value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"cmaq_full"


#' Moderate Resolution Imaging Spectro-radiometer (MODIS) Aerosol Optical Depth (AOD) For Atlanta Metro Area, June 2004
#'
#' Daily MODIS AOD, and covariates of interest in Atlanta Metro Area for June 2004.
#' Data included for every 12km x 12km grid square in study area for which AOD data is available. 
#' Each row of the data is a unique location/day. 
#'
#' @format ## `modis_full`
#' A data frame with 6531 rows and 20 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{aod}{Aerosol optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{cmaq}{CMT value spatio-temporal covariate (USEPA Models-3/Community Multiscale Air Quality)}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   ...
#' }
#' @source <https://pubmed.ncbi.nlm.nih.gov/31465992/>
"modis_full"

