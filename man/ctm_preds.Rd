% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ctm_preds}
\alias{ctm_preds}
\title{Chemical Transport Model Data For Atlanta Metro Area, June 2004}
\format{
\subsection{\code{ctm_preds}}{

A data frame with 72,000 rows and 15 columns:
\describe{
\item{time_id}{Week identifier}
\item{space_id}{Location identifier}
\item{spacetime_id}{Year identifier}
\item{ctm}{Chemical transport model value}
\item{elevation}{Elevation spatial covariate (US Geological Survey)}
\item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
\item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
\item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
\item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
\item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
\item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
\item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
\item{date}{Date of measurement}
\item{x}{Location x-coordinate}
\item{y}{Location y-coordinate}
...
}
}
}
\source{
\url{https://www.nature.com/articles/jes201390}
}
\usage{
ctm_preds
}
\description{
Daily CMT output, and covariates of interest in Atlanta Metro Area for June 2004.
Data included for every 12km x 12km grid square in study area.
Each row of the data is a unique location/day.
}
\keyword{datasets}
