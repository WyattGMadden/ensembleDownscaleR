% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_cv.R
\name{create_cv_from_previous}
\alias{create_cv_from_previous}
\title{Create Cross Validation ID's For New Dataset Based On Previously Created Cross Validation ID's}
\usage{
create_cv_from_previous(
  previous.cv.object,
  time.id,
  space.id,
  spacetime.id,
  coords
)
}
\arguments{
\item{previous.cv.object}{The cross-validation object created from the original dataset, used to determine the cross-validation assignment for the new dataset}

\item{time.id}{Temporal location ID vector (n)}

\item{space.id}{Spatial location ID vector (n)}

\item{spacetime.id}{ID vector of time points where spatial trends vary (n)}

\item{coords}{Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (n, 2)}
}
\value{
A named list containing a vector of cross validation ID's, a matrix of which observations to drop for each fold if the cv type is "spatial_buffered", and inputted objects.
}
\description{
This function creates creates a cross-validation assignment for a new dataset, based off a previously created cross-validation assignment
}
