
#' Generate weights from ensemble spatial fit
#'
#' This function generates weights using the output from the ensemble_spatial function
#'
#' @inheritParams grm
#'
#' @param ensemble.fit Fit from ensemble_spatial() function
#' @param coords Coordinate matrix for first variable predictions
#' @param space.id Space id for first primary variable predictions
#'
#' @return A list with:
#' \itemize{
#'   \item{\code{q}}{ : A matrix of imputed log-odds (weights) for each location, of dimension \code{(n.cell x n.iter)}.}
#'   \item{\code{locations}}{ : The \code{coords} and \code{space.id} used for prediction, in ordered form.}
#' }
#'
#' @examples
#' # weight_preds()
#' 
#' 
#' @export

weight_pred <- function(
    ensemble.fit, 
    coords,
    space.id,
    verbose = TRUE
    ) {

    ###################################
    ###         ARG CHECKS          ###
    ###################################

    # Single combined check that ensures 'ensemble.fit' is from ensemble_spatial()
    needed_in_ensemble <- c("q", "other", "locations")
    if (!is.list(ensemble.fit) ||
        !all(needed_in_ensemble %in% names(ensemble.fit)) ||
        !is.data.frame(ensemble.fit$q) ||
        !is.data.frame(ensemble.fit$other) ||
        !is.data.frame(ensemble.fit$locations)) {
        stop(
            "'ensemble.fit' must be a list output from ensemble_spatial(), containing at least:\n",
            " - 'q': a data frame of posterior samples for log-odds q,\n",
            " - 'other': a data frame with posterior samples for tau2, theta, dev, etc.,\n",
            " - 'locations': a data frame of the reference x,y coordinates."
        )
    }

    # Check coords
    if (!is.matrix(coords)) {
        stop("'coords' must be a matrix.")
    }
    if (ncol(coords) != 2 || !all(colnames(coords) %in% c("x","y"))) {
        stop("'coords' must have two columns named 'x' and 'y'.")
    }

    # Check space.id
    if (length(space.id) != nrow(coords)) {
        stop("'space.id' length must match the number of rows in 'coords'.")
    }

    q <- ensemble.fit$q
    theta <- ensemble.fit$other$theta
    tau2 <- ensemble.fit$other$tau2

    #Create unique location matrices
    locations.Y <- ensemble.fit$locations[, c("x", "y")]

    locations.pred <- unique(cbind(coords, space.id))
    locations.pred <- locations.pred[order(locations.pred$space.id), ]

    n.iter <- length(theta)
  
  
    #Some data processing
    q.mat <- as.matrix(q[, -1])
  
    N.mon <- nrow(locations.Y)
    N.cell <- nrow(locations.pred)
  
    D22 <- as.matrix(stats::dist(locations.Y, 
                                diag = TRUE, 
                                upper = TRUE))
    D12 <- calculate_distances(locations.Y, locations.pred[, c("x", "y")])
  
    q.pred <- matrix(NA, N.cell, n.iter)
  
    if (verbose == TRUE) {
        print("Imputing Spatial Weights") 
    }
  
    for (m in 1:n.iter) {
   
        if (verbose == TRUE & m%%100 == 0) {
            print(paste0("Imputing ", m, " of ", n.iter))
        } 
    
        q.m <- q.mat[, m]
        tau2.m <- tau2[m]
        theta.m <- theta[m]
        Sigma11.m <- tau2.m
        Sigma12.m <- tau2.m * exp(-1 / theta.m * D12)
        Sigma22.m <- tau2.m * exp(-1 / theta.m * D22)
        InvSigma22.m <- solve(Sigma22.m)
    
        q.mu.m <- t(Sigma12.m) %*% InvSigma22.m %*% q.m
        # In below for loop we are 
        # calculating diag_values = diag(t(Sigma12.m) %*% InvSigma22.m %*% Sigma12.m)
        # in a memory efficient fashion 
        diag_values <- numeric(ncol(Sigma12.m))
        for (i in 1:ncol(Sigma12.m)) {
            col_vector <- Sigma12.m[, i, drop = FALSE]
  
            temp_product <- t(col_vector) %*% InvSigma22.m %*% col_vector

            diag_values[i] <- temp_product
        }
        q.cov.m <- Sigma11.m - diag_values
        q.m.post <- stats::rnorm(N.cell, q.mu.m, sqrt(q.cov.m))
    
        q.pred[, m] <- q.m.post
    }
    # clean
    locations <- locations.pred[, c("x", "y", "space.id")]
  
    return(list(q = q.pred,
                locations = locations))
}
