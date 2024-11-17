
#' Generate weights from ensemble spatial fit
#'
#' This function generates weights using the output from the ensemble_spatial function
#'
#' @inheritParams grm
#'
#' @param ensemble.fit Fit from ensemble_spatial() function
#' @param coords.Y.1 Coordinate matrix for first primary variable
#' @param space.id.Y.1 Space id for first primary variable 
#' @param coords.pred.1 Coordinate matrix for first variable predictions
#' @param space.id.pred.1 Space id for first primary variable predictions
#'
#' @return A matrix containing weights
#'
#' @examples
#' # weight_preds()
#' 
#' 
#' @export

weight_pred <- function(ensemble.fit, 
                       coords.Y.1,
                       space.id.Y.1,
                       coords.pred.1,
                       space.id.pred.1,
                       verbose = TRUE) {

    q <- ensemble.fit$q
    theta <- ensemble.fit$other$theta
    tau2 <- ensemble.fit$other$tau2

    #Create unique location matrices
    locations.Y <- unique(cbind(coords.Y.1, space.id.Y.1))
    locations.Y <- locations.Y[order(locations.Y$space.id.Y), ]
    locations.Y <- locations.Y[, c("x", "y")]

    locations.pred <- unique(cbind(coords.pred.1, space.id.pred.1))
    locations.pred <- locations.pred[order(locations.pred$space.id.pred.1), ]
    locations.pred <- locations.pred[, c("x", "y")]

    n.iter <- length(theta)
  
  
    #Some data processing
    q.mat <- as.matrix(q[, -1])
  
    N.mon <- nrow(locations.Y)
    N.cell <- nrow(locations.pred)
  
    D22 <- as.matrix(stats::dist(locations.Y, 
                                diag = TRUE, 
                                upper = TRUE))
    D12 <- calculate_distances(locations.Y, locations.pred)
  
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
        q.cov.m <- Sigma11.m - diag(t(Sigma12.m) 
                                   %*% InvSigma22.m 
                                   %*% Sigma12.m)
        q.m.post <- stats::rnorm(N.cell, q.mu.m, sqrt(q.cov.m))
    
        q.pred[, m] <- q.m.post
    }
  
    return(q.pred)
}
  
