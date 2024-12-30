#' Make Predictions With Geostatistical Regression Model (GRM)
#'
#' This function makes predictions using a Bayesian Hierarchical Model (BHM) previously fit with grm()
#'
#' @inheritParams grm
#' @param grm.fit Fit object created with grm()
#' @param X Standardized primary independent variable vector for all predictions (n_pred)
#' @param L Standardized spatial covariate matrix for all predictions (n_pred, p1)
#' @param M Standardized spatio-temporal covariate matrix for all predictions (n_pred, p2)
#' @param coords Matrix of prediction x, y coodinates (n_pred, 2)
#' @param n.iter Number of iterations used in predictions. Must be <= than post-thined and burned iterations from grm.fit
#' @param in.sample False if predictions are being made at locations without observations in grm.fit (default is False)
#'
#' @return A data frame containing grm predictions
#'
#' @examples
#' # grm_pred()
#' 
#' 
#' @export

grm_pred <- function(grm.fit,
                    X, 
                    L = NULL, 
                    M = NULL, 
                    coords,
                    space.id, 
                    time.id, 
                    spacetime.id, 
                    include.additive.spatial.effect = T,
                    include.multiplicative.spatial.effect = T,
                    n.iter = 500,
                    verbose = TRUE,
                    in.sample = FALSE) {

    # assertions
    if (!in.sample & (nrow(grm.fit$others) < n.iter)) {
        stop("n.iter must be less than or equal to the number of iterations in the grm.fit object")
    }

    if (!in.sample) {

        ###Print some information
        cat("Preparing for Prediction\n")
  
        N <- length(X) #Total AOD observations
        N.space <- length(unique(space.id)) #Total number of prediction cells
        N.time <- max(time.id) #Maximum number of time interval (weeks) in prediction
        N.time.obs <- length(unique(time.id)) #Number of observed time interval (weeks)
        N.spacetime <- max(spacetime.id) #Time points where spatial trends vary by (year)
        
        N.Lmax <- 0
        if (!is.null(L)) {
            N.Lmax <- ncol(L) #Number of spatial predictors to use
        }
        N.Mmax <- 0
        if (!is.null(M)) {
            N.Mmax <- ncol(M) #Number of spatial-temporal predictors to use
        }

        cov_kern <- grm.fit$cov_kern

        ############################
        ###standardize X, L and M###
        ############################


        standardize.param <- grm.fit$standardize.param

        X <- (X - standardize.param[standardize.param$Type == "X", ]$Mean) / 
            standardize.param[standardize.param$Type == "X", ]$SD

        if (!is.null(L)) {

            L <- as.matrix(L)
            L.var <- as.character(standardize.param[standardize.param$Type == "L", ]$Name)

            for (l in L.var) {
                L[, colnames(L) == l] <- (L[, colnames(L) == l] - 
                                                   standardize.param[standardize.param$Name == l, ]$Mean) / 
                standardize.param[standardize.param$Name == l, ]$SD
            }
        }

        if (!is.null(M)) {

            M <- as.matrix(M)
            M.var <- as.character(standardize.param[standardize.param$Type == "M", ]$Name)

            for (m in M.var) {
                M[, colnames(M) == m] <- (M[, colnames(M) == m] - 
                                                   standardize.param[standardize.param$Name == m, ]$Mean) / 
                standardize.param[standardize.param$Name == m, ]$SD
            }

        }



        
        if (verbose == TRUE) {
            cat("####################################\n")
            cat("######## Preparing for MCMC ########\n")
            cat("####################################\n")
            cat(paste0("     Total number of prediction time points: ", 
                         N.time.obs, 
                         " (out of ", 
                         N.time,")\n") )
            cat(paste0("     Total number of prediction cells: ", 
                         N.space,
                         "\n"))
            cat(paste0("     Total number of spatial covariates: ", 
                         N.Lmax,
                         "\n"))
            cat(paste0("     Total number of spatial-temporal covariates: ", 
                         N.Mmax,
                         "\n"))
        }
        

        ###Create prediction distance matrix
        locations.Y <- grm.fit$locations[, c("x", "y")] 

        dist_dat_pred <- unique(cbind(space.id, coords))
        dist_dat_pred <- dist_dat_pred[order(dist_dat_pred$space.id), ]
        locations.pred <- dist_dat_pred[, c("x", "y")]

        N.mon <- nrow(locations.Y)
        N.cell <- nrow(locations.pred)

        
        ####Predict alpha and beta at grid cells
        D22 <- as.matrix(stats::dist(locations.Y, diag = TRUE, upper = TRUE))
        D12 <- calculate_distances(locations.Y, locations.pred)


        alpha_space_pred <- data.frame(expand.grid(1:N.space, 
                                                   1:N.spacetime))
        beta_space_pred <- data.frame(expand.grid(1:N.space, 
                                                 1:N.spacetime))

        names(alpha_space_pred) <- c("space.id", "spacetime.id")
        names(beta_space_pred) <- c("space.id", "spacetime.id")

        alpha_space_pred[paste0("Sample", 1:n.iter)] <- 0
        beta_space_pred[paste0("Sample", 1:n.iter)] <- 0





        
        #For alpha's
        if (include.additive.spatial.effect) {

          if (verbose == TRUE) {

            cat("Imputing Spatial Alphas\n") 

          }

          for (m in 1:n.iter) {
              tau.m <- grm.fit$others$tau.alpha[m]
              theta.m <- grm.fit$others$theta.alpha[m]
              Sigma11.m <- tau.m
              Sigma12.m <- tau.m * cov_kern(distance = D12, 
                                           theta = theta.m)

              Sigma22.m <- tau.m * cov_kern(distance = D22, 
                                            theta = theta.m)
              InvSigma22.m <- solve(Sigma22.m)
        
              for (j in 1:N.spacetime) {
                  alpha.m <- grm.fit$alpha.space[grm.fit$alpha.space$spacetime.id == j, 
                                                paste0("Sample", m)]
                  alpha.mu.m <- t(Sigma12.m) %*% InvSigma22.m %*% alpha.m

                  # In below for loop we are 
                  # calculating diag_values = diag(t(Sigma12.m) %*% InvSigma22.m %*% Sigma12.m)
                  # in a memory efficient fashion 
                  diag_values <- numeric(ncol(Sigma12.m))
                  for (i in 1:ncol(Sigma12.m)) {
                      col_vector <- Sigma12.m[, i, drop = FALSE]
  
                      temp_product <- t(col_vector) %*% InvSigma22.m %*% col_vector
  
                      diag_values[i] <- temp_product
                  }

                  alpha.cov.m <- Sigma11.m - diag_values
                  alpha.m.post <- stats::rnorm(N.cell, alpha.mu.m, sqrt(alpha.cov.m))
                  alpha_space_pred[alpha_space_pred$spacetime.id == j, 
                                   paste0("Sample",m)] = alpha.m.post
              } #End of locations
            
              if (verbose == TRUE) {
                  if (m %% (n.iter / 10) == 0) {
                      cat(paste("     Iteration", m, "of", n.iter, "\n"))
                  }
              }
          }
        }
         
        #For betas's
        if (include.multiplicative.spatial.effect) { if (verbose == TRUE) {
                cat("Imputing Spatial Betas\n") 
            }

            for (m in 1:n.iter) {
                tau.m <- grm.fit$others$tau.beta[m]
                theta.m <- grm.fit$others$theta.beta[m] 
                Sigma11.m <- tau.m
                Sigma12.m <- tau.m * cov_kern(distance = D12, 
                                             theta = theta.m)

                Sigma22.m <- tau.m * cov_kern(distance = D22, 
                                              theta = theta.m)
                InvSigma22.m <- solve(Sigma22.m)
          
                for (j in 1:N.spacetime) {
                    beta.m <- grm.fit$beta.space[grm.fit$beta.space$spacetime.id == j, 
                                        paste0("Sample", m)]
                    beta.mu.m <- t(Sigma12.m) %*% InvSigma22.m %*% beta.m

                    # In below for loop we are 
                    # calculating diag_values = diag(t(Sigma12.m) %*% InvSigma22.m %*% Sigma12.m)
                    # in a memory efficient fasion 
                    diag_values <- numeric(ncol(Sigma12.m))
                    for (i in 1:ncol(Sigma12.m)) {
                        col_vector <- Sigma12.m[, i, drop = FALSE]
  
                        temp_product <- t(col_vector) %*% InvSigma22.m %*% col_vector
  
                        diag_values[i] <- temp_product
                    }

                    beta.cov.m <- Sigma11.m - diag_values
                    beta.m.post <- stats::rnorm(N.cell, beta.mu.m, sqrt(beta.cov.m))
                    beta_space_pred[beta_space_pred$spacetime.id == j, 
                                    paste0("Sample", m)] = beta.m.post
                } #End of locations
                if (verbose == TRUE) {
                    if (m %% (n.iter / 10) == 0) {
                        cat(paste("     Iteration", m, "of", n.iter), "\n")
                    }
                }
        
            }#End of iterations
        }



    } else if (in.sample) {

        #in sample uses random effects from the model fit
        alpha_space_pred <- grm.fit$alpha.space
        beta_space_pred <- grm.fit$beta.space

        #standardize based on the model fit
        standardize.param <- grm.fit$standardize.param
    
        X <- (X - standardize.param[standardize.param$Type == "X", ]$Mean) / 
            standardize.param[standardize.param$Type == "X", ]$SD
    
        if (!is.null(L)) {
            L.var <- as.character(standardize.param[standardize.param$Type == "L", ]$Name)

            for (l in L.var) {
                L[, colnames(L) == l] <- (L[, colnames(L) == l] - 
                                                   standardize.param[standardize.param$Name == l, ]$Mean) / 
                standardize.param[standardize.param$Name == l, ]$SD
            }
        }

        if (!is.null(M)) {
            M.var <- as.character(standardize.param[standardize.param$Type == "M", ]$Name)

            for (m in M.var) {
                M[, colnames(M) == m] <- (M[, colnames(M) == m] - 
                                                   standardize.param[standardize.param$Name == m, ]$Mean) / 
                standardize.param[standardize.param$Name == m, ]$SD
            }
        }


    } else {

        stop("You must specify either in.sample = TRUE or in.sample = FALSE")

    }
  
    ####Make Predictions
    results <- data.frame(time.id, space.id, spacetime.id)
    results$estimate <- 0
    results$sd <- 0
    
    id.temp.pred <- paste0(alpha_space_pred$space.id, 
                          "_", 
                          alpha_space_pred$spacetime.id)

    id.temp <- paste0(space.id, "_", spacetime.id)
    
    if (!is.null(grm.fit$delta)) {
        delta <- as.matrix(grm.fit$delta)
    }
    if (!is.null(grm.fit$gamma)) {
        gamma <- as.matrix(grm.fit$gamma)
    }

    
    
    for (m in 1:n.iter) {
        intercept <- grm.fit$others$alpha0[m] + 
            grm.fit$alpha.time[time.id, m + 1] + 
            alpha_space_pred[match(id.temp, id.temp.pred), m + 2]
        slope <- grm.fit$others$beta0[m] + 
            grm.fit$beta.time[time.id, m + 1] + 
            beta_space_pred[match(id.temp, id.temp.pred), m + 2]

        if (!is.null(grm.fit$gamma)) {
            fix.L <- L %*% t(as.matrix(grm.fit$gamma[m, ]))
        }
        if (!is.null(grm.fit$delta)) {
            fix.M <- M %*% t(as.matrix(grm.fit$delta[m, ]))
        }

        pred.mu <- intercept + slope * X
        if (!is.null(grm.fit$gamma)) {
            pred.mu <- pred.mu + as.vector(fix.L)
        }
        if (!is.null(grm.fit$delta)) {
            pred.mu <- pred.mu + as.vector(fix.M)
        }

        pred.mu <- pred.mu + stats::rnorm(length(pred.mu), 
                                         0, 
                                         sqrt(grm.fit$others$sigma2[m])) 
        results$estimate <- results$estimate + pred.mu / n.iter
        results$sd <- results$sd + pred.mu^2 / n.iter
    }

    results$sd <- sqrt((results$sd - results$estimate^2))


    return(results)

}
