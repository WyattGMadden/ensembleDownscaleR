#' Make Predictions With Geostatistical Regression Model (GRM)
#'
#' This function makes predictions using a Bayesian Hierarchical Model (BHM) previously fit with grm()
#'
#' @inheritParams grm
#' @param grm.fit Fit object created with grm()
#' @param X.pred Standardized primary independent variable vector for all predictions (n_pred)
#' @param L.pred Standardized spatial covariate matrix for all predictions (n_pred, p1)
#' @param M.pred Standardized spatio-temporal covariate matrix for all predictions (n_pred, p2)
#' @param coords.Y Matrix of primary independent variable x, y coodinates (n_obs, 2)
#' @param space.id.Y Vector of space id numbers for primary variables, corresponding with coords.Y (n_obs)
#' @param coords.pred Matrix of prediction x, y coodinates (n_pred, 2)
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
                    X.pred, 
                    L.pred = NULL, 
                    M.pred = NULL, 
                    coords.Y, 
                    space.id.Y,
                    coords.pred, 
                    space.id, 
                    time.id, 
                    spacetime.id, 
                    include.additive.annual.resid = T,
                    include.multiplicative.annual.resid = T,
                    n.iter = 500,
                    verbose = TRUE,
                    in.sample = FALSE) {

    if (!in.sample) {

        ###Print some information
        cat("Preparing for Prediction\n")
  
        N <- length(X.pred) #Total AOD observations
        N.space <- max(space.id) #Total number of prediction cells
        N.time <- max(time.id) #Maximum number of time interval (weeks) in prediction
        N.time.obs <- length(unique(time.id)) #Number of observed time interval (weeks)
        N.spacetime <- max(spacetime.id) #Time points where spatial trends vary by (year)
        
        N.Lmax <- 0
        if (!is.null(L.pred)) {
            N.Lmax <- ncol(L.pred) #Number of spatial predictors to use
        }
        N.Mmax <- 0
        if (!is.null(M.pred)) {
            N.Mmax <- ncol(M.pred) #Number of spatial-temporal predictors to use
        }

        cov_kern <- grm.fit$cov_kern

        ############################
        ###standardize X, L and M###
        ############################
        X.pred <- scale(X.pred)

        if (!is.null(L.pred)) {
            L.pred <- as.matrix(L.pred)
            # scale unless variable is all one value 
            #(can happen in cross validation with binary variables)
            L.pred <- apply(X = L.pred, 
                            MARGIN = 2, 
                            FUN = function(x) ifelse(length(unique(x)) != 1,
                                                    scale(x),
                                                    x))
        }

        if (!is.null(M.pred)) {
            M.pred <- as.matrix(M.pred)
            # scale unless variable is all one value 
            #(can happen in cross validation with binary variables)
            M.pred <- apply(X = M.pred, 
                            MARGIN = 2, 
                            FUN = function(x) ifelse(length(unique(x)) <= 1,
                                                    scale(x),
                                                    x))
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
        dist_dat_Y <- unique(cbind(space.id.Y, coords.Y))
        dist_dat_Y <- dist_dat_Y[order(dist_dat_Y$space.id.Y), ]
        locations.Y <- dist_dat_Y[, c("x", "y")]

        dist_dat_pred <- unique(cbind(space.id, coords.pred))
        dist_dat_pred <- dist_dat_pred[order(dist_dat_pred$space.id), ]
        locations.pred <- dist_dat_pred[, c("x", "y")]

        N.mon <- nrow(locations.Y)
        N.cell <- nrow(locations.pred)

        
        ####Predict alpha and beta at grid cells
        XY <- rbind(locations.Y, locations.pred)
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
        if (include.additive.annual.resid) {

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
        if (include.multiplicative.annual.resid) { if (verbose == TRUE) {
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
    
        X.pred <- (X.pred - standardize.param[standardize.param$Type == "X", ]$Mean) / 
            standardize.param[standardize.param$Type == "X", ]$SD
    
        if (!is.null(L.pred)) {
            L.var <- as.character(standardize.param[standardize.param$Type == "L", ]$Name)

            for (l in L.var) {
                L.pred[, colnames(L.pred) == l] <- (L.pred[, colnames(L.pred) == l] - 
                                                   standardize.param[standardize.param$Name == l, ]$Mean) / 
                standardize.param[standardize.param$Name == l, ]$SD
            }
        }

        if (!is.null(M.pred)) {
            M.var <- as.character(standardize.param[standardize.param$Type == "M", ]$Name)

            for (m in M.var) {
                M.pred[, colnames(M.pred) == m] <- (M.pred[, colnames(M.pred) == m] - 
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
            fix.L <- L.pred %*% t(as.matrix(grm.fit$gamma[m, ]))
        }
        if (!is.null(grm.fit$delta)) {
            fix.M <- M.pred %*% t(as.matrix(grm.fit$delta[m, ]))
        }

        pred.mu <- intercept + slope * X.pred
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
