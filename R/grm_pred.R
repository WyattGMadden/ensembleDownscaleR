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

grm_pred <- function(
    grm.fit,
    X, 
    L = NULL, 
    M = NULL, 
    coords,
    space.id, 
    time.id, 
    spacetime.id, 
    n.iter = 500,
    verbose = TRUE,
    in.sample = FALSE
    ) {

    ###################################
    ###         ARG CHECKS          ###
    ###################################
    
    # Check if grm.fit is from grm()
    needed_in_grmfit <- c("others", "alpha.space", "beta.space", 
                          "alpha.time", "standardize.param", 
                          "locations", "cov_kern")
    if (!is.list(grm.fit) ||
        !all(needed_in_grmfit %in% names(grm.fit))) {
        stop("'grm.fit' must be a list output from grm(), containing at least the elements:",
             " 'others', 'alpha.space', 'beta.space', 'alpha.time', 'standardize.param', 'locations', 'cov_kern'.")
    }

    # X check
    if (!is.numeric(X) || !is.vector(X)) {
        stop("'X' must be a numeric vector.")
    }
    # L, M checks
    if (!is.null(L)) {
        if (!is.matrix(L) & !is.data.frame(L)) {
            stop("'L' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(L) != length(X)) {
            stop("Number of rows in 'L' must match length of 'Y'.")
        }
    }
    if (!is.null(M)) {
        if (!is.matrix(M) & !is.data.frame(M)) {
            stop("'M' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(M) != length(X)) {
            stop("Number of rows in 'M' must match length of 'Y'.")
        }
    }
    # coords checks
    if (!is.matrix(coords) & !is.data.frame(coords)) {
        stop("'coords' must be a matrix.")
    }
    if (nrow(coords) != length(X)) {
        stop("Number of rows in 'coords' must match length of 'X'.")
    }
    if (!all(colnames(coords) == c("x", "y"))) {
        stop("'coords' must have colnames 'x' and 'y'.")
    }
    # space.id/time.id/spacetime.id checks
    if (length(space.id) != length(X)) {
        stop("'space.id' must have the same length as 'X'.")
    }
    if (length(time.id) != length(X)) {
        stop("'time.id' must have the same length as 'X'.")
    }
    if (length(spacetime.id) != length(X)) {
        stop("'spacetime.id' must have the same length as 'X'.")
    }


    # n.iter checks
    if (!is.numeric(n.iter) || n.iter < 1) {
        stop("'n.iter' must be a positive numeric scalar.")
    }
    # in.sample checks
    if (!is.logical(in.sample) || length(in.sample) != 1) {
        stop("'in.sample' must be a single TRUE/FALSE value.")
    }
    # check if n.iter <= rows in grm.fit$others
    if (!in.sample && nrow(grm.fit$others) < n.iter) {
        stop("'n.iter' must be less than or equal to the number of MCMC samples in 'grm.fit$others'.")
    }
    if (!in.sample & (nrow(grm.fit$others) < n.iter)) {
        stop("n.iter must be less than or equal to the number of iterations in the grm.fit object")
    }
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
        

        if (is.null(grm.fit$nngp.info)) {

            ###Create prediction distance matrix
            locations.Y <- grm.fit$locations[, c("x", "y")] 

            dist_dat_pred <- unique(cbind(space.id, coords))
            dist_dat_pred <- dist_dat_pred[order(dist_dat_pred$space.id), ]
            locations.pred <- dist_dat_pred[, c("x", "y")]

            N.mon <- nrow(locations.Y)
            N.cell <- nrow(locations.pred)
            
            ####Predict alpha and beta at grid cells
            XY <- rbind(locations.Y, locations.pred)
            D22 <- as.matrix(stats::dist(locations.Y, diag = TRUE, upper = TRUE))
            D12 <- as.matrix(stats::dist(XY, diag = TRUE, upper = TRUE))[c(1:N.mon), -c(1:N.mon)]
            
            alpha_space_pred <- data.frame(expand.grid(1:N.space, 
                                                       1:N.spacetime))
            beta_space_pred <- data.frame(expand.grid(1:N.space, 
                                                     1:N.spacetime))

            names(alpha_space_pred) <- c("space.id", "spacetime.id")
            names(beta_space_pred) <- c("space.id", "spacetime.id")

            alpha_space_pred[paste0("Sample", 1:n.iter)] <- 0
            beta_space_pred[paste0("Sample", 1:n.iter)] <- 0





            
            #For alpha's
            if (grm.fit$include.additive.spatial.effect) {

              if (verbose == TRUE) {

                cat("Imputing Spatial Alphas\n") 

              }

              for (m in 1:n.iter) {
                  tau.m <- grm.fit$others$tau.alpha[m]
                  theta.m <- grm.fit$others$theta.alpha[m]
                  Sigma11.m <- tau.m
                  Sigma12.m <- tau.m * cov_kern(distance = D12, 
                                               theta = theta.m)

                  if (is.null(grm.fit$discrete.theta.alpha.info)) {
                      Sigma22.m <- tau.m * cov_kern(distance = D22, 
                                                   theta = theta.m)
                      InvSigma22.m <- solve(Sigma22.m)
                  } else {
                      which.theta.alpha <- grm.fit$discrete.theta.alpha.info$which.theta.alpha.discrete[m]
                      InvSigma22.m <- (1 / tau.m) * grm.fit$discrete.theta.alpha.info$kernals.inv.alpha[[which.theta.alpha]]
                  }
            
                  for (j in 1:N.spacetime) {
                      alpha.m <- grm.fit$alpha.space[grm.fit$alpha.space$spacetime.id == j, 
                                                    paste0("Sample", m)]
                      alpha.mu.m <- t(Sigma12.m) %*% InvSigma22.m %*% alpha.m
                      alpha.cov.m <- Sigma11.m - diag(t(Sigma12.m) %*% InvSigma22.m %*% Sigma12.m)
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
            if (grm.fit$include.multiplicative.spatial.effect) { 

                if (verbose == TRUE) {
                    cat("Imputing Spatial Betas\n") 
                }

                for (m in 1:n.iter) {
                    tau.m <- grm.fit$others$tau.beta[m]
                    theta.m <- grm.fit$others$theta.beta[m] 
                    Sigma11.m <- tau.m
                    Sigma12.m <- tau.m * cov_kern(distance = D12, 
                                                 theta = theta.m)

                    if (is.null(grm.fit$discrete.theta.beta.info)) {
                        Sigma22.m <- tau.m * cov_kern(distance = D22, 
                                                     theta = theta.m)
                        InvSigma22.m <- solve(Sigma22.m)
                    } else {
                      which.theta.beta <- grm.fit$discrete.theta.beta.info$which.theta.beta.discrete[m]
                      InvSigma22.m <- (1 / tau.m) * grm.fit$discrete.theta.beta.info$kernals.inv.beta[[which.theta.beta]]
                  }
              
                    for (j in 1:N.spacetime) {
                        beta.m <- grm.fit$beta.space[grm.fit$beta.space$spacetime.id == j, 
                                            paste0("Sample", m)]
                        beta.mu.m <- t(Sigma12.m) %*% InvSigma22.m %*% beta.m
                        beta.cov.m <- Sigma11.m - diag(t(Sigma12.m) %*% 
                                                      InvSigma22.m %*% 
                                                      Sigma12.m)
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

        } else if (!is.null(grm.fit$nngp.info)) {




            ###Create prediction distance matrix
            locations.Y <- grm.fit$locations[, c("x", "y")] 

            dist_dat_pred <- unique(cbind(space.id, coords))
            dist_dat_pred <- dist_dat_pred[order(dist_dat_pred$space.id), ]
            locations.pred <- dist_dat_pred[, c("x", "y")]

            N.mon <- nrow(locations.Y)
            N.cell <- nrow(locations.pred)

            nngp_info <- grm.fit$nngp.info
            ordered_coords_Y <- nngp_info$ordered.coords
            neighbors_pred <- get_neighbors_ref(ordered_coords = nngp_info$ordered.coords, 
                                                locations.pred, 
                                                nngp_info$number.neighbors)

            dist_mats_pred <- get_dist_matrices_ref(ordered_coords = nngp_info$ordered.coords, 
                                                    coords_pred = locations.pred,
                                                    neighbors = neighbors_pred)
            
            
            alpha_space_pred  <- data.frame(expand.grid(1:N.space, 
                                                       1:N.spacetime))
            beta_space_pred <- data.frame(expand.grid(1:N.space, 
                                                     1:N.spacetime))

            names(alpha_space_pred) <- c("space.id", "spacetime.id")
            names(beta_space_pred) <- c("space.id", "spacetime.id")

            alpha_space_pred[paste0("Sample", 1:n.iter)] <- 0
            beta_space_pred[paste0("Sample", 1:n.iter)] <- 0
            
            #For alpha's
            if (grm.fit$include.additive.spatial.effect) {

              if (verbose == TRUE) {

                cat("Imputing Spatial Alphas\n") 

              }
                for (m in 1:n.iter) {
                    tau.m <- grm.fit$others$tau.alpha[m]
                    theta.m <- grm.fit$others$theta.alpha[m]
                    alpha.m <- grm.fit$alpha.space[, paste0("Sample", m)]
                    alpha.space.id <- grm.fit$alpha.space$space.id
                    alpha.spacetime.id <- grm.fit$alpha.space$spacetime.id

                    for (j in 1:N.spacetime) {

                        alpha.m.j <- alpha.m[alpha.spacetime.id == j]
                        alpha.m.j.ord <- alpha.m.j[nngp_info$coord.ordering]
                        alpha.j.space.id.pred <- alpha_space_pred$space.id[alpha_space_pred$spacetime.id == j]
                        alpha.m.j.pred <- rep(0, length(alpha.j.space.id.pred))

                        for (i in 1:length(alpha.j.space.id.pred)) {
                            sigma.m.j.i <- tau.m * cov_kern(distance = dist_mats_pred[[i]], 
                                                            theta = theta.m)
                            inv.sigma.m.j.i <- solve(sigma.m.j.i[-1, -1])
                            alpha.mu.m.j.i <- sigma.m.j.i[1, -1] %*% 
                                inv.sigma.m.j.i %*% 
                                alpha.m.j.ord[neighbors_pred[[i]]]

                            alpha.cov.m.j.i <- sigma.m.j.i[1, 1] - 
                                sigma.m.j.i[1, -1] %*% 
                                inv.sigma.m.j.i %*% 
                                sigma.m.j.i[-1, 1]

                            alpha.m.j.pred[i] <- stats::rnorm(1, 
                                                              alpha.mu.m.j.i, 
                                                              sqrt(alpha.cov.m.j.i))

                        }
                        alpha_space_pred[alpha_space_pred$spacetime.id == j,
                                         paste0("Sample", m)] <- alpha.m.j.pred
                    } #End of locations

                    if (verbose == TRUE) {
                        if (m %% (n.iter / 10) == 0) {
                            cat(paste("     Iteration", m, "of", n.iter, "\n"))
                        }
                    }
                }
            }
             
            #For betas's
            if (grm.fit$include.multiplicative.spatial.effect) { 

                if (verbose == TRUE) {
                    cat("Imputing Spatial Betas\n") 
                }

                for (m in 1:n.iter) {

                    tau.m <- grm.fit$others$tau.beta[m]
                    theta.m <- grm.fit$others$theta.beta[m]
                    beta.m <- grm.fit$beta.space[, paste0("Sample", m)]
                    beta.space.id <- grm.fit$beta.space$space.id
                    beta.spacetime.id <- grm.fit$beta.space$spacetime.id

                    for (j in 1:N.spacetime) {

                        beta.m.j <- beta.m[beta.spacetime.id == j]
                        beta.m.j.ord <- beta.m.j[nngp_info$coord.ordering]
                        beta.j.space.id.pred <- unique(beta_space_pred$space.id[beta_space_pred$spacetime.id == j])
                        beta.m.j.pred <- rep(0, length(beta.j.space.id.pred))

                        for (i in 1:length(beta.j.space.id.pred)) {
                            sigma.m.j.i <- tau.m * cov_kern(distance = dist_mats_pred[[i]], 
                                                            theta = theta.m)
                            inv.sigma.m.j.i <- solve(sigma.m.j.i[-1, -1])
                            beta.mu.m.j.i <- sigma.m.j.i[1, -1] %*% 
                                inv.sigma.m.j.i %*% 
                                beta.m.j.ord[neighbors_pred[[i]]]

                            beta.cov.m.j.i <- sigma.m.j.i[1, 1] - 
                                sigma.m.j.i[1, -1] %*% 
                                inv.sigma.m.j.i %*% 
                                sigma.m.j.i[-1, 1]

                            beta.m.j.pred[i] <- stats::rnorm(1, 
                                                              beta.mu.m.j.i, 
                                                              sqrt(beta.cov.m.j.i))

                        }
                        beta_space_pred[beta_space_pred$spacetime.id == j,
                                         paste0("Sample", m)] <- beta.m.j.pred
                    } #End of locations
                    if (verbose == TRUE) {
                        if (m %% (n.iter / 10) == 0) {
                            cat(paste("     Iteration", m, "of", n.iter), "\n")
                        }
                    }

                }#End of iteration

            }

        }


    } else if (in.sample) {

        #in sample uses random effects from the model fit
        alpha_space_pred <- grm.fit$alpha.space
        beta_space_pred <- grm.fit$beta.space


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
