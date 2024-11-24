#' Fit the Spatial Ensemble model
#'
#' This function fits the spatial ensemble, using output from the grm_cv function applied to multiple data sets
#'
#' @inheritParams grm
#
#' @param grm.fit.cv.1 Output from grm_cv function for first primary variable
#' @param grm.fit.cv.2 Output from grm_cv function for second primary variable
#' @param tau.a First tau prior hyperparameter
#' @param tau.b Second tau prior hyperparameter
#' @param theta.tune Theta Metropolis-Hastings proposal tuning parameter
#' @param theta.a First theta prior hyperparameter
#' @param theta.b Second theta prior hyperparameter
#'
#' @return A data frame containing cross validation predictions
#'
#' @examples
#' # ensemble_spatial()
#' 
#' 
#' @export


ensemble_spatial <- function(grm.fit.cv.1,
                             grm.fit.cv.2,
                             n.iter = 25000, 
                             burn = 5000, 
                             thin = 4,
                             tau.a = 0.001, 
                             tau.b = 0.001, 
                             theta.tune = 0.2, 
                             theta.a = 5, 
                             theta.b = 0.05) {
  

    space.id.Y.1 <- grm.fit.cv.1$space.id
    coords.Y.1 <- grm.fit.cv.1[, c("x", "y")]


    # Remove NA's from the first and last time interval
    # and due to insufficient obs within space.id
    grm.fit.cv.1 <- grm.fit.cv.1[!is.na(grm.fit.cv.1$estimate), ]
    grm.fit.cv.2 <- grm.fit.cv.2[!is.na(grm.fit.cv.2$estimate), ]



    # Use only first primary variable
    # results with second primary variable observed
    # and only second primary variable results
    # with first primary variable observed
    grm.fit.cv.1$link_id <- paste(grm.fit.cv.1$time.id, 
                                  grm.fit.cv.1$space.id, 
                                  sep = "_")
    grm.fit.cv.2$link_id <- paste(grm.fit.cv.2$time.id, 
                                  grm.fit.cv.2$space.id, 
                                  sep = "_")
    grm.fit.cv.1 <- grm.fit.cv.1[grm.fit.cv.1$link_id %in% grm.fit.cv.2$link_id, ]
    grm.fit.cv.2 <- grm.fit.cv.2[grm.fit.cv.2$link_id %in% grm.fit.cv.1$link_id, ] #new
    grm.fit.cv.2$estimate_1 <- grm.fit.cv.1$estimate[match(grm.fit.cv.2$link_id, 
                                                             grm.fit.cv.1$link_id)]

    grm.fit.cv.2$sd_1 <- grm.fit.cv.1$sd[match(grm.fit.cv.2$link_id,
                                                 grm.fit.cv.1$link_id)]


    # Calculate densities
    d1 <- stats::dnorm(grm.fit.cv.2$obs, 
                       grm.fit.cv.2$estimate_1, 
                       grm.fit.cv.2$sd_1)
    d2 <- stats::dnorm(grm.fit.cv.2$obs, 
                       grm.fit.cv.2$estimate, 
                       grm.fit.cv.2$sd)

    S <- max(grm.fit.cv.2$space.id)
    n <- length(d1)

    ###Create distance matrix
    locs <- unique(cbind(space.id.Y.1, coords.Y.1))
    locs <- locs[order(locs$space.id.Y.1), ]
    dist.space.mat <- as.matrix(stats::dist(locs[, c("x", "y")], 
                                            diag = TRUE, 
                                            upper = TRUE))
  
    ##Results to save
    #Total number of samples in the end
    K <- (n.iter - burn) / thin
    q.save <- matrix(NA, ncol = S, nrow = K)
    theta.save <- rep(NA, K)
    tau2.save <- rep(NA, K)
    dev.save <- rep(NA, K)
    theta.acc <- 0
  
    ##Initial values
    #w <- rep (0.5, S)
    q <- rep(0, S)
    z <- stats::rbinom(n, 1, 0.5)
    tau2 <- 0.01
    theta <- 100
    Sigma <- tau2 * exp(-1 / theta * dist.space.mat)

    ##Update W
    for (i in 1:n.iter){

        if ((i %% 1000) == 0) {
            cat(paste("  Iteration", i, "of", n.iter, "\n"))
        }
 
        ##Update log-odds q
        m <- q[grm.fit.cv.2$space.id]
        omega <- BayesLogit::rpg.devroye(n, 1, m)
    
        O <- diag(tapply(omega, grm.fit.cv.2$space.id, sum))
        VVV <- solve(O + solve(Sigma))
        MMM <- VVV %*% (tapply(z - 0.5, grm.fit.cv.2$space.id, sum))  
        q <- MASS::mvrnorm(1, MMM, VVV)
    
        #update tau
        SSS <- t(q) %*% solve(exp(-dist.space.mat / theta)) %*% q
        tau2 <- 1 / stats::rgamma(1, S / 2 + tau.a, SSS / 2 + tau.b)
        Sigma <- tau2 * exp(-1 / theta * dist.space.mat)
    
        #Update theta
        theta.prop <- stats::rlnorm(1, log(theta), theta.tune)
        SSS.curr <- Sigma
        SSS.prop <- tau2 * exp(-dist.space.mat / theta.prop)
    
        lik.curr <- mvtnorm::dmvnorm(q, rep(0, S), SSS.curr, log = T)
        lik.prop <- mvtnorm::dmvnorm(q, rep(0, S), SSS.prop, log = T)
    
        ratio <- lik.prop + 
            stats::dgamma(theta.prop, theta.a, theta.b, log = T) + 
            log(theta.prop) -
            lik.curr - 
            stats::dgamma(theta, theta.a, theta.b, log = T) - 
            log(theta)
    
        if (log(stats::runif(1)) < ratio) {
            theta <- theta.prop
            theta.acc <- theta.acc + 1
        }
    
        ##Update indicator z
        w <- 1 / (1 + exp(-q))[grm.fit.cv.2$space.id]
        p <- w * d1 / (w * d1 + (1 - w) * d2)
        z <- stats::rbinom(n, 1, p)
    
        if (i > burn & ((i - burn) %% thin)  == 0) {
            k <- (i - burn) / thin
            q.save[k, ] <- q
            tau2.save[k] <- tau2
            theta.save[k] <- theta
            dev.save[k] <- -2 * sum(log(w * d1 + (1 - w) * d2))
        }
    }
  
    q.save <- data.frame(space.id = 1:S, t(q.save))
    names(q.save) <- c("space.id", paste0("Sample",1:K))
  
    other.save <- data.frame(tau2 = tau2.save, 
                            theta = theta.save, 
                            dev = dev.save)
  
    list(q = q.save, 
         other = other.save, 
         theta.acc = theta.acc / n.iter,
         locations = locs)

}

