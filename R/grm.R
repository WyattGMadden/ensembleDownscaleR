
#' Fit the Geostatistical Regression Model (GRM)
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M
#'
#' @param Y Dependent variable vector (n)
#' @param X Unstandardized primary independent variable vector (n)
#' @param L Unstandardized spatial covariate matrix (n, p1)
#' @param M Unstandardized spatio-temporal covariate matrix (n, p2)
#' @param coords Matrix of x y coordinates, with colnames(coords) == c("x", "y"), (n, 2)
#' @param space.id Spatial location ID vector (n)
#' @param time.id Temporal location ID vector (n)
#' @param spacetime.id ID vector of time points where spatial trends vary (n)
#' @param include.additive.temporal.effect Include additive temporal random effects.
#' @param include.additive.spatial.effect Include additive spatial random effects.
#' @param include.multiplicative.temporal.effect Include multiplicative temporal random effects.
#' @param include.multiplicative.spatial.effect Include multiplicative spatial random effects.
#' @param num_neighbors Number of nearest neighbors to use in NNGP
#' @param n.iter Number of iterations used in predictions. 
#' @param burn Number of pre-covergence simulations
#' @param thin Save every thin'th simulation
#' @param covariance Specify covariance function (from "exponential", "matern", "custom")
#' @param covariance.kernal Specify a custom covariance function if covariance = "custom". Must be a function with "distance" and "theta" parameters.
#' @param matern.nu Specify nu parameter for Matern covariance function if used (from 0.5, 1.5, and 2.5)
#' @param tau.alpha.a First tau alpha prior hyperparameter
#' @param tau.alpha.b Second tau alpha prior hyperparameter
#' @param tau.beta.a First tau beta prior hyperparameter
#' @param tau.beta.b Second tau beta prior hyperparameter
#' @param omega.a First omega prior hyperparameter
#' @param omega.b Second omega prior hyperparameter
#' @param theta.alpha.tune Theta alpha Metropolis-Hastings proposal tuning parameter
#' @param theta.alpha.a First theta alpha prior hyperparameter
#' @param theta.alpha.b Second theta alpha prior hyperparameter
#' @param theta.alpha.init Initial value for theta alpha
#' @param theta.beta.tune Theta beta Metropolis-Hastings proposal tuning parameter
#' @param theta.beta.a First theta beta prior hyperparameter
#' @param theta.beta.b Second theta beta prior hyperparameter
#' @param theta.beta.init Initial value for theta beta
#' @param rho.alpha.init Initial value for rho alpha
#' @param rho.beta.init Initial value for rho beta
#' @param sigma.a First sigma prior hyperparameter
#' @param sigma.b Second sigma prior hyperparameter
#' @param verbose Print MCMC output
#' @param verbose.iter print MCMC output step number each 'verbose.iter' iterations
#'
#' @return A list containing MCMC output 
#'
#' @examples
#' # grm()
#' 
#' 
#' @export
grm <- function(
    Y, 
    X, 
    L = NULL, 
    M = NULL, 
    coords, 
    space.id, 
    time.id, 
    spacetime.id,
    include.additive.temporal.effect = T,
    include.additive.spatial.effect = T,
    include.multiplicative.temporal.effect = T,
    include.multiplicative.spatial.effect = T,
    num_neighbors = 10,
    n.iter = 25000,
    burn = 5000,
    thin = 4,
    covariance = "exponential",
    covariance.kernal = NULL,
    matern.nu = 1.5,
    tau.alpha.a = 0.5,
    tau.alpha.b = 0.005,
    tau.beta.a = 0.5,
    tau.beta.b = 0.005,
    omega.a = 0.5,
    omega.b = 0.005,
    theta.alpha.tune = 0.2, 
    theta.alpha.a = 5, 
    theta.alpha.b = 0.05,
    theta.alpha.init = 100,
    theta.beta.tune = 0.2, 
    theta.beta.a = 5, 
    theta.beta.b = 0.05,
    theta.beta.init = 100,
    rho.alpha.init = 0.9999,
    rho.beta.init = 0.9999,
    sigma.a = 0.001, 
    sigma.b = 0.001,
    verbose = TRUE,
    verbose.iter = 1000
    ) {


    ########################
    ### Argument Checks ####
    ########################

    # Check basic types and lengths
    if (!is.numeric(Y) || !is.vector(Y)) {
        stop("'Y' must be a numeric vector.")
    }
    if (!is.numeric(X) || !is.vector(X)) {
        stop("'X' must be a numeric vector.")
    }
    if (length(Y) != length(X)) {
        stop("'Y' and 'X' must have the same length.")
    }
    if (!is.null(L)) {
        if (!is.matrix(L) | !is.data.frame(L)) {
            stop("'L' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(L) != length(Y)) {
            stop("Number of rows in 'L' must match length of 'Y'.")
        }
    }
    if (!is.null(M)) {
        if (!is.matrix(M) | !is.data.frame(M)) {
            stop("'M' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(M) != length(Y)) {
            stop("Number of rows in 'M' must match length of 'Y'.")
        }
    }
    # coords checks
    if (!is.matrix(coords)) {
        stop("'coords' must be a matrix.")
    }
    if (nrow(coords) != length(Y)) {
        stop("Number of rows in 'coords' must match length of 'Y'.")
    }
    if (!all(colnames(coords) == c("x", "y"))) {
        stop("'coords' must have colnames 'x' and 'y'.")
    }
    # space.id / time.id / spacetime.id checks
    if (length(space.id) != length(Y)) {
        stop("'space.id' must have the same length as 'Y'.")
    }
    if (length(time.id) != length(Y)) {
        stop("'time.id' must have the same length as 'Y'.")
    }
    if (length(spacetime.id) != length(Y)) {
        stop("'spacetime.id' must have the same length as 'Y'. If you do not wish to include spacetime effects, input a vector of 1's with the same length as space.id and time.id.")
    }
    # Boolean checks
    if (!is.logical(include.additive.temporal.effect) ||
        !is.logical(include.additive.spatial.effect) ||
        !is.logical(include.multiplicative.temporal.effect) ||
        !is.logical(include.multiplicative.spatial.effect)) {
        stop("'include.*.effect' arguments must be TRUE/FALSE.")
    }
    # numeric / integer checks
    if (!is.numeric(num_neighbors) || length(num_neighbors) != 1 || num_neighbors < 1) {
        stop("'num_neighbors' must be a positive numeric scalar.")
    }
    if (!is.numeric(n.iter) || n.iter <= 0) {
        stop("'n.iter' must be a positive numeric scalar.")
    }
    if (!is.numeric(burn) || burn < 0) {
        stop("'burn' must be a nonnegative numeric scalar.")
    }
    if (!is.numeric(thin) || thin < 1) {
        stop("'thin' must be a positive numeric scalar.")
    }
    if (!is.character(covariance) || !(covariance %in% c("exponential", "matern", "custom"))) {
        stop("'covariance' must be one of 'exponential', 'matern', or 'custom'.")
    }
    if (covariance == "custom" && !is.function(covariance.kernal)) {
        stop("If 'covariance' is 'custom', 'covariance.kernal' must be a valid function.")
    }
    if (covariance == "matern" && !(matern.nu %in% c(0.5, 1.5, 2.5))) {
        stop("'matern.nu' must be 0.5, 1.5, or 2.5 when 'covariance' is 'matern'.")
    }
    # prior hyperparameters (non-exhaustive)
    if (!is.numeric(sigma.a) || !is.numeric(sigma.b) || sigma.a <= 0 || sigma.b <= 0) {
        stop("'sigma.a' and 'sigma.b' must be positive numeric values.")
    }
    if (!is.numeric(theta.alpha.a) || theta.alpha.a <= 0) {
        stop("'theta.alpha.a' must be positive.")
    }
    if (!is.numeric(theta.alpha.b) || theta.alpha.b <= 0) {
        stop("'theta.alpha.b' must be positive.")
    }
    if (!is.numeric(theta.beta.a) || theta.beta.a <= 0) {
        stop("'theta.beta.a' must be positive.")
    }
    if (!is.numeric(theta.beta.b) || theta.beta.b <= 0) {
        stop("'theta.beta.b' must be positive.")
    }
    if (!is.numeric(verbose.iter) || verbose.iter < 1) {
        stop("'verbose.iter' must be a positive integer.")
    }
    # checks for initial values
    if (!is.numeric(theta.alpha.init) || !is.numeric(theta.beta.init)) {
        stop("Initial values for theta.alpha.init and theta.beta.init must be numeric.")
    }
    if (!is.numeric(rho.alpha.init) || !is.numeric(rho.beta.init)) {
        stop("Initial values for rho.alpha.init and rho.beta.init must be numeric.")
    }

 

    ###############################
    ### Get summary statistics ####
    ###############################
  
    N = length(Y) #Total AOD-PM linked pairs
    N.space <- max(space.id) #Total number of monitors in training
    N.time <- max(time.id) #Maximum number of time interval in training
    N.time.obs <- length(unique(time.id)) #Number of observed time interval
    N.spacetime <- max(spacetime.id) #Time points where spatial trends vary by (year)
    
    N.Lmax <- ncol(L) #Number of spatial predictors to use
    N.Mmax <- ncol(M) #Number of spatial-temporal predictors to use
    
    if (verbose == TRUE) {
      print("#################################### ")
      print("######## Preparing for MCMC ######## ")
      print("#################################### ")
      print(paste0("     Total number of available days: ", 
                    N.time.obs, 
                    " (out of ", N.time,")")
      )
      print(paste0("     Total number of monitors: ", 
                   N.space)
      )
      print(paste0("     Total number of spatial covariates: ", 
                   ncol(L))
      )
      print(paste0("     Total number of spatial-temporal covariates: ", 
                   ncol(M))
      )
    }

    #################################
    ### Specify Covariance Kernal ###
    #################################

    if (covariance == "exponential") {

      cov_kern <- exponential_kernal

    } else if (covariance == "matern") {

      cov_kern <- function(distance, theta) {
          matern_kernal(distance = distance, 
                        theta = theta, 
                        nu = matern.nu)
      }
    } else if (covariance == "custom") {

      cov_kern <- covariance.kernal
      
    } else {

      stop("Covariance function not recognized")

    }

    ##############################
    ### Standardize X, L and M ###
    ##############################
    
    X.mean <- mean(X)
    X.sd <- stats::sd(X)
    X <- (X - X.mean) / X.sd

    if (!is.null(L)) {
        L.temp <- covariate_matrix_standardize(L)
        L <- L.temp$x
        L.mean <- L.temp$x.mean
        L.sd <- L.temp$x.sd
    }

    if (!is.null(M)) {
        M.temp <- covariate_matrix_standardize(M)
        M <- M.temp$x
        M.mean <- M.temp$x.mean
        M.sd <- M.temp$x.sd
    }


    ######################################
    ### Create Spatial Distance Matrix ###
    ######################################

    locs <- unique(cbind(space.id, coords))
    locs <- locs[order(locs$space.id), ]
    dist.space.mat <- as.matrix(stats::dist(locs[, c("x", "y")], 
                                            diag = TRUE, 
                                            upper = TRUE))

    
    #######################################
    ### Create Temporal Distance Matrix ###
    #######################################
    
    ### Create full temporal distance matrix
    A <- spam::as.spam(matrix(0, ncol=N.time, N.time)) #Adjacency matrix
    for (i in 1:(N.time - 1)) {
        A[i, i + 1] <- 1
    } 

    for (i in 2:(N.time)) { 
        A[i, i - 1] <- 1 
    }

    Time.obs <- 1:N.time %in% time.id ## T/F for which days are observed
    
    Q <- spam::diag.spam(apply(A, 2, sum)) #Row sums of A 


    
    ####################################################
    ### Pre-calculate quantities to save computation ###
    ####################################################

    ### space assignment to spacetime vector
    space_to_spacetime_assign <- rep(unique(spacetime.id), each = N.space)
    
    ###Indicator Matrix for Time
    Gamma_time <- matrix(0, nrow = N, ncol = N.time)
    for (i in 1:N){
        Gamma_time[i, time.id[i]] <- 1
    }
    Gamma_time <- spam::as.spam(Gamma_time)
    GtG_time <- spam::diag.spam(t(Gamma_time) %*% Gamma_time)
    
    ###Indicator Matrix for SpaceTime
    Gamma_space <- matrix(0, nrow = N, ncol = N.spacetime * N.space)
    for (i in 1:N) {
        Gamma_space[i, (spacetime.id[i]-1) * N.space + space.id[i]] <- 1
    }
    Gamma_space <- spam::as.spam(Gamma_space)
    GtG_space <- spam::diag.spam(t(Gamma_space) %*% Gamma_space)
    Z_ID <- (spacetime.id-1) * N.space + space.id
    
    
    ### Design matrix for temporal random effects
    X_W <- NULL
    #GammatX_W <- NULL
    for (time.i in (1:N.time)) {
        X.i <- X[time.id == time.i]
      
        if (length (X.i)>0){
            X_W <- c(X_W, t(X.i) %*% (X.i))
        } else {
            X_W <- c(X_W, 0)
        }
    }
    
    ### Design matrix for spatial random effects
    X_S <- NULL
    GammatX_S <- NULL
    Space.perm <- NULL
    SpaceTime.avail <- NULL
    
    for (spacetime.i in 1:N.spacetime) {
        for (mon.i in 1:N.space) {
            use <- which(space.id == mon.i & spacetime.id == spacetime.i) 
            SpaceTime.avail <- rbind(SpaceTime.avail, 
                                    c(spacetime.i, mon.i, length(use)))
          
            if (length(use) != 0) {
                Space.perm <- c(Space.perm, use)
                X.i <- as.matrix(X[use])
                X_S <- c(X_S, t(X.i) %*% X.i)
            } else {
                X_S <- c(X_S, 0)
            }
        }
    }
    
    ### Spatial and Spatial-temporal coefficients
    if (!is.null(M)) MtM <- t(M) %*% M
    if (!is.null(L)) LtL <- t(L) %*% L
    XtX <- t(cbind(1, X)) %*% cbind(1, X)
    XtX_inv <- solve(XtX)


    ### Discretize CAR parameters
    nrho <- 2000
    jun <- 1 / sqrt(Q)
    dt <- eigen(jun %*% A %*% jun)$values
    canrho <- detpart <- stats::qbeta(seq(0, 0.9999, length = nrho), 1, 1)
    for (j in 1:nrho) {
        detpart[j] <- 1 / 2 * sum(log(1 - canrho[j] * dt)) 
    }
    
    
    #####################
    ### Initial Value ###
    #####################
    
    ###Initialize alpha0, beta0, gamma, and delta
    ###Initialize lambda if needed
    if (!is.null(L) & !is.null(M)) {
        fit <- stats::lm(Y ~ cbind(X, L, M))
        alpha0 <- stats::coef(fit)[1]
        beta0 <- stats::coef(fit)[2]
        gamma <- stats::coef(fit)[3:(ncol(L) + 2)]
        delta <- stats::coef(fit)[(3 + ncol(L)):length(stats::coef(fit))]
        lambda_gamma <- stats::var(gamma)
        lambda_delta <- stats::var(delta)
        mu <- alpha0 + beta0 * X + L %*% gamma + M %*% delta
    } 
    
    if (is.null(L) & !is.null(M)) {
        fit <- stats::lm (Y ~ cbind(X, M))
        alpha0 <- stats::coef(fit)[1]
        beta0 <- stats::coef(fit)[2]
        delta <- stats::coef(fit)[3:length(stats::coef(fit))]
        lambda_delta <- stats::var(delta)
        mu <- alpha0 + beta0 * X + M %*% delta
    }
    
    if (!is.null(L) & is.null(M)) {
        fit <- stats::lm (Y ~ cbind(X, L))
        alpha0 <- stats::coef(fit)[1]
        beta0 <- stats::coef(fit)[2]
        gamma <- stats::coef(fit)[3:(ncol(L)+2)]
        lambda_gamma <- stats::var(gamma)
        mu <- alpha0 + beta0 * X + L %*% gamma 
    }
    
    if (is.null(L) & is.null(M)) {
        fit <- stats::lm (Y ~ X)
        alpha0 <- stats::coef(fit)[1]
        beta0 <- stats::coef(fit)[2]
        mu <- alpha0 + beta0 * X 
    } 
    
    
    ###Initialize alpha_time and its CAR variance omega_alpha
    alpha_time <- rep(0, N.time)
    omega_alpha <- 0
    if (include.additive.temporal.effect) { 
        alpha_time <- (1 / GtG_time) * t(Gamma_time) %*% (Y - mu)
        omega_alpha <- as.numeric(stats::var(alpha_time, na.rm = T))
    }
    
    ###Initialize beta_time and its CAR variance omega_beta
    beta_time <- rep(0, N.time)
    omega_beta <- 0
    if (include.multiplicative.temporal.effect) { 
        RRR <- Y - mu - alpha_time[time.id]
        beta_time <-  1 / X_W * t(Gamma_time) %*% (X * RRR)
        omega_beta <- as.numeric(stats::var(beta_time, na.rm = T))
    }
    
    ###Initialize alpha_spacetime and its spatial variance tau_alpha
    alpha_space <- rep(0, N.spacetime * N.space)
    tau_alpha <- 0
    if (include.additive.spatial.effect) {
        RRR <- Y - mu - alpha_time[time.id] - beta_time[time.id] * X
        alpha_space <- as.numeric((1 / GtG_space) * t(Gamma_space) %*% RRR)
        tau_alpha <- as.numeric(stats::var(alpha_space, na.rm = T))
    }
    
    ###Initialize beta_spacetime and its spatial variance tau_beta
    beta_space <- rep(0, N.spacetime * N.space)
    tau_beta <- 0
    if (include.multiplicative.spatial.effect) {
        RRR <- Y - mu - alpha_time[time.id] - beta_time[time.id] * X - alpha_space[Z_ID]
        beta_space <- as.numeric(1 / X_S * t(Gamma_space) %*% (X * RRR))
        tau_beta <- as.numeric(stats::var(beta_space, na.rm = T))
    }
    
    ###Initialize mean
    ###Initialize mean
    MMM <- alpha0 + beta0 * X +
        alpha_time[time.id] +
        beta_time[time.id] * X +
        alpha_space[Z_ID] +
        beta_space[Z_ID] * X

    if (!is.null(L)) {
        MMM <- MMM + L %*% gamma
    }
    if (!is.null(M)) {
        MMM <- MMM + M %*% delta
    }
    
    ###Initializae sigma2
    sigma2 <- stats::var(as.numeric(Y - MMM), na.rm = T)
    theta_alpha <- theta.alpha.init
    theta_beta <- theta.beta.init
    rho_alpha <- rho.alpha.init
    rho_beta <- rho.beta.init
    
    #############################
    ### MCMC saved parameters ###
    #############################
    ###Declare  matrices and vectors to save results
    
    ###Total number of samples in the end
    K <- (n.iter - burn) / thin
    
    alpha0.save <- rep(NA, K)
    beta0.save <- rep(NA, K)
    
    gamma.save <- NULL
    if (!is.null(L)) {
        gamma.save <- matrix(NA, nrow = K, ncol = length(gamma))
        lambda_gamma.save <- rep(NA, K)
    }

    delta.save <- NULL
    if (!is.null(M)) {
        delta.save <- matrix(NA, nrow = K, ncol = length(delta))
        lambda_delta.save <- rep(NA, K)
    }
    
    alpha_time.save <- matrix(NA, nrow = K, ncol = N.time)
    beta_time.save <- matrix(NA, nrow = K, ncol = N.time)
    
    alpha_space.save <- matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    beta_space.save <- matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    
    sigma2.save <- rep(NA, K)
    theta_alpha.save <- theta_beta.save <- rep(NA, K)
    rho_alpha.save <- rho_beta.save <- rep(NA, K)
    
    omega_alpha.save <- omega_beta.save <- rep(NA, K)
    lambda_alpha.save <- lambda_beta.save <- rep(NA, K)
    tau_alpha.save <- tau_beta.save <- rep(NA, K)
    
    LL.save <- rep(NA, K)
    
    Y.hat <- rep(0, length(Y))
    
    tau.acc <- c(0,0)
    theta.acc <- c(0,0)


    
    ###########################
    ### BEGIN MCMC SAMPLING ###
    ###########################
    if (verbose == TRUE) {
      print ("#################################### ")
      print ("########    MCMC Begins!    ######## ")
      print ("#################################### ")
    }

    for (i in 1:n.iter) {

        if (verbose == TRUE) {
            if ((i %% verbose.iter) == 0) print(paste("     Iteration", i, "of", n.iter))
        }
    
        ###Update alpha0 and beta0
        MMM <-  MMM - alpha0 - beta0 * X 
        RRR <- Y - MMM
        XXX <- t(cbind(1,X)) %*% RRR
        VVV <- XtX_inv
        alpha0_beta0 <- matrix(mvnfast::rmvn(1, VVV %*% XXX, sigma2 * VVV), ncol = 1)
        alpha0 <- alpha0_beta0[1]
        beta0 <- alpha0_beta0[2]
        MMM <- MMM + alpha0 + beta0 * X
      
        #Update gamma
        if (!is.null(L)) {
            MMM <- MMM - L %*% gamma
            RRR <-  Y - MMM
            XXX <- 1 / sigma2 * t(L) %*% RRR
            VVV <- 1 / sigma2 * LtL + 1 / lambda_gamma * diag(N.Lmax)
            VVV <- solve(VVV)
            gamma <- matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
            MMM <- MMM + L %*% gamma

            ##Update lambda_gamma
            lambda_gamma <- 1 / stats::rgamma(1,  
                                             length(gamma) / 2 + sigma.a, 
                                             sum(gamma^2)/2 + sigma.b)
        }
      
        #Update delta
        if (!is.null(M)) {
            MMM <- MMM - M %*% delta
            RRR <-  Y - MMM
            XXX <- 1 / sigma2 * t(M) %*% RRR
            VVV <- 1 / sigma2 * MtM + 1 / lambda_delta * diag(N.Mmax)
            VVV <- solve(VVV)
            delta <- matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
            MMM <- MMM + M %*% delta
        
            ##Update lambda_delta
            lambda_delta <- 1 / stats::rgamma(1, 
                                             length(delta) / 2 + sigma.a, 
                                             sum(delta ^ 2) / 2 + sigma.b)
        }
      
        #Update residual error sigma2
        RRR <- Y - MMM
        sigma2 <- 1 / stats::rgamma(1, length(RRR) / 2 + sigma.a, sum(RRR ^ 2) / 2 + sigma.b)
       

        if (include.additive.spatial.effect) {

            MMM <- MMM - alpha_space[Z_ID]
            RRR <- Y - MMM
            XXX <- 1 / sigma2 * t(Gamma_space) %*% RRR
            kern <- cov_kern(distance = dist.space.mat, theta = theta_alpha)
            kern_inv <- solve(kern)
            SSS <- tau_alpha * kern
            SSS_inv <- (1 / tau_alpha) * kern_inv
            for (st in unique(spacetime.id)) {
                GtG_space_st <- GtG_space[space_to_spacetime_assign == st]
                XXX_st <- XXX[space_to_spacetime_assign == st]
                VVV_st <- diag(1 / sigma2 * GtG_space_st) + SSS_inv
                VVV_inv_st <- solve(VVV_st)
                alpha_space_st <- as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 
                                                         VVV_inv_st))
                alpha_space[space_to_spacetime_assign == st] <- alpha_space_st
            }
            MMM <- MMM + alpha_space[Z_ID]
      
            #update tau_alpha
            SSS <- 0
            for (st in unique(spacetime.id)) {
                alpha_space_st <- alpha_space[space_to_spacetime_assign == st]
                SSS_st <- t(alpha_space_st) %*% kern_inv %*% alpha_space_st
                SSS <- SSS + SSS_st
            }
            tau_alpha <- 1 / stats::rgamma(1, 
                                          N.space * N.spacetime / 2 + tau.alpha.a, 
                                          SSS / 2 + tau.alpha.b)
      
            #Update theta_alpha
            theta.prop <- stats::rlnorm(1, 
                                       log(theta_alpha), 
                                       theta.alpha.tune)
            SSS.curr <- tau_alpha * kern
            SSS.prop <- tau_alpha * cov_kern(distance = dist.space.mat, 
                                            theta = theta.prop)
      
            lik.curr <- 0
            lik.prop <- 0

            for (st in unique(spacetime.id)) {
                alpha_space_st <- alpha_space[space_to_spacetime_assign == st]
                lik.prop_st <- mvtnorm::dmvnorm(alpha_space_st, 
                                               rep(0, N.space), 
                                               SSS.prop, 
                                               log = T)
                lik.curr_st <- mvtnorm::dmvnorm(alpha_space_st,
                                               rep(0, N.space), 
                                               SSS.curr, 
                                               log = T)
                lik.prop <- lik.prop + lik.prop_st
                lik.curr <- lik.curr + lik.curr_st
            }
      
            ratio <- lik.prop + 
                stats::dgamma(theta.prop, 
                              theta.alpha.a, 
                              theta.alpha.b, 
                              log = T) + 
                    log(theta.prop) -
                    lik.curr - 
                    stats::dgamma(theta_alpha, 
                                  theta.alpha.a, 
                                  theta.alpha.b, 
                                  log = T) - 
                    log(theta_alpha)
            if (log(stats::runif(1)) < ratio) {
                theta_alpha <- theta.prop
                theta.acc[1] <- theta.acc[1] + 1
            }

        }

      
        #Update spatial coefficent for AOD if GP
        if (include.multiplicative.spatial.effect) {

                    MMM <- MMM - beta_space[Z_ID] * X
                    RRR <- Y - MMM
                    XXX <- 1 / sigma2 * t(Gamma_space) %*% (X * RRR)
                    kern <- cov_kern(distance = dist.space.mat, 
                                    theta = theta_beta)
                    kern_inv <- solve(kern)
                    SSS <- tau_beta * kern
                    SSS_inv <- (1 / tau_beta) * kern_inv
                    for (st in unique(spacetime.id)) {
                        X_S_st <- X_S[space_to_spacetime_assign == st]
                        XXX_st <- XXX[space_to_spacetime_assign == st]
                        VVV_st <- diag(1 / sigma2 * X_S_st)
                        VVV_st <- VVV_st + SSS_inv
                        VVV_inv_st <- solve(VVV_st)
                        beta_space_st <- as.vector(mvnfast::rmvn(1, VVV_inv_st %*% XXX_st, 
                                                                VVV_inv_st))
                        beta_space[space_to_spacetime_assign == st] <- beta_space_st
                    }
                    MMM <- MMM + beta_space[Z_ID] * X

              
                    #update tau_beta
                    SSS <- 0
                    for (st in unique(spacetime.id)) {
                        beta_space_st <- beta_space[space_to_spacetime_assign == st]
                        SSS_st <- t(beta_space_st) %*% 
                            kern_inv %*% 
                            beta_space_st
                        SSS <- SSS + SSS_st
                    }
                    tau_beta <- 1 / stats::rgamma(1, N.space * N.spacetime /2 + tau.beta.a, SSS / 2 + tau.beta.b)
              
                    #Update theta_beta
                    theta.prop <- stats::rlnorm(1, log(theta_beta), theta.beta.tune)
                    SSS.curr <- tau_beta * kern
                    SSS.prop <- tau_beta * cov_kern(distance = dist.space.mat, 
                                                   theta = theta.prop)

                    lik.curr <- 0
                    lik.prop <- 0
                    for (st in unique(spacetime.id)) {
                        beta_space_st <- beta_space[space_to_spacetime_assign == st]
                        lik.prop_st <- mvtnorm::dmvnorm(beta_space_st, 
                                                       rep(0, N.space), 
                                                       SSS.prop, 
                                                       log = T)
                        lik.curr_st <- mvtnorm::dmvnorm(beta_space_st,
                                                       rep(0, N.space), 
                                                       SSS.curr, 
                                                       log = T)
                        lik.prop <- lik.prop + lik.prop_st
                        lik.curr <- lik.curr + lik.curr_st
                    }
              
                    ratio <- lik.prop + 
                        stats::dgamma(theta.prop,  
                                      theta.beta.a, 
                                      theta.beta.b, 
                                      log = T) + 
                        log(theta.prop) -
                        lik.curr - 
                        stats::dgamma(theta_beta, 
                                      theta.beta.a, 
                                      theta.beta.b, 
                                      log = T) - 
                        log(theta_beta)

                    if (log(stats::runif(1)) < ratio) {
                        theta_beta <- theta.prop
                        theta.acc[2] <- theta.acc[2] + 1
                    }

        }
      
        #Update temporal intercepts and its parameters
        if (include.additive.temporal.effect) {
            MMM = MMM - alpha_time[time.id]
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_time) %*% RRR
            VVV = diag(1 / sigma2 * GtG_time) + 1 / omega_alpha * (Q - rho_alpha * A)
            VVV = solve(VVV)
            alpha_time = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1,]
            MMM = MMM + alpha_time[time.id]
      
            #Update omega_alpha
            SS1 = alpha_time %*% Q %*% alpha_time
            SS2 = alpha_time %*% A %*% alpha_time
            omega_alpha = 1 / stats::rgamma(1, N.time / 2 + omega.a, (SS1 - rho_alpha * SS2) / 2 + omega.b)
            
            #Update rho_alpha
            R = detpart + 0.5 * omega_alpha * canrho * c(SS2)
            rho_alpha = sample(canrho, 1, prob = exp(R - max(R)))
        }

        #Update temporal AOD coefficients
        if (include.multiplicative.temporal.effect) {
            MMM = MMM - beta_time[time.id] * X
            RRR = Y - MMM
            XXX = 1 / sigma2 * t(Gamma_time) %*% (X * RRR)
            VVV = diag(1 / sigma2 * X_W) + 1 / omega_beta * (Q - rho_beta * A)
            VVV = solve(VVV)
            beta_time = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1,]
            MMM = MMM + beta_time[time.id] * X
        
            #Update omega_beta
            SS1 = beta_time %*% Q %*% beta_time
            SS2 = beta_time %*% A %*% beta_time
            omega_beta = 1 / stats::rgamma(1, N.time / 2 + omega.a, (SS1 - rho_beta * SS2) / 2 + omega.b)
            
            #Update rho_beta
            R = detpart + 0.5 * omega_beta * canrho * c(SS2)
            rho_beta = sample(canrho, 1, prob = exp(R - max(R)))
        }
    
     
        ###Save Samples 
     
        if (i > burn & ((i - burn) %% thin == 0)) {
            k = (i - burn) / thin

            #Save statistics
            LL.save[k] = sum(-2 * stats::dnorm(Y, MMM, sqrt(sigma2), log = T))
  
            alpha0.save[k] = alpha0
            beta0.save[k] = beta0

            if (!is.null(L)) {
                gamma.save[k, ] <- gamma
                lambda_gamma.save[k] <- lambda_gamma
            }

            if (!is.null(M)) {
                delta.save[k, ] <- delta
                lambda_delta.save[k] <- lambda_delta
            }
       
       
            alpha_time.save[k,] = alpha_time
            beta_time.save[k,] = beta_time
       
            alpha_space.save[k,] = alpha_space
            beta_space.save[k,] = beta_space
       
            sigma2.save[k] = sigma2
            theta_alpha.save[k] = theta_alpha
            theta_beta.save[k] = theta_beta
            rho_alpha.save[k] = rho_alpha
            rho_beta.save[k] = rho_beta
       
            omega_alpha.save[k] = omega_alpha
            omega_beta.save[k] = omega_beta
            tau_alpha.save[k] = tau_alpha
            tau_beta.save[k] = tau_beta
       
            Y.hat = Y.hat + MMM / K

        }

    } ## End of MCMC iterations
    
    if (verbose == TRUE) {
        print("#################################### ")
    }

    ##Process for output
    if (!is.null(L)) {
        gamma.save <- data.frame(gamma.save)
        names(gamma.save) <- colnames(L)
    }

    if (!is.null(M)) {
        delta.save <- data.frame(delta.save)
        names(delta.save) <- colnames(M)
    }
    
    alpha_time.save = data.frame(time.id = 1:N.time, t(alpha_time.save))
    names(alpha_time.save) = c("time.id", 
                               paste0("Sample", 1:K)) 
    row.names(alpha_time.save) = NULL
    beta_time.save = data.frame(time.id = 1:N.time, 
                                t(beta_time.save))
    names(beta_time.save) = c("time.id", 
                              paste0("Sample", 1:K)) 
    row.names(alpha_time.save) = NULL
    
    alpha_space.save = data.frame(space.id = rep(1:N.space, N.spacetime),
                                  spacetime.id = rep(1:N.spacetime, each = N.space),  
                                  t(alpha_space.save))
    names(alpha_space.save) = c("space.id", 
                                "spacetime.id", 
                                paste0("Sample", 1:K)) 
    row.names(alpha_space.save) = NULL
    beta_space.save = data.frame(space.id = rep(1:N.space, N.spacetime),
                                 spacetime.id = rep(1:N.spacetime, each = N.space),  
                                t(beta_space.save))
    names(beta_space.save) = c("space.id", "spacetime.id", paste0("Sample", 1:K))
    row.names(beta_space.save) = NULL
    
    other.save = data.frame(alpha0 = alpha0.save, 
                            beta0 = beta0.save,
                            sigma2 = sigma2.save,
                            theta.alpha = theta_alpha.save, 
                            theta.beta = theta_beta.save, 
                            tau.alpha = tau_alpha.save, 
                            tau.beta = tau_beta.save, 
                            rho.alpha = rho_alpha.save, 
                            rho.beta = rho_beta.save, 
                            omega.alpha = omega_alpha.save, 
                            omega.beta = omega_beta.save,
                            dev = LL.save)

    if (!is.null(L)) {
        other.save$lambda.gamma <- lambda_gamma.save
    }
    if (!is.null(M)) {
        other.save$lambda.delta <- lambda_delta.save
    }   
    
    standardize.param = rbind(data.frame(Type = "X", 
                                         Name = "AOD/CTM", 
                                         Mean = X.mean, 
                                         SD = X.sd))
    if (!is.null(L)) {
        standardize.param <- rbind(standardize.param, 
                                  data.frame(Type = "L", 
                                             Name = colnames(L), 
                                             Mean = L.mean, 
                                             SD = L.sd))
    }
    if (!is.null(M)) {
        standardize.param <- rbind(standardize.param, 
                                  data.frame(Type = "M", 
                                             Name = colnames(M), 
                                             Mean = M.mean, 
                                             SD = M.sd))
    }

    row.names(standardize.param) = NULL

    
    list(delta = delta.save, 
         gamma = gamma.save,
         alpha.time = alpha_time.save, 
         beta.time = beta_time.save,
         alpha.space = alpha_space.save, 
         beta.space = beta_space.save,
         others = other.save, 
         Y = Y.hat, 
         standardize.param = standardize.param,
         theta.acc = theta.acc,
         tau.acc = tau.acc,
         cov_kern = cov_kern,
         locations = locs)
}


