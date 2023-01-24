
#' Fit the Geostatistical Regression Model (GRM)
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M
#'
#' @param Y Dependent variable vector (n)
#' @param X Unstandardized primary independent variable vector (n)
#' @param L Unstandardized spatial covariate matrix (n, p1)
#' @param M Unstandardized spatio-temporal covariate matrix (n, p2)
#' @param dist.space.mat Spatial distance matrix 
#' @param space.id Spatial location ID vector (n)
#' @param time.id Temporal location ID vector (n)
#' @param spacetime.id ID vector of time points where spatial trends vary (n)
#' @param include.additive.weekly.resid Include additive temporal (weekly) residual bias not explained by other inputes
#' @param include.additive.annual.resid Include additive spatial (annual) residual bias not explained by other inputes
#' @param include.multiplicative.weekly.resid Include multiplicative temporal (weekly) residual bias not explained by other inputes
#' @param include.multiplicative.annual.resid Include multiplicative spatial (weekly) residual bias not explained by other inputes
#' @param n.iter Number of MCMC simulations
#' @param burn Number of pre-covergence simulations
#' @param thin Save every thin'th simulation
#' @param tau.a First tau prior hyperparameter
#' @param tau.b Second tau prior hyperparameter
#' @param omega.a First omega prior hyperparameter
#' @param omega.b Second omega prior hyperparameter
#' @param theta.tune Theta Metropolis-Hastings proposal tuning parameter
#' @param theta.a First theta prior hyperparameter
#' @param theta.b Second theta prior hyperparameter
#' @param sigma.a First sigma prior hyperparameter
#' @param sigma.b Second sigma prior hyperparameter
#' @param verbose Print MCMC output
#'
#' @return A list containing MCMC output 
#'
#' @examples
#' # grm()
#' 
#' 
#' @export
grm = function(Y, 
               X, 
               L = NULL, 
               M = NULL, 
               dist.space.mat, 
               space.id, 
               time.id, 
               spacetime.id,
               include.additive.weekly.resid = T,
               include.additive.annual.resid = T,
               include.multiplicative.weekly.resid = T,
               include.multiplicative.annual.resid = T,
               n.iter = 25000,
               burn = 5000,
               thin = 4,
               tau.a = 0.5,
               tau.b = 0.005,
               omega.a = 0.5,
               omega.b = 0.005,
               theta.tune = 0.2, 
               theta.a = 5, 
               theta.b = 0.05,
               sigma.a = 0.001, 
               sigma.b = 0.001,
               verbose = TRUE) {
 
    ####################################
    ### Get some summary statistics ####
    ####################################
  
    N = length(Y) #Total AOD-PM linked pairs
    N.space = max(space.id) #Total number of monitors in training
    N.time = max(time.id) #Maximum number of time interval (weeks) in training
    N.time.obs = length(unique(time.id)) #Number of observed time interval (weeks)
    N.spacetime = max(spacetime.id) #Time points where spatial trends vary by (year)
    
    Space.labels = 1:N.space 
    Day.labels = 1:N.time
    
    N.Lmax = ncol(L) #Number of spatial predictors to use
    N.Mmax = ncol(M) #Number of spatial-temporal predictors to use
    
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

    ##############################
    ### Standardize X, L and M ###
    ##############################
    
    ### Calculate means and standard deviation
    X.mean = mean(X)
    if (!is.null(M)) M.mean = apply(M, 2, mean)
    if (!is.null(L)) L.mean = apply(L, 2, mean)
    
    X.sd = stats::sd(X)
    if (!is.null(M)) M.sd = apply(M, 2, stats::sd)
    if (!is.null(L)) L.sd = apply(L, 2, stats::sd)
    
    ### Standardize
    X = (X-mean (X))/X.sd
    if (!is.null(M) & !any(M.sd==0)) { 
        M = sweep(sweep(M, 2, M.mean, "-"), 2, M.sd, "/") 
    }

    if (!is.null(L) & !any(L.sd==0)){ 
        L = sweep(sweep(L, 2, L.mean, "-"), 2, L.sd, "/") 
    }
    
    #######################################
    ### Create Temporal Distance Matrix ###
    #######################################
    
    ### Create full temporal distance matrix
    A = spam::as.spam(matrix(0, ncol=N.time, N.time)) #Adjacency matrix
    for (i in 1:(N.time-1)){
        A[i, i + 1] = 1
    } 

    for (i in 2:(N.time)){ 
        A[i, i - 1] = 1 
    }

    Time.obs =  1:N.time %in% time.id ## T/F for which days are observed
    
    Q = spam::diag.spam(apply(A, 2, sum)) #Row sums of A 
    
    ####################################################
    ### Pre-calculate quantities to save computation ###
    ####################################################
    
    ###Indicator Matrix for Time
    Gamma_time = matrix(0, nrow = N, ncol = N.time)
    for (i in 1:N){
        Gamma_time[i, time.id[i]] <- 1
    }
    Gamma_time = spam::as.spam(Gamma_time)
    GtG_time = diag(t(Gamma_time) %*% Gamma_time) 
    
    ###Indicator Matrix for SpaceTime
    Gamma_space = matrix(0, nrow = N, ncol = N.spacetime * N.space)
    for (i in 1:N) {
        Gamma_space[i, (spacetime.id[i]-1) * N.space + space.id[i]] <- 1
    }
    Gamma_space = spam::as.spam(Gamma_space)
    GtG_space = diag(t(Gamma_space) %*% Gamma_space)
    Z_ID = (spacetime.id-1) * N.space + space.id
    
    
    ### Design matrix for temporal random effects
    X_W = NULL
    #GammatX_W = NULL
    for (time.i in (1:N.time)){
      X.i = X[time.id == time.i]
      
      if (length (X.i)>0){
        X_W = c(X_W, t(X.i) %*% (X.i))
        } else {
        X_W = c(X_W, 0)
      }
    }
    
    ### Design matrix for spatial random effects
    X_S = NULL
    GammatX_S = NULL
    Space.perm = NULL
    SpaceTime.avail = NULL
    
    for (spacetime.i in 1:N.spacetime) {
        for (mon.i in 1:N.space) {
          use = which (space.id == mon.i & spacetime.id == spacetime.i) 
          SpaceTime.avail = rbind(SpaceTime.avail, 
                                  c(spacetime.i, mon.i, length(use)))
          
          if (length(use) != 0) {
            Space.perm = c(Space.perm, use)
            X.i = as.matrix(X[use])
            X_S = c(X_S, t(X.i) %*% X.i)
          } else {
            X_S = c(X_S, 0)
          }
        }
    }
    
    ### Spatial and Spatial-temporal coefficients
    if (!is.null(M)) MtM = t(M) %*% M
    if (!is.null(L)) LtL = t(L) %*% L
    XtX = t(cbind(1, X)) %*% cbind(1, X)

    ### Discretize CAR parameters
    nrho = 2000
    jun = 1 / sqrt(Q)
    dt = eigen(jun %*% A %*% jun)$values
    canrho = detpart = stats::qbeta(seq(0, 0.9999, length = nrho), 1, 1)
    for (j in 1:nrho){ 
        detpart[j] = 1 / 2 * sum(log(1 - canrho[j] * dt)) 
    }
    
    
    #####################
    ### Initial Value ###
    #####################
    
    ###Initialize alpha0, beta0, gamma, and delta
    ###Initialize lambda if needed
    if (!is.null(L) & !is.null(M)) {
      fit = stats::lm(Y ~ cbind(X, L, M))
      alpha0 = stats::coef(fit)[1]
      beta0 = stats::coef(fit)[2]
      gamma = stats::coef(fit)[3:(ncol(L) + 2)]
      delta = stats::coef(fit)[(3 + ncol(L)):length(stats::coef(fit))]
      lambda_gamma = stats::var(gamma)
      lambda_delta = stats::var(delta)
      mu = alpha0 + beta0 * X + L %*% gamma + M %*% delta
    } 
    
    if (is.null(L) & !is.null(M)) {
      fit = stats::lm (Y ~ cbind(X, M))
      alpha0 = stats::coef(fit)[1]
      beta0 = stats::coef(fit)[2]
      delta = stats::coef(fit)[3:length(stats::coef(fit))]
      lambda_delta = stats::var(delta)
      mu = alpha0 + beta0 * X + M %*% delta
    }
    
    if (!is.null(L) & is.null(M)) {
      fit = stats::lm (Y ~ cbind(X, L))
      alpha0 = stats::coef(fit)[1]
      beta0 = stats::coef(fit)[2]
      gamma = stats::coef(fit)[3:(ncol(L)+2)]
      lambda_gamma = stats::var(gamma)
      mu = alpha0 + beta0 * X + L %*% gamma 
    }
    
    if (is.null(L) & is.null(M)) {
      fit = stats::lm (Y ~ X)
      alpha0 = stats::coef(fit)[1]
      beta0 = stats::coef(fit)[2]
      mu = alpha0 + beta0 * X 
    } 
    
    
    ###Initialize alpha_time and its CAR variance omega_alpha
    alpha_time = rep(0, N.time)
    omega_alpha = 0
    if (include.additive.weekly.resid) { 
      alpha_time = (1 / GtG_time) * t(Gamma_time) %*% (Y - mu)
      omega_alpha = as.numeric(stats::var(alpha_time, na.rm = T))
    }
    
    ###Initialize beta_time and its CAR variance omega_beta
    beta_time = rep(0, N.time)
    omega_beta = 0
    if (include.multiplicative.weekly.resid) { 
      RRR = Y - mu - alpha_time[time.id]
      beta_time=  1 / X_W * t(Gamma_time) %*% (X * RRR)
      omega_beta = as.numeric(stats::var(beta_time, na.rm = T))
    }
    
    ###Initialize alpha_spacetime and its spatial variance tau_alpha
    alpha_space = rep(0, N.spacetime * N.space)
    tau_alpha = 0
    if (include.additive.annual.resid) {
      RRR = Y - mu - alpha_time[time.id] - beta_time[time.id] * X
      alpha_space = as.numeric((1 / GtG_space) * t(Gamma_space) %*% RRR)
      tau_alpha = as.numeric(stats::var(alpha_space, na.rm = T))
    }
    
    ###Initialize alpha_spacetime and its spatial variance tau_beta
    beta_space = rep(0, N.spacetime * N.space)
    tau_beta = 0
    if (include.multiplicative.annual.resid) {
      RRR = Y - mu - alpha_time[time.id] - beta_time[time.id] * X - alpha_space[Z_ID]
      beta_space = as.numeric(1 / X_S * t(Gamma_space) %*% (X * RRR))
      tau_beta = as.numeric(stats::var(beta_space, na.rm = T))
    }
    
    ###Initializae mean
    MMM = alpha0 + beta0 * X + L %*% gamma+ M %*% delta + alpha_time[time.id] + beta_time[time.id] * X + alpha_space[Z_ID] + beta_space[Z_ID] * X
    
    ###Initializae sigma2
    sigma2 = stats::var(as.numeric(Y - MMM), na.rm = T)
    
    theta_alpha = 100
    theta_beta = 100
    rho_alpha = 0.99999
    rho_beta = 0.99999
    
    #############################
    ### MCMC saved parameters ###
    #############################
    ###Declare lots of matrices and vectors to save results
    
    ###Total number of samples in the end
    K = (n.iter - burn) / thin
    
    alpha0.save = rep(NA, K)
    beta0.save = rep(NA, K)
    
    gamma.save = matrix(NA, nrow = K, ncol = length(gamma))
    delta.save = matrix(NA, nrow = K, ncol = length(delta))
    
    lambda_gamma.save = rep(NA, K)
    lambda_delta.save = rep(NA, K)
    
    alpha_time.save = matrix(NA, nrow = K, ncol = N.time)
    beta_time.save = matrix(NA, nrow = K, ncol = N.time)
    
    alpha_space.save = matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    beta_space.save = matrix(NA, nrow = K, ncol = N.spacetime * N.space)
    
    sigma2.save = rep(NA, K)
    theta_alpha.save = theta_beta.save = rep(NA, K)
    rho_alpha.save = rho_beta.save = rep(NA, K)
    
    omega_alpha.save = omega_beta.save = rep(NA, K)
    lambda_alpha.save = lambda_beta.save = rep(NA, K)
    tau_alpha.save = tau_beta.save = rep(NA, K)
    
    LL.save = rep(NA, K)
    
    Y.hat = rep(0, length(Y))
    
    theta.acc = c(0,0)
    
    
    ###########################
    ### BEGIN MCMC SAMPLING ###
    ###########################
    if (verbose == TRUE){
      print ("#################################### ")
      print ("########    MCMC Begins!    ######## ")
      print ("#################################### ")
    }

    for (i in 1:n.iter){

      if (verbose == TRUE) {
        if ((i %% 1000) == 0) print(paste("     Iteration", i, "of", n.iter))
      }
    
      ###Update alpha0 and beta0
      MMM =  MMM - alpha0 - beta0 * X 
      RRR = Y - MMM
      XXX = t(cbind(1,X)) %*% RRR
      VVV = solve(XtX)
      junk = matrix(mvnfast::rmvn(1, VVV %*% XXX, sigma2 * VVV), ncol = 1)
      alpha0 = junk[1]
      beta0 = junk[2]
      MMM = MMM + alpha0 + beta0 * X
      
      #Update gamma
      if (!is.null(L)) {
        MMM = MMM - L %*% gamma
        RRR =  Y - MMM
        XXX = 1 / sigma2 * t(L) %*% RRR
        VVV = 1 / sigma2 * LtL + 1 / lambda_gamma * diag(N.Lmax)
        VVV = solve(VVV)
        gamma = matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
        MMM = MMM + L %*% gamma
        
        ##Update lambda_gamma
        lambda_gamma = 1 / stats::rgamma(1,  
                                  length(gamma) / 2 + sigma.a, 
                                  sum(gamma^2)/2 + sigma.b)
      }
      
      #Update delta
      if (!is.null(M)) {
        MMM = MMM - M %*% delta
        RRR =  Y - MMM
        XXX = 1 / sigma2 * t(M) %*% RRR
        VVV = 1 / sigma2 * MtM + 1 / lambda_delta * diag(N.Mmax)
        VVV = solve(VVV)
        delta = matrix(mvnfast::rmvn(1, VVV %*% XXX, VVV), ncol = 1)
        MMM = MMM + M %*% delta
        
        ##Update lambda_delta
        lambda_delta = 1 / stats::rgamma(1, 
                                  length(delta) / 2 + sigma.a, 
                                  sum(delta ^ 2) / 2 + sigma.b)
      }
      
      #Update residual error sigma2
      RRR = Y - MMM
      sigma2 = 1 / stats::rgamma(1, length(RRR) / 2 + sigma.a, sum(RRR ^ 2) / 2 + sigma.b)
       
      #Update spatial intercepts and parameters
      if (include.additive.annual.resid) {
         MMM = MMM - alpha_space[Z_ID]
         RRR = Y - MMM
         XXX = 1 / sigma2 * t(Gamma_space) %*% RRR
         SSS = tau_alpha * exp(-dist.space.mat / theta_alpha)
         VVV = diag(1 / sigma2 * GtG_space) + kronecker(solve(SSS), diag(N.spacetime))
         VVV = solve(VVV)
         alpha_space = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1,]
         MMM = MMM + alpha_space[Z_ID]
      
         #update tau_alpha
         SSS = t(alpha_space) %*% 
             kronecker(solve(exp(-dist.space.mat / theta_alpha)), 
                       diag(N.spacetime)) %*% 
            alpha_space
         tau_alpha = 1 / stats::rgamma(1, 
                                N.space * N.spacetime / 2 + tau.a, 
                                SSS / 2 + tau.b)
      
         #Update theta_alpha
         theta.prop = stats::rlnorm(1, 
                             log(theta_alpha), 
                             theta.tune)
         SSS.curr = tau_alpha * kronecker((exp(-dist.space.mat / theta_alpha)), 
                                          diag(N.spacetime))
         SSS.prop = tau_alpha * kronecker((exp(-dist.space.mat / theta.prop)), 
                                          diag(N.spacetime))
      
         lik.prop = mvtnorm::dmvnorm(alpha_space, 
                            rep(0, N.space * N.spacetime), 
                            SSS.curr, 
                            log = T)
         lik.curr = mvtnorm::dmvnorm(alpha_space, 
                            rep(0, N.space * N.spacetime), 
                            SSS.prop, 
                            log = T)
      
         ratio = lik.prop + 
             stats::dgamma(theta.prop, 
                    theta.a, 
                    theta.b, 
                    log = T) + 
             log(theta.prop) -
             lik.curr - stats::dgamma(theta_alpha, theta.a, theta.b, log = T) - log(theta_alpha)
         if (log(stats::runif(1)) < ratio) {
           theta_alpha = theta.prop
           theta.acc[1] = theta.acc[1] + 1
         }
      }
      
      # #Update spatial coefficent for AOD
      if (include.multiplicative.annual.resid) {
         MMM = MMM - beta_space[Z_ID] * X
         RRR = Y - MMM
         XXX = 1 / sigma2 * t(Gamma_space) %*% (X * RRR)
         SSS = tau_beta * exp(-dist.space.mat / theta_beta)
         VVV = diag(1 / sigma2 * X_S) + kronecker(solve(SSS), diag(N.spacetime))
         VVV = solve(VVV)
         beta_space = mvnfast::rmvn(1, VVV %*% XXX, VVV)[1, ]
         MMM = MMM + beta_space[Z_ID] * X
      
         #update tau_beta
         SSS = t(beta_space) %*% 
             kronecker(solve(exp(-dist.space.mat / theta_beta)), 
                       diag(N.spacetime)) %*% 
             beta_space
         tau_beta = 1 / stats::rgamma(1, N.space * N.spacetime /2 + tau.a, SSS / 2 + tau.b)
      
         #Update theta_beta
         theta.prop = stats::rlnorm(1, log(theta_beta), theta.tune)
         SSS.curr = tau_beta * kronecker((exp(-dist.space.mat/theta_beta)), 
                                         diag(N.spacetime))
         SSS.prop = tau_beta * kronecker((exp(-dist.space.mat/theta.prop)), 
                                         diag(N.spacetime))
      
         lik.prop = mvtnorm::dmvnorm(beta_space, 
                            rep(0,N.space*N.spacetime), 
                            SSS.curr, 
                            log = T)
         lik.curr = mvtnorm::dmvnorm(beta_space, 
                            rep(0, N.space * N.spacetime), 
                            SSS.prop, 
                            log = T)
      
         ratio = lik.prop + stats::dgamma(theta.prop,  
                                   theta.a, 
                                   theta.b, 
                                   log = T) + 
            log (theta.prop) -
            lik.curr - 
            stats::dgamma(theta_beta, 
                   theta.a, 
                   theta.b, 
                   log = T) - 
            log(theta_beta)

         if(log(stats::runif(1)) < ratio) {
           theta_beta = theta.prop
           theta.acc[2] = theta.acc[2] + 1
         }
       }
      
      #Update temporal intercepts and its parameters
       if (include.additive.weekly.resid) {
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
      if (include.multiplicative.weekly.resid){
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
     
     if (i > burn & i %% thin == 0){
         k = (i - burn) / thin

         #Save statistics
         LL.save[k] = sum(-2 * stats::dnorm(Y, MMM, sqrt(sigma2), log = T))
  
         alpha0.save[k] = alpha0
         beta0.save[k] = beta0
         gamma.save[k, ] = gamma
         delta.save[k, ] = delta
       
         lambda_gamma.save[k] = lambda_gamma
         lambda_delta.save[k] = lambda_delta
       
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
    delta.save = data.frame(delta.save)
    names(delta.save) = colnames(M)
    gamma.save = data.frame(gamma.save);
    names(gamma.save) = colnames(L)
    
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
                            lambda.delta = lambda_delta.save, 
                            lambda.gamma = lambda_gamma.save,
                            theta.alpha = theta_alpha.save, 
                            theta.beta = theta_beta.save, 
                            tau.alpha = tau_alpha.save, 
                            tau.beta = tau_beta.save, 
                            rho.alpha = rho_alpha.save, 
                            rho.beta = rho_beta.save, 
                            omega.alpha = omega_alpha.save, 
                            omega.beta = omega_beta.save,
                            dev = LL.save)
    
    standardize.param = rbind(data.frame(Type = "X", 
                                         Name = "AOD/CTM", 
                                         Mean = X.mean, 
                                         SD = X.sd),
                               data.frame(Type = "L", 
                                          Name = colnames(L), 
                                          Mean = L.mean, 
                                          SD = L.sd), 
                               data.frame(Type = "M", 
                                          Name = colnames (M), 
                                          Mean = M.mean, SD = M.sd))
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
         theta.acc = theta.acc)
}
