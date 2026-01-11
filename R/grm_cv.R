#' Fit the Geostatistical Regression Model (GRM) With Cross Validation
#'
#' This function fits Bayesian Hierarchical Model (BHM) in the form of Y ~ beta X + gamma L + delta M with cross-validation
#'
#' @inheritParams grm_pred
#' @inheritParams grm
#' @param cv.object A named list containing cv.id, num.folds, and type. Can be created with create_cv function. 
#' @param fit.i Fold number integer to only fit this CV fold, useful for parallelization.
#'
#' @return A data frame containing cross validation predictions
#'
#' @examples
#' # grm_cv()
#' 
#' 
#' @export
grm_cv <- function(
    Y, 
    X, 
    cv.object,
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
    nngp = F,
    number.neighbors = 10,
    discrete.theta.alpha.values = NULL,
    discrete.theta.beta.values = NULL,
    discrete.theta.gibbs = T,
    n.iter = 25000,
    burn = 5000,
    thin = 4,
    covariance = "exponential",
    covariance.kernel = NULL,
    matern.nu = 1.5,
    tau.alpha.a = 0.5,
    tau.alpha.b = 0.005,
    tau.alpha.tune = 0.2, 
    tau.beta.a = 0.5,
    tau.beta.b = 0.005,
    tau.beta.tune = 0.2, 
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
    sigma.fix.iter.num = 0,
    verbose = TRUE,
    verbose.iter = 1000,
    fit.i = NULL
    ) {

    #################################
    ###         ARG CHECKS         ###
    #################################

    # Basic checks for Y, X
    if (!is.numeric(Y) || !is.vector(Y)) {
        stop("'Y' must be a numeric vector.")
    }
    if (!is.numeric(X) || !is.vector(X)) {
        stop("'X' must be a numeric vector.")
    }
    if (length(Y) != length(X)) {
        stop("'Y' and 'X' must have the same length.")
    }
    
    # Checks for cv.object
    if (!is.list(cv.object)) {
        stop("'cv.object' must be a list.")
    }
    needed_names <- c("cv.id", "num.folds", "type")
    if (!all(needed_names %in% names(cv.object))) {
        stop("'cv.object' must be list with elements: 'cv.id', 'num.folds', and 'type'. This can be created with create_cv()")
    }
    if (!is.numeric(cv.object$num.folds) || length(cv.object$num.folds) != 1 || cv.object$num.folds < 1) {
        stop("'cv.object$num.folds' must be a positive numeric scalar.")
    }
    if (!is.character(cv.object$type) || 
        !(cv.object$type %in% c("spatial", "ordinary", "spatial_clustered", "spatial_buffered"))) {
        stop("'cv.object$type' must be one of 'spatial', 'ordinary', 'spatial_clustered', or 'spatial_buffered'.")
    }
    if (!is.numeric(cv.object$cv.id) && !is.integer(cv.object$cv.id)) {
        stop("'cv.object$cv.id' must be numeric or integer.")
    }
    if (length(cv.object$cv.id) != length(Y)) {
        stop("Length of 'cv.object$cv.id' must match length of 'Y'.")
    }
    if (cv.object$type == "spatial_buffered" && !"drop.matrix" %in% names(cv.object)) {
        stop("'cv.object$drop.matrix' is required if cv.object$type == 'spatial_buffered'.")
    }
    
    # L, M checks
    if (!is.null(L)) {
        if (!is.matrix(L) & !is.data.frame(L)) {
            stop("'L' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(L) != length(Y)) {
            stop("Number of rows in 'L' must match length of 'Y'.")
        }
    }
    if (!is.null(M)) {
        if (!is.matrix(M) & !is.data.frame(M)) {
            stop("'M' must be a matrix, data.frame, or NULL.")
        }
        if (nrow(M) != length(Y)) {
            stop("Number of rows in 'M' must match length of 'Y'.")
        }
    }
    # coords checks
    if (!is.matrix(coords) & !is.data.frame(coords)) {
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
    if (!is.numeric(number.neighbors) || length(number.neighbors) != 1 || number.neighbors < 1) {
        stop("'number.neighbors' must be a positive numeric scalar.")
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
    if (covariance == "custom" && !is.function(covariance.kernel)) {
        stop("If 'covariance' is 'custom', 'covariance.kernel' must be a valid function.")
    }
    if (covariance == "matern" && !(matern.nu %in% c(0.5, 1.5, 2.5))) {
        stop("'matern.nu' must be 0.5, 1.5, or 2.5 when 'covariance' is 'matern'.")
    }
    # prior hyperparameters 
    if (!is.numeric(sigma.a) || !is.numeric(sigma.b) || sigma.a <= 0 || sigma.b <= 0) {
        stop("'sigma.a' and 'sigma.b' must be positive numeric values.")
    }
    #sigma.fix.iter.num must be nonnegative integer
    if (!is.numeric(sigma.fix.iter.num) || sigma.fix.iter.num < 0 || sigma.fix.iter.num != round(sigma.fix.iter.num)) {
        stop("'sigma.fix.iter.num' must be a nonnegative integer.")
    }

    if (!is.numeric(tau.alpha.a) || tau.alpha.a <= 0) {
        stop("'tau.alpha.a' must be positive and numeric.")
    }
    if (!is.numeric(tau.alpha.b) || tau.alpha.b <= 0) {
        stop("'tau.alpha.b' must be positive and numeric.")
    }
    if (!is.numeric(tau.alpha.tune) || tau.alpha.tune <= 0) {
        stop("'tau.alpha.tune' must be positive and numeric.")
    }
    if (!is.numeric(tau.beta.a) || tau.beta.a <= 0) {
        stop("'tau.beta.a' must be positive and numeric.")
    }
    if (!is.numeric(tau.beta.b) || tau.beta.b <= 0) {
        stop("'tau.beta.b' must be positive and numeric.")
    }
    if (!is.numeric(tau.beta.tune) || tau.beta.tune <= 0) {
        stop("'tau.beta.tune' must be positive and numeric.")
    }
    if (!is.numeric(theta.alpha.a) || theta.alpha.a <= 0) {
        stop("'theta.alpha.a' must be positive and numeric.")
    }
    if (!is.numeric(theta.alpha.b) || theta.alpha.b <= 0) {
        stop("'theta.alpha.b' must be positive and numeric.")
    }
    if (!is.numeric(theta.alpha.tune) || theta.alpha.tune <= 0) {
        stop("'theta.alpha.tune' must be positive and numeric.")
    }
    if (!is.numeric(theta.beta.a) || theta.beta.a <= 0) {
        stop("'theta.beta.a' must be positive and numeric.")
    }
    if (!is.numeric(theta.beta.b) || theta.beta.b <= 0) {
        stop("'theta.beta.b' must be positive and numeric.")
    }
    if (!is.numeric(theta.beta.tune) || theta.beta.tune <= 0) {
        stop("'theta.beta.tune' must be positive and numeric.")
    }
    if (!is.numeric(omega.a) || omega.a <= 0) {
        stop("'omega.a' must be positive and numeric.")
    }
    if (!is.numeric(omega.b) || omega.b <= 0) {
        stop("'omega.b' must be positive and numeric.")
    }
    if (!is.logical(verbose)) {
        stop("'verbose' must be TRUE/FALSE.")
    }
    if (!is.numeric(verbose.iter) || verbose.iter < 1) {
        stop("'verbose.iter' must be a positive integer.")
    }
    # checks for initial values
    if (!is.numeric(rho.alpha.init) || rho.alpha.init <= 0 || rho.alpha.init >= 1) {
        stop("'rho.alpha.init' must be between 0 and 1 and numeric.")
    }
    if (!is.numeric(rho.beta.init) || rho.beta.init <= 0 || rho.beta.init >= 1) {
        stop("'rho.beta.init' must be between 0 and 1 and numeric.")
    }
    if (!is.numeric(theta.alpha.init) || theta.alpha.init <= 0) {
        stop("'theta.alpha.init' must be positive and numeric.")
    }
    if (!is.numeric(theta.beta.init) || theta.beta.init <= 0) {
        stop("'theta.beta.init' must be positive and numeric.")
    }

    if (min(table(space.id, spacetime.id)) < 2) {

        stop("Number of observations in each space.id/spacetime.id combination must be at least 2.")

    }


    if (!is.null(fit.i)) {
        if (!is.numeric(fit.i) || fit.i < 1 || fit.i != round(fit.i) || fit.i > cv.object$num.folds) {
            stop("'fit.i' must be an integer between 1 and cv.object$num.folds.")
        }
    }

    cv.id <- cv.object$cv.id

    Y.cv <- data.frame(time.id = time.id, 
                       space.id = space.id, 
                       spacetime.id = spacetime.id,
                       obs = Y, 
                       estimate = NA, 
                       sd = NA,
                       x = coords[, c("x")],
                       y = coords[, c("y")])
  
  
    for (cv.i in 1:cv.object$num.folds) {
    
        #fit.i override
        if (!is.null(fit.i)) cv.i <- fit.i

        print(paste0("Performing CV Experiment ---- Fold ", cv.i))

        if (cv.object$type == "spatial_buffered") {

            train.id <- cv.id != cv.i & cv.id != 0 & (!cv.object$drop.matrix[, cv.i])

        } else {

            train.id <- cv.id != cv.i & cv.id != 0

        }

        test.id.temp <- cv.id == cv.i
        #remove any test observations that are not within the training observations time range
        #these will be NA's in final cv predictions
        test.id.remove <- (min(time.id[train.id]) > time.id | max(time.id[train.id]) < time.id) & test.id.temp
        test.id <- test.id.temp & !test.id.remove

        time.id.train <- time.id[train.id]
        time.id.test <- time.id[test.id]

        Y.train <- Y[train.id]
        Y.test <- Y[test.id]

        X.train <- X[train.id]
        X.test <- X[test.id]
    
        #Subset of L matrix based on variable s
        L.train <- NULL
        L.test <- NULL
        if (!is.null(L)) {
            L <- as.matrix(L)
            L.train <- L[train.id, , drop = FALSE]
            L.test <- L[test.id, , drop = FALSE]
        }
        M.train <- NULL
        M.test <- NULL
        if (!is.null(M)) {
            M <- as.matrix(M)
            M.train <- M[train.id, , drop = FALSE]
            M.test <- M[test.id, , drop = FALSE]
        }
    
        space.id.train <- space.id[train.id]
        space.id.test <- space.id[test.id]

        #grm requires space.id to be from 1:max(space_id)
        #spatial cross validation breaks this assumption (missing space_id values)
        #here we create temporary space.id values that are from 1:length(unique(space.id))
        space.id.train.key <- sort(unique(space.id.train))
        space.id.test.key <- sort(unique(space.id.test))
        space.id.train.temp <- sapply(space.id.train, 
                                     function(x) which(space.id.train.key == x))
        space.id.test.temp <- sapply(space.id.test,
                                    function(x) which(space.id.test.key == x))
        spacetime.id.train <- spacetime.id[train.id]
        spacetime.id.test <- spacetime.id[test.id]
        coords.train <- coords[train.id, ]
        coords.test <- coords[test.id, ]
   
        fit.cv <- grm(
            Y = Y.train, 
            X = X.train, 
            L = L.train, 
            M = M.train, 
            coords = coords.train,
            space.id = space.id.train.temp, 
            time.id = time.id.train, 
            spacetime.id = spacetime.id.train, 
            include.additive.temporal.effect = include.additive.temporal.effect,
            include.additive.spatial.effect = include.additive.spatial.effect,
            include.multiplicative.temporal.effect = include.multiplicative.temporal.effect,
            include.multiplicative.spatial.effect = include.multiplicative.spatial.effect,
            nngp = nngp,
            number.neighbors = number.neighbors,
            discrete.theta.alpha.values = discrete.theta.alpha.values,
            discrete.theta.beta.values = discrete.theta.beta.values,
            discrete.theta.gibbs = discrete.theta.gibbs,
            n.iter = n.iter,
            burn = burn,
            thin = thin,
            covariance = covariance,
            covariance.kernel = covariance.kernel,
            matern.nu = matern.nu,
            tau.alpha.tune = tau.alpha.tune,
            tau.alpha.a = tau.alpha.a,
            tau.alpha.b = tau.alpha.b,
            tau.beta.tune = tau.beta.tune,
            tau.beta.a = tau.beta.a,
            tau.beta.b = tau.beta.b,
            omega.a = omega.a,
            omega.b = omega.b,
            theta.alpha.tune = theta.alpha.tune, 
            theta.alpha.a = theta.alpha.a, 
            theta.alpha.b = theta.alpha.b,
            theta.alpha.init = theta.alpha.init,
            theta.beta.tune = theta.beta.tune, 
            theta.beta.a = theta.beta.a, 
            theta.beta.b = theta.beta.b,
            theta.beta.init = theta.beta.init,
            rho.alpha.init = rho.alpha.init,
            rho.beta.init = rho.beta.init,
            sigma.a = sigma.a, 
            sigma.b = sigma.b,
            sigma.fix.iter.num = sigma.fix.iter.num,
            verbose = verbose,
            verbose.iter = verbose.iter
        )
    
        in_sample <- if (cv.object$type == "ordinary") {

            TRUE

        } else if (cv.object$type %in% c("spatial", "spatial_clustered", "spatial_buffered")) {

            FALSE

        } else {

            stop("cv.object$type must be either 'ordinary', 'spatial', 'spatial_clustered', or 'spatial_buffered'")

        }


        cv.results <- grm_pred(grm.fit = fit.cv, 
                              n.iter = (n.iter - burn) / thin,
                              X = X.test, 
                              L = L.test, 
                              M = M.test, 
                              coords = coords.test,
                              space.id = space.id.test.temp,
                              time.id = time.id.test,
                              spacetime.id = spacetime.id.test,
                              in.sample = in_sample)
                              



        Y.cv$estimate[test.id] <- cv.results$estimate
        Y.cv$sd[test.id] <- cv.results$sd

        if (!is.null(fit.i)) break

    }
 
    Y.cv$upper.95 <- Y.cv$estimate + 1.96 * Y.cv$sd
    Y.cv$lower.95 <- Y.cv$estimate - 1.96 * Y.cv$sd
  
    return(Y.cv)
}
