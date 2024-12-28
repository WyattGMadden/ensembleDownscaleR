
#' Create Cross Validate ID's
#'
#' This function creates a vector of cross validation ID's for a given number of folds and type of cross validation. 
#'
#' @inheritParams grm
#' @param num.folds Number of folds used in the cross validation process (default = 10)
#' @param type Type of cross validation to be performed. Options are "spatial", "ordinary", "spatial_clustered", or "spatial_buffered". (default = "spatial")
#' @param buffer.size Radius of buffer size, if type = "spatial_buffered" (default = 0)    
#'
#' @return A named list containing a vector of cross validation ID's, number of folds, and cross validation type. 
#'
#' @examples
#' # grm_cv()
#' 
#' 
#' @export
create_cv <- function(space.id,
                      time.id,
                      num.folds = 10,
                      type = "spatial",
                      coords = NULL,
                      buffer.size = NULL) {

    #remove first/last time observations
    time_obs_1 <- sum(time.id == 1)
    time_obs_max <- sum(time.id == max(time.id))

    space_id_cv <- space.id[(time_obs_1 + 1):(length(space.id) - time_obs_max)]
    cv_id <- rep(NA, length(space.id) - time_obs_1 - time_obs_max)

    #later overwritten if type = "spatial_buffered"
    drop_matrix <- NULL
   
    if (!is.null(coords)) {

        if (sum(colnames(coords) == c("x", "y")) != 2) {

            stop("Column names of coords must be 'x' and 'y'")

        }
        #remove first/last time observations
        coords_cv <- coords[(time_obs_1 + 1):(nrow(coords) - time_obs_max), ]


    }


    if (type == "spatial") {

        cv_id <- get_spat_cv_id(space_id_cv, num.folds)


    } else if (type == "ordinary") {

        for (i in 1:max(space_id_cv)) {

            # make sure no spatial location is not in some but not all folds
            if (sum(space_id_cv == i) < num.folds) {

                cv_id[space_id_cv == i] <- rep(0, sum(space_id_cv == i))

            } else {


                #number of observations in site i for cv
                obs_for_site_i <- sum(space_id_cv == i) 

                #unshuffled cv id's with (almost) equal number of observations in each fold
                cv_id_i <- (1:num.folds)[(1:obs_for_site_i) %% num.folds + 1]

                #shuffle cv id's
                cv_id_i <- sample(cv_id_i, replace = F)
            
                cv_id[space_id_cv == i] <- cv_id_i
            }
        
        }
            

    } else if (type == "spatial_clustered") {

        if (is.null(coords)) {
            stop("coords must be provided for spatial_clustered cross validation")
        }

        
        cv_id <- stats::kmeans(x = coords_cv, 
                               centers = num.folds)$cluster |>
            unname()



    } else if (type == "spatial_buffered") {

        if (is.null(coords)) {
            stop("coords must be provided for spatial_buffered cross validation")
        }
        if (is.null(buffer.size)) {
            stop("buffer.size must be provided for spatial_buffered cross validation")
        }

        cv_id <- get_spat_cv_id(space_id_cv, num.folds)

        locs <- unique(cbind(space_id_cv, coords_cv, cv_id))
        locs <- locs[order(locs[, 1]), ]
        dist_locs <- as.matrix(stats::dist(locs[, c("x", "y")]))

        drop_matrix <- matrix(0, 
                              nrow = length(space_id_cv), 
                              ncol = num.folds)

        for (i in 1:num.folds) {
            cv_space_id_i <- locs[, "space_id_cv"][locs[, "cv_id"] == i]
            ncv_space_id_i <- locs[, "space_id_cv"][locs[, "cv_id"] != i]
            dist_cv_ncv <- dist_locs[cv_space_id_i, ncv_space_id_i]
            spat_id_to_drop <- ncv_space_id_i[which(apply(dist_cv_ncv, 2, function(x) sum(x < buffer.size) > 0))]
            drop_matrix[, i] <- space_id_cv %in% spat_id_to_drop

        }


    } else {

        stop("type must be either 'spatial', 'ordinary', 'spatial_clustered', or spatial_buffered")

    }

    cv_id_full <- c(rep(0, time_obs_1), 
                    cv_id, 
                    rep(0, time_obs_max))

    if (type == "spatial_buffered") {
        drop_matrix = rbind(matrix(0, nrow = time_obs_1, ncol = num.folds), 
                            drop_matrix, 
                            matrix(0, nrow = time_obs_max, ncol = num.folds))
    }

    return(list(cv.id = cv_id_full, 
                num.folds = num.folds, 
                type = type,
                drop.matrix = drop_matrix))
}

#' Get Spatial Cross Validate ID's For Regular Spatial Cross Validation
#'
#' This function creates a vector of spatial cross validation ID's, used in both spatial cv and spatial_buffered cv, in main create_cv() function.
#'
#' @param space_id_cv Spatial ID's for cross validation 
#' @param num.folds Number of folds used in the cross validation process (default = 10)
#'
#' @return A vector of cross validation ID's
get_spat_cv_id <- function(space_id_cv, num.folds) {
    #unshuffled cv spatial id's with (almost) equal number of observations in each fold
    cv_spat_id <- (1:num.folds)[(1:max(space_id_cv)) %% num.folds + 1]

    #shuffle cv spatial id's
    cv_spat_id <- sample(cv_spat_id, replace = F)
    
    #assign spatial id's to cv spatial id's
    cv_id <- sapply(space_id_cv, function(i) cv_spat_id[i])
    return(cv_id)
}
