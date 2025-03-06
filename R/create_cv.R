
#' Create Cross Validate ID's
#'
#' This function creates a vector of cross validation ID's for a given number of folds and type of cross validation. 
#'
#' @inheritParams grm
#' @param num.folds Number of folds used in the cross validation process (default = 10)
#' @param type Type of cross validation to be performed. Options are "spatial", "ordinary", "spatial_clustered", or "spatial_buffered". (default = "spatial")
#' @param buffer.size Radius of buffer size, if type = "spatial_buffered" (default = 0)    
#' @param create.from Optional cross-validation object, used to determine the cross-validation assignment for the new dataset. space.id (and time.id for "ordinary") cross-validation fold assignments will match those for the provided CV object, if present. (default = NULL)
#'
#' @return A named list containing a vector of cross validation ID's, a matrix of which observations to drop for each fold if the cv type is "spatial_buffered", and inputted objects.   
#'
#' @export
create_cv <- function(time.id,
                      space.id,
                      spacetime.id,
                      num.folds = 10L,
                      type = "ordinary",
                      coords = NULL,
                      buffer.size = NULL,
                      create.from = NULL) {
    # if not "ordinary", must supply coords
    if (type %in% c("spatial_clustered", "spatial_buffered") && is.null(coords)) {
        stop("coords must be provided for spatial_clustered, or spatial_buffered cross validation")
    }
    # space.id and time.id must be same length
    if (length(space.id) != length(time.id)) {
        stop("space.id and time.id must be the same length")
    }
    # num folds must be greater than 1 and integer
    if (!is.numeric(num.folds) || num.folds < 2) {
        stop("num.folds must be an integer greater than 1")
    }
    # coords checks
    if (!is.null(coords)) {
        if (!is.matrix(coords) & !is.data.frame(coords)) {
            stop("'coords' must be a matrix or data.frame.")
        }
    }
    # buffer.size must be numeric, and one value
    if (!is.null(buffer.size) && (!is.numeric(buffer.size) || length(buffer.size) != 1)) {
        stop("buffer.size must be a numeric value")
    }

    if (is.null(create.from)) {
        cv.output <- create_cv_original(
            time.id = time.id,
            space.id = space.id,
            spacetime.id = spacetime.id,
            num.folds = num.folds,
            type = type,
            coords = coords,
            buffer.size = buffer.size
        )
    } else {
        # error if unique space ids not shared
        if (!all(sort(unique(space.id)) == sort(unique(create.from$space.id)))) {
            stop("Unique space IDs in new data do not match those in the CV object.")
        }
        # verify create.from is cv object (list with correct names)
        if (!all(c("cv.id", "num.folds", "type", "drop.matrix", "time.id", "space.id", "coords", "buffer.size") %in% names(create.from))) {
            stop("create.from must be a cross-validation object created by create_cv()")
        }
        if (!create.from$type %in% c("spatial", "spatial_clustered", "spatial_buffered")) {
            stop("Can only use create.from for spatial-type CV (spatial, spatial_clustered, spatial_buffered).")
        }
        if (create.from$type == "spatial_buffered" & is.null(coords)) {
          stop("Must provide 'coords' for spatial_buffered re-creation.")
        }


        cv.output <- create_cv_from_previous(
            previous.cv.object = create.from,
            time.id = time.id,
            space.id = space.id,
            spacetime.id = spacetime.id,
            coords = coords
        )
    }
    return(cv.output)
}


#' Create Cross Validate ID's
#'
#' This function creates a vector of cross validation ID's for a given number of folds and type of cross validation. 
#'
#' @inheritParams grm
#' @param num.folds Number of folds used in the cross validation process (default = 10)
#' @param type Type of cross validation to be performed. Options are "spatial", "ordinary", "spatial_clustered", or "spatial_buffered". (default = "spatial")
#' @param buffer.size Radius of buffer size, if type = "spatial_buffered" (default = 0)    
#'
#' @return A named list containing a vector of cross validation ID's, a matrix of which observations to drop for each fold if the cv type is "spatial_buffered", and inputted objects.   
#'
create_cv_original <- function(time.id,
                               space.id,
                               spacetime.id,
                               num.folds = 10,
                               type = "spatial",
                               coords = NULL,
                               buffer.size = NULL) {

    #remove first/last time observations
    time_obs_1 <- sum(time.id == 1)
    time_obs_max <- sum(time.id == max(time.id))

    space_id_cv <- space.id[(time_obs_1 + 1):(length(space.id) - time_obs_max)]
    spacetime_id_cv <- spacetime.id[(time_obs_1 + 1):(length(spacetime.id) - time_obs_max)]
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

        space_spacetime_key <- paste(space_id_cv, spacetime_id_cv)
        space_spacetime_id <- unique(space_spacetime_key)

        for (i in space_spacetime_id) {

            # make sure no space/spacetime is not in some but not all folds
            if (sum(space_spacetime_key == i) < num.folds) {

                cv_id[space_spacetime_key == i] <- rep(0, sum(space_spacetime_key == i))

            } else {

                #number of observations in site i for cv
                obs_for_site_i <- sum(space_spacetime_key == i) 

                #unshuffled cv id's with (almost) equal number of observations in each fold
                cv_id_i <- (1:num.folds)[(1:obs_for_site_i) %% num.folds + 1]

                #shuffle cv id's
                cv_id_i <- sample(cv_id_i, replace = F)
            
                cv_id[space_spacetime_key == i] <- cv_id_i
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
                drop.matrix = drop_matrix,
                time.id = time.id,
                space.id = space.id,
                coords = coords,
                buffer.size = buffer.size
                ))
}


#' Create Cross Validation ID's For New Dataset Based On Previously Created Cross Validation ID's
#'
#' This function creates creates a cross-validation assignment for a new dataset, based off a previously created cross-validation assignment
#'
#' @inheritParams grm
#' @param previous.cv.object The cross-validation object created from the original dataset, used to determine the cross-validation assignment for the new dataset
#'
#' @return A named list containing a vector of cross validation ID's, a matrix of which observations to drop for each fold if the cv type is "spatial_buffered", and inputted objects.   
create_cv_from_previous <- function(previous.cv.object,
                                    time.id,
                                    space.id,
                                    spacetime.id,
                                    coords) {
  type <- previous.cv.object$type
  buffer.size <- previous.cv.object$buffer.size
  num.folds <- previous.cv.object$num.folds
  old_cv_id <- previous.cv.object$cv.id
  old_dropmat <- previous.cv.object$drop.matrix


  # remove zeros from old cv
  zero_inds <- which(old_cv_id == 0)
  old_space_id_nozero <- previous.cv.object$space.id[-zero_inds]
  old_cv_id_nozero <- old_cv_id[-zero_inds]
  if (type == "spatial_buffered") {
    old_dropmat_nozero <- old_dropmat[-zero_inds, ]
  }

  # get unique space ids and corresponding cv id's
  unique_space_id <- sort(unique(old_space_id_nozero))
  first_ocurrence_index <- match(unique_space_id, old_space_id_nozero)
  cv_id_unique <- old_cv_id_nozero[first_ocurrence_index]


  # build new cv id's
  new_inds <- match(space.id, unique_space_id)
  new_cv_id <- cv_id_unique[new_inds]
  new_dropmat <- NULL
  if (type == "spatial_buffered") {
    new_dropmat <- old_dropmat_nozero[new_inds, ]
  }

  # set first and last time observations to 0
  new_cv_id[time.id == min(time.id)] <- 0
  new_cv_id[time.id == max(time.id)] <- 0


  list(
    cv.id = new_cv_id,
    num.folds = num.folds,
    type = type,
    drop.matrix = new_dropmat,
    time.id = time.id,
    space.id = space.id,
    coords = coords,
    buffer.size = buffer.size
  )
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
