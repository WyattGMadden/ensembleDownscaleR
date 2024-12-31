
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
create_cv <- function(space.id,
                      time.id,
                      num.folds = 10L,
                      type = "ordinary",
                      coords = NULL,
                      buffer.size = NULL,
                      create.from = NULL) {
    # if not "ordinary", must supply coords
    if (type != "ordinary" && is.null(coords)) {
        stop("coords must be provided for spatial, spatial_clustered, or spatial_buffered cross validation")
    }
    # space.id and time.id must be same length
    if (length(space.id) != length(time.id)) {
        stop("space.id and time.id must be the same length")
    }
    # num folds must be greater than 1 and integer
    if (!is.numeric(num.folds) || num.folds < 2) {
        stop("num.folds must be an integer greater than 1")
    }
    # buffer.size must be numeric, and one value
    if (!is.null(buffer.size) && (!is.numeric(buffer.size) || length(buffer.size) != 1)) {
        stop("buffer.size must be a numeric value")
    }

    if (is.null(create.from)) {
        cv.output <- create_cv_original(
            space.id = space.id,
            time.id = time.id,
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


        cv.output <- create_cv_from_previous(
            previous.cv.object = create.from,
            space.id = space.id,
            time.id = time.id,
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
create_cv_original <- function(space.id,
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
                                    space.id,
                                    time.id,
                                    coords = NULL) {

  type <- previous.cv.object$type
  buffer.size <- previous.cv.object$buffer.size
  num.folds <- previous.cv.object$num.folds
  old_cv_id <- previous.cv.object$cv.id
  old_dropmat <- previous.cv.object$drop.matrix  # may be NULL if not "spatial_buffered"

  df_old <- data.frame(
    space_id = previous.cv.object$space.id,
    time_id = previous.cv.object$time.id,
    cv_id = previous.cv.object$cv.id
  )
  df_new <- data.frame(
    space_id = space.id,
    time_id  = time.id
  )

  # Merge old assignments into new
  df_new <- merge(
    df_new,
    df_old,
    by    = c("space_id","time_id"),
    all.x = TRUE,
    sort  = FALSE
  )
  # df_new now has columns: space_id, time_id, cv_id (NA if new combo)

  # Only "ordinary" depends on both time and space.
  # For new combos under "ordinary", replicate the original logic:
  #   - group by space_id
  #   - if # combos < num.folds => cv_id=0
  #   - else => random shuffle of folds 1..num.folds
  if (type == "ordinary") {
    new_rows <- which(is.na(df_new$cv_id))
    if (length(new_rows) > 0) {
      site_groups <- split(new_rows, df_new$space_id[new_rows])
      for (s in names(site_groups)) {
        idx   <- site_groups[[s]]
        n_obs <- length(idx)
        if (n_obs < num.folds) {
          df_new$cv_id[idx] <- 0
        } else {
          folds_seq         <- (1:num.folds)[(1:n_obs) %% num.folds + 1]
          folds_seq         <- sample(folds_seq, replace = FALSE)
          df_new$cv_id[idx] <- folds_seq
        }
      }
    }
  } else {
    # For "spatial", "spatial_clustered", "spatial_buffered", etc.,
    # time does not affect assignment. All space IDs appear in old + new,
    # so new combos simply inherit the old space's fold if it exists,
    # which was already merged in df_new. (No new site IDs.)
  }

  # Zero out min/max time
  all_time_ids  <- df_new$time_id
  idx_min_time  <- which(all_time_ids == 1)
  idx_max_time  <- which(all_time_ids == max(all_time_ids))
  df_new$cv_id[idx_min_time] <- 0
  df_new$cv_id[idx_max_time] <- 0

  # If "spatial_buffered", recalc drop matrix for the new data
  new_dropmat <- NULL
  if (type == "spatial_buffered") {
    # Middle rows only
    rows_mid <- setdiff(seq_len(nrow(df_new)), c(idx_min_time, idx_max_time))
    if (length(rows_mid) == 0) {
      # No middle rows
      new_dropmat <- matrix(0, nrow = nrow(df_new), ncol = num.folds)
    } else {
      # Distances among unique space IDs in the middle chunk
      coords_mid   <- coords[rows_mid, , drop=FALSE]
      space_mid    <- df_new$space_id[rows_mid]
      cv_mid       <- df_new$cv_id[rows_mid]
      locs <- unique(cbind(space_mid, coords_mid, cv_mid))
      colnames(locs) <- c("space_id","x","y","cv_id_mid")
      locs <- locs[order(locs[,"space_id"]), ]

      dist_mat       <- as.matrix(stats::dist(locs[,c("x","y")]))
      new_dropmat_mid<- matrix(0, nrow = nrow(locs), ncol = num.folds)

      for (fold_i in seq_len(num.folds)) {
        fold_sites     <- locs[locs[,"cv_id_mid"] == fold_i, "space_id"]
        non_fold_sites <- locs[locs[,"cv_id_mid"] != fold_i, "space_id"]
        if (length(fold_sites)==0 || length(non_fold_sites)==0) next
        dist_cv_ncv    <- dist_mat[fold_sites, non_fold_sites, drop=FALSE]
        spat_id_to_drop<- non_fold_sites[
          apply(dist_cv_ncv, 2, function(x) sum(x < buffer.size) > 0)
        ]
        new_dropmat_mid[, fold_i] <- locs[,"space_id"] %in% spat_id_to_drop
      }

      # Expand to full data
      new_dropmat <- matrix(0, nrow = nrow(df_new), ncol = num.folds)
      site_map <- match(space_mid, locs[,"space_id"])
      for (i in seq_along(rows_mid)) {
        row_idx  <- rows_mid[i]
        site_idx <- site_map[i]
        new_dropmat[row_idx, ] <- new_dropmat_mid[site_idx, ]
      }
    }
  }

  list(
    cv.id = df_new$cv_id,
    num.folds = num.folds,
    type = type,
    drop.matrix = new_dropmat,
    time.id = df_new$time_id,
    space.id = df_new$space_id,
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
