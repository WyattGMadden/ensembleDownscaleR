#' Make predictions for locations with two observations
#'
#' This functions uses the weights and predictions from the previous separate grm's to make ensemble predictions
#'
#' @param grm.pred.1 Output from grm_pred() function for first primary variable
#' @param grm.pred.2 Output from grm_pred() function for second primary variable
#' @param weights Output from weight_pred() or ensemble_fit(), depending on if you want ensemble predictions at full prediction locations or observation locations respectively.
#'
#' @return A list including ensemble predictions and standard deviations
#'
#' @examples
#' # gap_fill()
#' 
#' 
#' @export

gap_fill <- function(
    grm.pred.1,
    grm.pred.2,
    weights
    ) {


    ###################################
    ###         ARG CHECKS          ###
    ###################################

    # Single combined check that each of grm.pred.* is from grm_pred()
    #   i.e., must be a data frame with columns: time.id, space.id, spacetime.id, estimate, sd.
    needed_cols <- c("time.id", "space.id", "spacetime.id", "estimate", "sd")
    if (!all(sapply(list(grm.pred.1, grm.pred.2), function(df) {
        is.data.frame(df) && all(needed_cols %in% names(df))
    }))) {
        stop(
            "Both 'grm.pred.1' and 'grm.pred.2' must be data frames from grm_pred(),\n",
            "containing at least columns: 'time.id', 'space.id', 'spacetime.id', 'estimate', 'sd'."
        )
    }

    # Check 'weights' is from weight_pred() or ensemble_spatial()
    # Minimal check: must be a list with an element named 'q'.
    if (!is.list(weights) || !"q" %in% names(weights)) {
        stop("'weights' must be a list output from weight_pred() or ensemble_spatial(), containing an element named 'q'.")
    }
    if (!is.matrix(weights$q) && !is.data.frame(weights$q)) {
        stop("'weights$q' must be a matrix or data frame of posterior samples (log-odds).")
    }



    space.id <- grm.pred.1$space.id
    #Merge second variable predictions onto first variable predictions
    grm.pred.1$link_id <- paste(grm.pred.1$time.id, 
                                grm.pred.1$space.id, 
                                sep = "_")
    grm.pred.2$link_id <- paste(grm.pred.2$time.id, 
                                grm.pred.2$space.id, 
                                sep = "_")

    Y.pred.1 <- grm.pred.1$estimate
    Y.sd.1 <- grm.pred.1$sd

    Y.pred.2 <- grm.pred.2$estimate[match(grm.pred.1$link_id,
                                          grm.pred.2$link_id)]
    Y.sd.2 <- grm.pred.2$sd[match(grm.pred.1$link_id, 
                                  grm.pred.2$link_id)]


  
    weights <- weights$q
    W <- 1 / (exp(-weights) + 1)
    W.mean <-  apply(W, 1, mean)
    W.sd <-  apply(1 / (exp(-weights) + 1), 1, stats::sd)
      
    which.use <- which(!is.na(Y.pred.2))
      
    W.space <- W.mean[space.id[which.use]] 
      
    Est <- (Y.pred.1[which.use]) * W.space + Y.pred.2[which.use] * (1 - W.space)
    Est.SD <- sqrt((Y.sd.1[which.use])^2 * W.space + 
                  Y.sd.2[which.use]^2 * (1 - W.space) + 
                  (Y.pred.1[which.use]^2 * W.space + 
                   Y.pred.2[which.use]^2 * (1 - W.space) - Est^2))
      
      
    Est.ensemb <- Y.pred.1
    Est.ensemb[which.use] <- Est
    SD.ensemb <- Y.sd.1
    SD.ensemb[which.use] <- Est.SD
      
    return(
        data.frame(
            ensemble.estimate = Est.ensemb,
            ensemble.sd = SD.ensemb,
            time.id = grm.pred.1$time.id,
            space.id = grm.pred.1$space.id,
            spacetime.id = grm.pred.1$spacetime.id
            )
        )
}
