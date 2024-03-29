#' Make predictions for locations with two observations
#'
#' This functions uses the weights and predictions from the previous separate grm's to make ensemble predictions
#'
#' @inheritParams grm
#' @param grm.pred.1 Output from grm_pred() function for first primary variable
#' @param grm.pred.2 Output from grm_pred() function for second primary variable
#' @param date.pred.1 Date (or other unique day identifier) from first variable prediction dataset
#' @param date.pred.2 Date (or other unique day identifier) from second variable prediction dataset
#' @param weights Ensemble weights output from weight_pred 
#'
#' @return A list including ensemble predictions and standard deviations
#'
#' @examples
#' # gap_fill()
#' 
#' 
#' @export

gap_fill <- function(grm.pred.1,
                    grm.pred.2,
                    date.pred.1,
                    date.pred.2,
                    space.id, 
                    weights) {


    grm.pred.1$date <- date.pred.1
    grm.pred.2$date <- date.pred.2

    #Merge second variable predictions onto first variable predictions
    grm.pred.1$link_id <- paste(grm.pred.1$date, 
                                grm.pred.1$space.id, 
                                sep = "_")
    grm.pred.2$link_id <- paste(grm.pred.2$date, 
                                grm.pred.2$space.id, 
                                sep = "_")

    Y.pred.1 <- grm.pred.1$estimate
    Y.sd.1 <- grm.pred.1$sd

    Y.pred.2 <- grm.pred.2$estimate[match(grm.pred.1$link_id,
                                          grm.pred.2$link_id)]
    Y.sd.2 <- grm.pred.2$sd[match(grm.pred.1$link_id, 
                                  grm.pred.2$link_id)]


  
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
      
    list(ensemble.estimate = Est.ensemb,
         ensemble.sd = SD.ensemb)
}
