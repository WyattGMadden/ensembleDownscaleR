
#standardize covariate matrices in beginning of grm function
covariate_matrix_standardize <- function(x) {
  x <- as.matrix(x)
  x.mean <- apply(x, 2, mean)
  x.sd <- apply(x, 2, stats::sd)
  x[, x.sd != 0] <- sweep(sweep(x[, x.sd != 0], 
                                2, 
                                x.mean[x.sd != 0], 
                                "-"),
                          2,
                          x.sd[x.sd != 0],
                          "/") 
  return(list(x = x, x.mean = x.mean, x.sd = x.sd))
}

# Calculate distances 
calculate_distances <- function(coords1, coords2) {
  # Expand the coordinates in each matrix
  x1 <- matrix(rep(coords1[, "x"], 
                   each = nrow(coords2)), 
               nrow = nrow(coords1), 
               ncol = nrow(coords2))
  y1 <- matrix(rep(coords1[, "y"], 
                   each = nrow(coords2)), 
               nrow = nrow(coords1), 
               ncol = nrow(coords2))
  x2 <- matrix(rep(coords2[, "x"], 
                   times = nrow(coords1)), 
               nrow = nrow(coords1), 
               ncol = nrow(coords2))
  y2 <- matrix(rep(coords2[, "y"], 
                   times = nrow(coords1)), 
               nrow = nrow(coords1), 
               ncol = nrow(coords2))
  
  # Compute distances
  sqrt((x1 - x2)^2 + (y1 - y2)^2)
}
