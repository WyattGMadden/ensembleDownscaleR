
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

# Function to calculate Euclidean distances between two sets of coordinates
calculate_distances <- function(coords1, coords2) {
  # Ensure that coords1 and coords2 are matrices
  if (!is.matrix(coords1)) {
    coords1 <- as.matrix(coords1)
  }
  if (!is.matrix(coords2)) {
    coords2 <- as.matrix(coords2)
  }

  # Number of points in each set
  n1 <- nrow(coords1)
  n2 <- nrow(coords2)
  
  # Initialize the distance matrix
  distance_matrix <- matrix(0, nrow = n1, ncol = n2)
  
  # Compute distances
  for (i in 1:n1) {
    for (j in 1:n2) {
      distance_matrix[i, j] <- sqrt(sum((coords1[i, ] - coords2[j, ]) ^ 2))
    }
  }
  
  return(distance_matrix)
}
