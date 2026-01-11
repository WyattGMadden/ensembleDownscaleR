

euc_dist <- function(x1, y1, x2, y2) {
    sqrt((x1 - x2)^2 + (y1 - y2)^2)
}

order_coords <- function(coords, space_id) {
    coords_dist <- euc_dist(coords[, "x"], 
                           coords[, "y"], 
                           max(coords[, "x"]), 
                           max(coords[, "y"]))
    coord_ordering <- order(coords_dist)
    coords_ordered <- coords[coord_ordering, ]
    coords_dist_ordered <- coords_dist[coord_ordering]
    space_id_ordered <- space_id[coord_ordering]

    coord_reverse_ordering <- match(space_id, space_id_ordered)

    coord_list <- list("ordered_coords" = coords_ordered,
                       "coord_ordering" = coord_ordering,
                       "coord_reverse_ordering" = coord_reverse_ordering)
    return(coord_list)
}

#get neighbors for each location
get_neighbors <- function(ordered_coords, m) {
    neighbors <- list()
    for (j in 1:(nrow(ordered_coords) - 1)) {
        dist_to_i <- euc_dist(ordered_coords[j, "x"], 
                              ordered_coords[j, "y"], 
                              ordered_coords[(j + 1):nrow(ordered_coords), "x"], 
                              ordered_coords[(j + 1):nrow(ordered_coords), "y"])  
        #get the m neighbor index wrt full coords
        neighbors[[j]] <- order(dist_to_i)[1:(min(m, length(dist_to_i)))] + j
        
    }
    neighbors[[nrow(ordered_coords)]] <- numeric(0)
    return(neighbors)
}

#get locations for which each location is a neighbor
get_neighbors_inverse <- function(neighbors) {
    neighbors_inverse <- list(length(neighbors))
    for (i in 1:length(neighbors)) {
        neighbors_inverse[[i]] <- which(sapply(neighbors, function(x) i %in% x))
    }
    return(neighbors_inverse)
}

### get neighbors position
### get list of position in neighbors for each location/neighbor set
### to avoid repeated which(neighbors[[tt]] == s) in later spatial
### effects fitting
### ie list of vectors of length N.space,
### where each vector indicates zero of that location is not a neighbor
### and is the neighbor number if it is a neighbor
get_neighbor_positions <- function(neighbors){
    N.space <- length(neighbors)
    pos_in_neighbors <- vector("list", length(N.space))

    for (s in seq_len(N.space)) {
        pos <- rep(0, N.space)
        neighbor_s <- neighbors[[s]]
        if (length(neighbor_s) > 0) {
            pos[neighbor_s] <- seq_along(neighbor_s)
        }
        pos_in_neighbors[[s]] <- pos
    }
    return(pos_in_neighbors)
}

get_F_and_B <- function(kernels, tau, neighbors, kernels_partial_inverse){

    F_s <- rep(0, length(kernels))
    B_s <- list()

    for (s in 1:length(kernels)) {

        SSS_s <- tau * kernels[[s]]
        F_s[s] <- SSS_s[1, 1]
        B_s[[s]] <- rep(0, length(neighbors[[s]]))

        #if there are neighbors
        if (length(neighbors[[s]]) > 0) {

            B_s[[s]] <- as.numeric(kernels[[s]][1, -1, drop = FALSE] %*% kernels_partial_inverse[[s]])
            F_s[s] <- tau * (kernels[[s]][1,1] - B_s[[s]] %*% kernels[[s]][-1,1])

        } else {

            B_s[[s]] <- numeric(0)
            F_s[s] <- tau * kernels[[s]][1,1]

        }
    }
    return(list("B" = B_s, "F" = F_s))
}

#get neighbors from a reference set of ordered coordinates
get_neighbors_ref <- function(ordered_coords, pred_coords, m){
    neighbors_pred <- list()
    ordered_coords_temp <- ordered_coords[1:20, ]
    pred_coords_temp <- pred_coords[1:20, ]

    for (i in 1:nrow(pred_coords)) {
        dist_to_i <- euc_dist(pred_coords[i, "x"], 
                              pred_coords[i, "y"], 
                              ordered_coords[, "x"], 
                              ordered_coords[, "y"])  
        #get the m neighbor index wrt full coords
        neighbors_pred[[i]] <- order(dist_to_i)[1:(min(m, length(dist_to_i)))]
    }

    return(neighbors_pred)
}

get_dist_matrices <- function(coords, neighbors) {
    dist_matrices <- list()
    for (i in 1:nrow(coords)) {
        neighbor_i <- neighbors[[i]]
        neighbor_coords <- coords[c(i, neighbor_i), , drop = FALSE]
        dist_matrices[[i]] <- as.matrix(stats::dist(neighbor_coords, 
                                                    upper = TRUE, 
                                                    diag = TRUE))
    }
    return(dist_matrices)
}
get_dist_matrices_ref <- function(ordered_coords, coords_pred, neighbors) {
    dist_matrices <- list()
    for (i in 1:nrow(coords_pred)) {
        neighbor_i <- neighbors[[i]]
        neighbor_coords <- rbind(coords_pred[i, , drop = FALSE], 
                                 ordered_coords[neighbor_i, , drop = FALSE])
        dist_matrices[[i]] <- as.matrix(stats::dist(neighbor_coords, 
                                                    upper = TRUE, 
                                                    diag = TRUE))
    }
    return(dist_matrices)
}
rnngp <- function(ordered_coords, neighbors, phi, r, cov_kern) {
    y <- rep(0, nrow(ordered_coords))
    y[nrow(ordered_coords)] <- stats::rnorm(1, 0, sqrt(phi * cov_kern(distance = 0, 
                                                                      theta = r)))
    for (i in (nrow(ordered_coords) - 1):1) {
        neighbor_i <- neighbors[[i]]
        neighbor_coords <- ordered_coords[neighbor_i, , drop = FALSE]
        neighbor_y <- y[neighbor_i]
        y_cov <- phi * cov_kern(distance = 0, theta = r)
        cross_cov <- phi * cov_kern(distance = euc_dist(neighbor_coords[, "x"], 
                                                        neighbor_coords[, "y"], 
                                                        ordered_coords[i, "x"], 
                                                        ordered_coords[i, "y"]),
                                    theta = r)
        neighbor_cov <- phi * cov_kern(distance = as.matrix(stats::dist(neighbor_coords, 
                                                                        upper = TRUE, 
                                                                        diag = TRUE)),
                                       theta = r)
        cross_neighbor_cov <- t(cross_cov) %*% solve(neighbor_cov) 
        conditional_mean <- cross_neighbor_cov %*% neighbor_y
        conditional_cov <- y_cov - cross_neighbor_cov %*% cross_cov
        y[i] <- stats::rnorm(1, 
                             conditional_mean,
                             sqrt(conditional_cov))
    }
    return(y)
}



dnngp <- function(y, neighbors, dist_matrices, phi, r, cov_kern, log = FALSE) {
    d <- rep(0, length(y))
    d[length(y)] <- stats::dnorm(y[length(y)], 
                                            0, 
                                            sqrt(phi * cov_kern(distance = 0,
                                                                theta = r)), 
                                            log = log)
    for (i in 1:(length(y) - 1)) {
        neighbor_i <- neighbors[[i]]
        dist_space_mat_i <- dist_matrices[[i]]

        joint_cov <- phi * cov_kern(distance = dist_space_mat_i, 
                                    theta = r)

        neighbor_y <- y[neighbor_i]
        cross_neighbor_cov <- t(joint_cov[1, -1]) %*% solve(joint_cov[-1, -1])

        conditional_mean <- cross_neighbor_cov %*% neighbor_y
        conditional_cov <- joint_cov[1, 1] - cross_neighbor_cov %*% joint_cov[1, -1]
        d[i] <- stats::dnorm(y[i], 
                             conditional_mean,
                             sqrt(conditional_cov),
                             log = log)
    }
    return(d)
}

dnngp_discrete_theta <- function(y, neighbors, dist_matrices, phi, which_theta, kerns, kerns_partial_inv, log = FALSE) {
    d <- rep(0, length(y))
    d[length(y)] <- stats::dnorm(y[length(y)], 
                                            0, 
                                            sqrt(phi * kerns[[length(kerns)]]),
                                            log = log)
    for (i in 1:(length(y) - 1)) {
        neighbor_i <- neighbors[[i]]
        dist_space_mat_i <- dist_matrices[[i]]

        joint_cov <- phi * kerns[[i]]

        neighbor_y <- y[neighbor_i]
        cross_neighbor_cov <- t(joint_cov[1, -1]) %*% ((1 / phi) * kerns_partial_inv[[i]])

        conditional_mean <- cross_neighbor_cov %*% neighbor_y
        conditional_cov <- joint_cov[1, 1] - cross_neighbor_cov %*% joint_cov[1, -1]
        d[i] <- stats::dnorm(y[i], 
                             conditional_mean,
                             sqrt(conditional_cov),
                             log = log)
    }
    return(d)
}



mcmc_draw_spatial_nngp <- function(
    spatial_effect, 
    sigma2,
    neighbors,
    neighbors_inverse,
    pos_in_neighbors,
    B_s,
    F_s,
    covariate_sums,
    resid_sums
    ){

    for (s in 1:length(spatial_effect)) {

        sum_B_F_inv_B <- 0
        sum_B_F_inv_a <- 0
        for (tt in neighbors_inverse[[s]]) {
            B_tt <- B_s[[tt]]
            B_tts <- B_tt[pos_in_neighbors[[tt]][s]]
            F_tt <- F_s[tt]

            a_tt_s <- spatial_effect[tt]
            for (l in neighbors[[tt]]) {
                if (l != s) {
                    a_tt_s <- a_tt_s - B_tt[pos_in_neighbors[[tt]][l]] * spatial_effect[l]
                }
            }

            sum_B_F_inv_B <- sum_B_F_inv_B + B_tts^2 / F_tt
            sum_B_F_inv_a <- sum_B_F_inv_a + B_tts * (1 / F_tt) * a_tt_s
        }

        V_s <- (1 / sigma2 * covariate_sums[s] + (1 / F_s[s]) + sum_B_F_inv_B) ^ (-1)
        mu_s <- resid_sums[s] + sum_B_F_inv_a
        if (length(neighbors[[s]]) > 0) {
            mu_s <- mu_s + (1 / F_s[s]) * B_s[[s]] %*% spatial_effect[neighbors[[s]]]
        }

        spatial_effect[s] <- stats::rnorm(1, V_s * mu_s, sqrt(V_s))

    }

    return(spatial_effect)

}

