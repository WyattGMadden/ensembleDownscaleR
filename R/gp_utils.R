
# covariance kernels and functions


exponential_kernel <- function(distance, theta) {
    return(exp(-distance / theta))
}


matern_kernel <- function(distance, theta, nu) {

    if (!(nu %in% c(0.5, 1.5, 2.5))) {
        stop("nu must be one of 0.5, 1.5, or 2.5")
    }
    
    root_2_nu_distance_theta <- sqrt(2 * nu) * distance / theta

    if (nu == 0.5) {

        kern <- exp(-root_2_nu_distance_theta)

    } else if (nu == 1.5) {

        kern <- (1 + root_2_nu_distance_theta) * exp(-root_2_nu_distance_theta)

    } else if (nu == 2.5) {
        
        kern <- (1 + root_2_nu_distance_theta + root_2_nu_distance_theta^2 / 3) * exp(-root_2_nu_distance_theta)

    }

    return(kern)
}





