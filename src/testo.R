library(dplyr)
library(tidyr)
set.seed(69)
ddat <- read.csv("./dat/ddat_og.csv")
print(head(ddat))

# -------------------------------------------------------------------------------- #
# ------------------------------------- Code ------------------------------------- #
# -------------------------------------------------------------------------------- #

# Central difference formula
dfdx <- function(f, x, dx = 1e-8, n = 0, ...) {
    if (n == 0) {
        return(f(...))
    } else {
        # Recursive call to compute the (n-1)th derivative at x + dx and x - dx
        df_forward <- dfdx(f = f, x = (x + dx), dx = dx, n = (n - 1), ...)
        df_backward <- dfdx(f = f, x = (x - dx), dx = dx, n = (n - 1), ...)
        # Central difference approximation for the nth derivative
        return((df_forward - df_backward) / (2 * dx))
    }
}

# Log-likelihood function for the exponential distribution
log_likely <- function(col, theta) {
    return(length(col) * log(theta) + (theta - 1) * sum(log(col)))
}

newton_raphson_mle <- function(x, theta = 1, improve = 1e-8, max_iter = 1000) {

    x <- x[!is.na(x) & x > 0]

    for (i in 1:max_iter) {
        score <- dfdx(log_likely, x = theta, theta = theta, col = x, n = 1)
        hessian <- dfdx(log_likely, x = theta, theta = theta, col = x, n = 2)
        theta_new <- theta - score / hessian

        if (is.nan(theta_new)) stop("Theta is NaN")

        if (abs(theta_new - theta) < improve) return(theta_new)
        theta <- theta_new
    }
  
    warning("Maximum iterations reached without convergence")
    return(theta)

}

ci <- function(x, mle, conf_level = 0.95, theta_range = c(0.1, 10), tolerance = 1e-8) {

    if (is.na(mle)) return(NA_real_)
    alpha <- 1 - conf_level

    normal <- {
        hessian <- dfdx(log_likely, x, theta = mle, n = 2)
        se <- sqrt(-1 / hessian)
        z <- qnorm(1 - alpha / 2)
        lower <- mle - z * se
        upper <- mle + z * se
        return(c(lower = lower, upper = upper))
    }

    ll <- {
        max_ll <- log_likely(mle, x)
        crit_value <- qchisq(conf_level, df = 1) / 2
        root_func <- function(theta) log_likely(x, theta) - (max_ll - crit_value)
        # Find lower and upper bounds
        lower <- uniroot(root_func, c(theta_range[1], mle), tol = tolerance)$root
        upper <- uniroot(root_func, c(mle, theta_range[2]), tol = tolerance)$root
        return(c(lower = lower, upper = upper))

    }

    return(list(normal = normal, ll = ll))

}

# ddat_ci <- ddat %>%
#     summarise(across(everything(), ~ {
#         x <- .x
#         if (all(x > 0)) {
#             mle_estimate <- newton_raphson_mle(x = x)
#             confidence_intervals <- ci(x, mle_estimate)
#             confidence_intervals$ll
#         } else NA}))


# print(head(ddat_ci))


# Example usage
mle_estimate <- newton_raphson_mle(x = ddat$bmi)
confidence_intervals <- ci(ddat$bmi, mle_estimate)

print(paste("MLE estimate:", round(mle_estimate, 4)))
print("Normal approximation CI:")
print(paste("[", round(confidence_intervals$normal[1], 4), ",", round(confidence_intervals$normal[2], 4), "]"))
print("Likelihood ratio CI:")
print(paste("[", round(confidence_intervals$ll[1], 4), ",", round(confidence_intervals$ll[2], 4), "]"))
