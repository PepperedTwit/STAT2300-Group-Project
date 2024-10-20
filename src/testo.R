set.seed(69)
ddat <- read.csv("./dat/diabetes.csv")
# head(ddat)

# -------------------------------------------------------------------------------- #
# ------------------------------------- Code ------------------------------------- #
# -------------------------------------------------------------------------------- #

# # Central difference formula
# dfdx <- function(f, x, dx = 1e-8, n = 0) {
#     if (n == 0) {
#         return(f(x))
#     } else {
#         # Recursive call to compute the (n-1)th derivative at x + dx and x - dx
#         df_forward <- dfdx(f, x + dx, dx, n - 1)
#         df_backward <- dfdx(f, x - dx, dx, n - 1)
#         # Central difference approximation for the nth derivative
#         return((df_forward - df_backward) / (2 * dx))
#     }
# }

# # Log-likelihood function for the exponential distribution
# log_likely <- function(theta, x) {
#     return(length(x) * log(theta) + (theta - 1) * sum(log(x)))
# }

# newton_raphson_mle <- function(x, theta = 1, improve = 1e-8, max_iter = 1000) {

#     for (i in 1:max_iter) {
#         score <- dfdx(log_likely, theta, n = 1)
#         hessian <- dfdx(log_likely, theta, n = 2)
#         theta_new <- theta - score / hessian
#         if (abs(theta_new - theta) < improve) return(theta_new)
#         theta <- theta_new
#     }
  
#     warning("Maximum iterations reached without convergence")
#     return(theta)

# }

# ci <- function(x, mle, conf_level = 0.95, theta_range = c(0.1, 10), tolerance = 1e-8) {

#     alpha <- 1 - conf_level

#     normal <- {
#         hessian <- dfdx(log_likely, mle, n = 2)
#         se <- sqrt(-1 / hessian)
#         z <- qnorm(1 - alpha / 2)
#         lower <- mle - z * se
#         upper <- mle + z * se
#         return(c(lower = lower, upper = upper))
#     }

#     ll <- {
#         max_ll <- log_likely(mle, x)
#         crit_value <- qchisq(conf_level, df = 1) / 2
#         root_func <- function(theta) log_likely(theta, x) - (max_ll - crit_value)
#         # Find lower and upper bounds
#         lower <- uniroot(root_func, c(theta_range[1], mle), tol = tolerance)$root
#         upper <- uniroot(root_func, c(mle, theta_range[2]), tol = tolerance)$root
#         return(c(lower = lower, upper = upper))

#     }

#     return(list(normal = normal, ll = ll))

# }

# Example usage
# mle_estimate <- newton_raphson_mle(ddat)
# confidence_intervals <- ci(ddat, mle_estimate)

# print(paste("MLE estimate:", round(mle_estimate, 4)))
# print("Normal approximation CI:")
# print(paste("[", round(confidence_intervals$normal[1], 4), ",", round(confidence_intervals$normal[2], 4), "]"))
# print("Likelihood ratio CI:")
# print(paste("[", round(confidence_intervals$ll[1], 4), ",", round(confidence_intervals$ll[2], 4), "]"))
