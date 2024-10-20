library(dplyr)
library(tidyr)
library(MASS)
set.seed(69)
ddat <- read.csv("./dat/diabetes.csv")

# -------------------------------------------------------------------------------- #
# ------------------------------------- Code ------------------------------------- #
# -------------------------------------------------------------------------------- #

# Log-likelihood function for the normal distribution
log_likely <- function(col, mu) return(- (1 / (2 * var(col))) * sum((col - mu)^2))

# First derivative of a of the log-likelihood function
score <- function(col, mu) return(sum((col - mu) / var(col)))

# Second derivative of a of the log-likelihood function
hessian <- function(col) return(-length(col) / var(col))

newton_raphson_mle <- function(col, θ = 1, tolerance = 1e-8, max_iter = 1000) {
    h <- hessian(col)
    for (i in 1:max_iter) {
        s <- score(col, θ)  # Call the score function
        theta_new <- θ - solve(h, s) # Update θ
        if (max(abs(theta_new - θ)) < tolerance) return(theta_new)
        θ <- theta_new
    }
    warning("Newton-Raphson algorithm did not converge")
    return(θ)
}

conf_int <- function(col, mle, conf_level = 0.95, range = c(-1, 1), tolerance = 1e-8) {

    h <- hessian(col)
    alpha <- 1 - conf_level

    # Normal-based confidence interval
    se <- sqrt(-1 / h)
    z <- qnorm(1 - alpha / 2)

    print(mle)

    nm_min <- mle - z * se
    nm_max <- mle + z * se

    # Likelihood-based confidence interval
    max_ll <- log_likely(col, mle)
    crit_value <- qchisq(conf_level, df = 1) / 2
    root_find <- function(theta) log_likely(col, theta) - (max_ll - crit_value)
    ll_min <- uniroot(root_find, c(range[1], mle), tol = tolerance)$root
    ll_max <- uniroot(root_find, c(mle, range[2]), tol = tolerance)$root

    return(list(nm = c(min = nm_min, max = nm_max), ll = c(min = ll_min, max = ll_max)))

}


bootstrap <- function(data, lvls, α = 0.05, studentized = FALSE, debug = FALSE) {

    par_mean <- mean(data)
    lvls_len <- length(lvls)

    recurse <- function(data, lvls, α, par_mean) {

        len <- length(lvls)
        if (len == 0) return(list(μ = mean(data), std_err = sd(data)/sqrt(length(data))))

        n <- length(data)
        iter <- lvls[1]
        boot <- data.frame(
            μ = numeric(iter),
            se = numeric(iter),
            t_stat = numeric(iter)
        )

        for (i in 1:iter) {# Loop through each iteration at the current depth

            # Bootstrap resample the data
            bsample <- sample(data, n, replace = TRUE)

            # Recurse to the next depth, passing the reduced slice of levels
            result <- recurse(bsample, lvls[-1], α, par_mean)

            # Store the mean and standard error
            boot$μ[i] <- result$μ
            boot$se[i] <- result$std_err

            # Store the studentized t-statistic if requested
            if (studentized) boot$t_stat[i] <- (boot$μ[i] - par_mean) / boot$se[i]

            if (debug) cat("Iteration: ", i, "\n", sep = "")

        }

        # This only happens once at the top level
        if (len == lvls_len) return(list(μ = mean(boot$μ), lvls = lvls, α = α, par_mean = par_mean)) 
        # This happens at every other level
        else return(list(μ = mean(boot$μ), std_err = sd(boot$μ)))

    }

    return(recurse(data, lvls, α, par_mean))

}


for (name in names(ddat)) {

    if (name == "target") next
    col <- ddat[[name]]
    boot <- bootstrap(ddat$age, c(100, 100, 5), studentized = TRUE)
    mle <- boot$μ
    ci <- conf_int(col, mle)

    cat(

        "\n\nColumn:", name, "\n",
        "Bootstrap MLE: ", mle, "\n",
        "Normal CI: [", ci$nm[1], ", ", ci$nm[2], "]\n",
        "Likelihood CI: [", ci$ll[1], ", ", ci$ll[2], "]\n",
        sep = ""

    )
}