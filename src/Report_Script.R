#STAT2300 Report

diabetes = read.csv("diabetes.csv", header = TRUE)
head(diabetes)

neg_log_likelihood = function(params, data) {
  # Extract the parameters
  beta0 = params[1]        # Intercept
  beta1 = params[2]        # Slope
  sigma = params[3]        # Standard Deviation
  
  # Extract dependent and independent variables
  x = data$bmi
  y = data$target
  
  # Predicted value of a linear model
  y_pred = beta0 + beta1 * x
  
  # Log-likelihood for normal distribution
  norm_log_likelihood = -sum(dnorm(y, mean = y_pred, sd = sigma, log = TRUE))
  
  return(norm_log_likelihood)
}
# Now to use optim() to minimise the negative log-likelihood and estimate params
int_params = c(beta0 = 0, beta1 = 0, sigma = 1)

mle = optim(
  par = int_params,           # Initial guess
  fn = neg_log_likelihood,    # The function being minimised
  data = diabetes,            # Data to be used
  method = "BFGS",            # Optimisation method (less fragile and more common)
  hessian = TRUE              # Return Hessian for variance estimation
)

# Results
mle$par

# Predicted values from the model put back into likelihood function for LP
y_pred <- mle$par[1] + mle$par[2] * diabetes$bmi

# Residuals
residuals <- diabetes$target - y_pred

# Plot residuals to check for patterns
plot(diabetes$bmi, residuals, main = "Residuals vs BMI", xlab = "BMI", ylab = "Residuals")
# Data is the same as another plot I found online
# https://rowannicholls.github.io/python/data/sklearn_datasets/diabetes.html

# beta0 being 21.86, indicates that is your progression if BMI was 0, this is
# unrealistic as no one has this, the slope being beta1 and 0.31 indicates
# every BMI increment of 1 would increase progression by 0.31 indicating a small 
# positive association which is reasonable, and standard deviation being 4102.63 
# is unrealistically high indicating that BMI alone is not enough to determine 
# the progression of diabetes based on the model or there is outlier variance
# that is just not explained by BMI. This is reasonable given the data is 10
# variables and you would assume that the other 9 are not all redundant and
# cannot be predicted by BMI alone even though there is some correlation.

