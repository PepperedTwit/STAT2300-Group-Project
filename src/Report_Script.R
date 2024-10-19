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
y_pred = mle$par[1] + mle$par[2] * diabetes$bmi

# Residuals
residuals = diabetes$target - y_pred

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

# Fitting Model 1: Progression ~ BMI
model1 = lm(target ~ bmi, data = diabetes)

# Fitting Model 2: Progression ~ BMI + Age
model2 = lm(target ~ bmi + age, data = diabetes)

# Summarise the models
summary(model1)
summary(model2)

# Extract the log-likelihood for both models
logLik_model1 = logLik(model1)
logLik_model2 = logLik(model2)

# Print the log-likelihoods (optional)
logLik_model1
logLik_model2

# Calculate the likelihood ratio test statistic
LRT_stat = 2 * (logLik_model2 - logLik_model1)

# Degrees of freedom is the difference in the number of parameters
df = df.residual(model1) - df.residual(model2)

# Compute the p-value for the test
p_value = pchisq(LRT_stat, df = df, lower.tail = FALSE)

# Print the result
cat("Likelihood Ratio Test Statistic:", LRT_stat, "\n")
cat("p-value:", p_value, "\n")

# Likelihood Ratio Test Statistic: 4.413884 
# p-value: 0.03564759 
# Given the p-value is such a small number being less than 0.05, it can be
# concluded that the model with BMI and Age is a better fit than the model
# with BMI as the only constituent and the likelihood ratio test backs up
# this data and indicates a large improvement. This was almost assumed given
# the data in the set has 10 variables and would be poorly constructed if 
# so many were redundant which concludes through testing and general reasoning
# that adding additional variables will contribute to the accuracy of the model
