#  MTH210 Statistical Computing Project

#########################################
# NAME      :   Jani Ravi Kailash
# Roll no.  :   240486
# Data File :   fort.11
##########################################

# Observed sample values
x <- c(1.417, 5.066, 0.413, 2.688, 5.922, 1.473, 0.517, 1.047, 0.878, 2.755,
       2.417, 1.270, 0.732, 0.314, 2.694, 0.697, 2.600, 1.025, 0.419, 0.602,
       1.174, 0.962, 3.025, 0.723, 1.314)

n <- length(x)   # number of observations in the sample

# Newton-Raphson method to estimate parameters (alpha, lambda)
nr_method <- function(data, init = c(1, 1), tol = 1e-6, max_iter = 100) {
  
  # Initial guess is taken as (1,1)
  # This is a reasonable starting point since both parameters must be positive
  theta <- init   # start from initial values
  
  for (iter in 1:max_iter) {
    
    alpha <- theta[1]
    lambda <- theta[2]
    
    # Compute exponential term once to avoid repetition
    exp_val <- exp(-lambda * data)
    
    # These expressions appear repeatedly in derivatives
    t1 <- (data * exp_val) / (1 - exp_val)
    t2 <- (data^2 * exp_val) / ((1 - exp_val)^2)
    
    # Score functions (tell us direction of improvement)
    d_alpha <- (n / alpha) + sum(log(1 - exp_val))
    d_lambda <- (n / lambda) - sum(data) + (alpha - 1) * sum(t1)
    
    score <- c(d_alpha, d_lambda)
    
    # Hessian matrix gives curvature (used for better updates)
    h11 <- -n / (alpha^2)
    h12 <- sum(t1)
    h22 <- -n / (lambda^2) - (alpha - 1) * sum(t2)
    
    H <- matrix(c(h11, h12, h12, h22), 2, 2)
    
    # Small diagonal is added to avoid numerical issues (near-singular matrix)
    theta_new <- as.numeric(theta - solve(H + diag(1e-6, 2)) %*% score)
    
    # Ensure parameters remain positive after update
    theta_new <- pmax(theta_new, 1e-6)
    
    # Stop when change between iterations becomes very small
    if (max(abs(theta_new - theta)) < tol) {
      return(list(est = theta_new, iter = iter))
    }
    
    theta <- theta_new   # move to next iteration
  }
  
  # If not converged, return last computed value
  return(list(est = theta, iter = max_iter))
}

# Compute MLE using NR method
mle <- nr_method(x, c(1, 1))
alpha_hat <- mle$est[1]
lambda_hat <- mle$est[2]

# Function to simulate data from fitted distribution
# Uses inverse transform sampling
rgenexp <- function(n, alpha, lambda) {
  u <- runif(n)
  -log(1 - u^(1/alpha)) / lambda
}

# Number of bootstrap samples
B <- 1000

# Vectors to store bootstrap estimates
alpha_param <- numeric(B)
lambda_param <- numeric(B)
alpha_nonparam <- numeric(B)
lambda_nonparam <- numeric(B)

set.seed(210)   # ensures same results every time code is run

for (i in 1:B) {
  
  # Parametric bootstrap: simulate data using estimated parameters
  x_param <- rgenexp(n, alpha_hat, lambda_hat)
  fit_param <- try(nr_method(x_param, c(alpha_hat, lambda_hat)), silent = TRUE)
  
  # Store results only if estimation succeeds
  if (class(fit_param) != "try-error") {
    alpha_param[i] <- fit_param$est[1]
    lambda_param[i] <- fit_param$est[2]
  }
  
  # Non-parametric bootstrap: resample from original data
  x_np <- sample(x, n, replace = TRUE)
  fit_np <- try(nr_method(x_np, c(alpha_hat, lambda_hat)), silent = TRUE)
  
  if (class(fit_np) != "try-error") {
    alpha_nonparam[i] <- fit_np$est[1]
    lambda_nonparam[i] <- fit_np$est[2]
  }
}

# Remove failed cases before computing confidence intervals
alpha_param <- na.omit(alpha_param)
lambda_param <- na.omit(lambda_param)
alpha_nonparam <- na.omit(alpha_nonparam)
lambda_nonparam <- na.omit(lambda_nonparam)

# Final output
cat("MLE Estimates:\n")
cat("Alpha =", round(alpha_hat, 4), "\n")
cat("Lambda =", round(lambda_hat, 4), "\n")
cat("Iterations =", mle$iter, "\n\n")

cat("Parametric Bootstrap 95% CI:\n")
cat("Alpha =", round(quantile(alpha_param, c(0.025, 0.975)), 4), "\n")
cat("Lambda =", round(quantile(lambda_param, c(0.025, 0.975)), 4), "\n\n")

cat("Non-Parametric Bootstrap 95% CI:\n")
cat("Alpha =", round(quantile(alpha_nonparam, c(0.025, 0.975)), 4), "\n")
cat("Lambda =", round(quantile(lambda_nonparam, c(0.025, 0.975)), 4), "\n")