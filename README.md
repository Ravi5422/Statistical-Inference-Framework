# Statistical-Inference-Framework


# MTH210 – Statistical Computing Project
# NAME : Jani Ravi Kailash
#ROLL NO : 240486
File Used: fort.11
# 1. Introduction
In this project, I am given a dataset generated from the following probability density function:
f(x; α, λ) = αλ e^{-λx} (1 - e^{-λx})^{α - 1} , x > 0
[
]
where α > 0 and λ > 0.
The objective is to estimate the unknown parameters α and λ using the Newton-Raphson
method, and then construct 95% confidence intervals using both parametric and
non-parametric bootstrap methods.
## 2. Methodology
# 2.1 Maximum Likelihood Estimation
The likelihood function was constructed from the given PDF and converted into a log-likelihood
form for easier computation.
The Newton-Raphson method was applied using:
●
●
First derivatives (score functions)
Second derivatives (Hessian matrix)
The update rule used was:
θ
_
new = θ - H^{-1} U
where U is the gradient vector and H is the Hessian matrix.
# 2.2 Initial Guess
The initial values used were:
(α, λ) = (1, 1)
# 2.3 Stopping Criterion
The algorithm was terminated when:
max(|θ
_
new - θ|) < 10^{-6}
# 2.4 Number of Iterations
The Newton-Raphson method converged in: 6 iterations.
## 3. Bootstrap Methods
# 3.1 Parametric Bootstrap
In this method, new samples were generated using the estimated values of α and λ. The inverse
transform method was used for simulation, and parameters were re-estimated for each
generated sample.
# 3.2 Non-Parametric Bootstrap
In this case, samples were drawn with replacement from the original dataset. The parameters α
and λ were estimated for each resampled dataset.
# 3.3 Number of Replications
The number of bootstrap replications used was:
B = 1000
## 4. Results
# 4.1 MLE Estimates
α̂ = 2.0361
λ̂ = 0.906
# 4.2 Parametric Bootstrap 95% Confidence Intervals
α: (1.2971 , 4.8004)
λ: (0.6353 , 1.5228)
# 4.3 Non-Parametric Bootstrap 95% Confidence Intervals
α: (1.3819 , 4.0621)
λ: (0.5836 , 1.6993)
## 5. Conclusion
The parameters α and λ were successfully estimated using the Newton-Raphson method. The
bootstrap methods provided reliable confidence intervals for both parameters.
Both parametric and non-parametric approaches gave consistent results, indicating stability of
the estimation procedure.
## 6. Remarks
●
●
●
The Newton-Raphson method converged efficiently with the chosen initial values
Bootstrap methods helped in assessing the variability of the estimates
The entire implementation was carried out in R and ve
