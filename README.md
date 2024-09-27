# PLCO: Mixed Modeling Characterizing Genetic Effects in a Longitudinal Phenotype 
Approaches for estimating genetic effects at the individual level often focus on analyzing phenotypes at a single time point, with less attention given to longitudinal phenotypes. This R code introduces a mixed modeling approach that includes both random genetic and subject-specific random effects, designed to estimate genetic effects on both the baseline and the rate of change in a longitudinal trajectory. The approach tackles challenges, especially the need to account for the crossed nature of genetic and subject-specific random effects, which induce dependence among longitudinal measurements across all subjects.
 
Utilizing the Average Information Restricted Maximum Likelihood (AI-ReML) algorithm within a crossed mixed model framework for the estimation of variance components, this code provides an example to explore differential genetic effects on baseline levels and rates of change in a longitudinal phenotype by mimicing prostate-specific antigen (PSA) trajectories among participants from prostate cancer-free in the Prostate, Lung, Colorectal, and Ovarian (PLCO) Cancer Screening Trial. 

Building on the measurement error model by Thompson et al. (2008), our R code employs the AI-ReML algorithm and procedure of best linear unbiased predictors for random effects. This R code estimates:
- All variance components for genetic- and subject- specific random effects and residual error term.
- The coefficients of fixed effects.   
- The two types of effect sizes of genetic effects on baseline and slope, respectively.
- The separate genetic effects on baseline and slope, repsectively.


# Usage Examples
For real data analysis, we use PLCO.R. Among PLCO.R: 
- `generate_data` function
- `generate_information` function
- `AI-ReML` function
- `est_fixed_effects` function
- `est_genetic_effects` function
where
- The `generate_data` function mimic the PLCO datasets to provide a list of observed data: n, Z, t and y. 
- The `generate_information` function provides a list of obaserved data and required matrices: n, Z, t, y, A, S, G, W, and H.
- The `AI-ReML` function provides estimates for variance components in a linear mixed model by integrating genetic effects in a longitudinal phenotype via AI-ReML algorithm based on inputed dataset. 
- The `est_fixed_effects` function generates estimated coefficients of fixed effects and related estimated variance-covariance matric of them. 
- The `est_genetic_effects` function provides estimated effect sizes and genetic effects on both baseline and slope in the longitudinal phenotype. 
For instance: 
```r
# Set a specific seed for reproducibility
set.seed(1)
# Define number of SNPs, sample size, etc.
#  P: number of SNPs
P <- 100
#  N: number of subjects
N <- 2000
#  J: maximum of expected measurements per person among N individuals
J <- 6
# theta: column vector of variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
theta = c(0.005, 0.002, 0.5, 0.2, 0.2)
# xi: coefficients of fixed effects beta_0, beta_1, mu_alpha, mu_eta
xi = c(0.5, 0.05, 0, 0)
# Call function to generate data: n,Z,t,y
data0 <- generate_data(P = P, N = N, J = J, theta, xi)
# call function to obtain information n,Z,t,y, and matrices A, S, G, W,H based on known n,Z,t,y 
data = generate_information(n= data0$n, Z= data0$Z, t= data0$t, y=data0$y)
# Perform AI-ReML algorithm for estimation of unknown variance parameters
# For instance, we choose an arbitrary input for AI-ReML algorithm
initial.par = theta
result <- AI_ReML(par = initial.par, l_ReML,  maxit = 1, maxtol = 1e-4, data = data, f_V, AI_DL, woodbury_inverse_V)
# print results from AI-ReML algorithm
print(result)
# Estimated variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
est.par= as.matrix(result$par, ncol=1)
# estimation for coefficient of fixed effects xi
est_fixed_effects(est.par, data, f_V, woodbury_inverse_V)
# Prediction of effect sizes (alpha_p and eta_p) and genetic effects (g_i and g_i^*)
est_genetic_effects(est.par, data)
```
where 
- data0 has to be structured as a list of vectors:
- - n : A N x 1 column vector represents total number of measurements for each subject. 
- - Z: A N x P matrix with standardized genotypic values. 
- - t: A sum(n) x 1 column vector of the time variable.
- - y: A sum(n) x 1 column vector of the longitudinal response. 
- data has to be structured as a list of vectors and matrices:
- - n: A N x 1 column vector represents total number of measurements for each subject. 
- - Z: A N x P matrix with standardized genotypic values. 
- - t: A sum(n) x 1 column vector of the time variable.
- - y: A sum(n) x 1 column vector of the longitudinal response.
- - A: A sum(n) x 2 covariate matrix of fixed effects beta_0 and beta_1.
- - S: A sum(n) x 4 covariate matrix of fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta).
- - G: A N x N matrix as P times genetic relationship matrix.
- - W: A sum(n) x (2P+2N) covariate matrix of random effects.
- - H: composite matrix for AI matrix and DL matrix.

# Usage Notes
1. We recommend transforming response (e.g., using a log transformation) to approximate a normal distribution before applying our functions.
2. We recommend transforming time variable (e.g., using estimated age -55 for PLCO prostate cancer dataset) to control the scale of time variable before applying our functions. 
3. Our calculations assume there is unbalanced longitudinal data structure.
4. For sim_II.R and sim_II_genetic_effects.R, both of them are for simulation studies for a specific scenario with 1000 repetitions, we conducted them on Biowulf by setting value of a from 1 to 1000.  
   
# References
Johnson, DL and Thompson, Robin (1995) Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information.  Journal of dairy science, 78(2): 449-456. 
Yang, Jian and Lee, S Hong and Goddard, Michael E and Visscher, Peter M (2011) GCTA: a tool for genome-wide complex trait analysis. The American Journal of Human Genetics, 88(1): 76-82.
   
