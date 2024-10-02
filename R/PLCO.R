#' Mixed Modeling Characterizing Genetic Effects in a Longitudinal Phenotype 
#' @details
#'
#' A package about estimating variance components via AI-ReML algorithm for a linear mixed model integrating separate genetic effects on both baseline and slope of a longitudinal response.  
#' Also, this package provides the best linear unbiased predictors for two types of effects sizes and genetic effects on intercept and slope of a logitudinal phenotype. 
#' 
#'
#'
#' @param data0 data0 has to be structured as a list of vectors 
#' @param n: A N x 1 column vector represents total number of measurements for each subject. 
#' @param Z: A N x P matrix with standardized genotypic values. 
#' @param t: A sum(n) x 1 column vector of the time variable.
#' @param y: A sum(n) x 1 column vector of the longitudinal response. 
#' @param data data has to be structured as a list of vectors and matrices.
#' @param n: A N x 1 column vector represents total number of measurements for each subject. 
#' @param Z: A N x P matrix with standardized genotypic values. 
#' @param t: A sum(n) x 1 column vector of the time variable.
#' @param y: A sum(n) x 1 column vector of the longitudinal response. 
#' @param A: A sum(n) x 2 covariate matrix of fixed effects beta_0 and beta_1.
#' @param S: A sum(n) x 4 covariate matrix of fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta).
#' @param G: A N x N matrix as P times genetic relationship matrix.
#' @param W: A sum(n) x (2P+2N) covariate matrix of random effects.
#' @param H: composite matrix for AI matrix and DL matrix.
#'
#' 
#'
#' @return
#' \item{par}{Estimated variance components via AI-ReML algorithm for a linear mixed model}
#' \item{convergence}{A binary value to indicate whether the iteration process converged or not}
#' \item{tol}{Absolute value of difference between two restricted log likelihoods }
#' \item{est.xi}{Estimated coefficients of fixed effects}
#' \item{var.est.xi}{Estimated variance-covariance matrix for Estimated coefficients of fixed effects}
#' \item{est.alpha}{Estimated effect sizes on baseline}
#' \item{est.eta}{Estimated effect sizes on slope}
#' \item{est.g.baseline}{Estimated genetic effects on baseline}
#' \item{est.g.slope}{Estimated genetic effects on slope}
#'
#'
#'
#' @importFrom stats sd 
#' @importFrom MASS mvrnorm
#' @examples
#' ## 
#' Set a specific seed for reproducibility
#' a=1
#' set.seed(a)
#' Define number of SNPs, sample size, etc.
#' P: number of SNPs
#' N: number of subjects
#' J: maximum of expected measurements per person among N individuals
#' P <- 100
#' N <- 2000
#' J <- 6
#' # theta: column vector of variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
#' theta = c(0.005, 0.002, 0.5, 0.2, 0.2)
#' # xi: coefficients of fixed effects beta_0, beta_1, mu_alpha, mu_eta
#' xi = c(0.5, 0.05, 0, 0)
#' # Call function to generate data: n,Z,t,y
#' data0 <- generate_data(P = P, N = N, J = J, theta, xi)
#' # call function to obtain information n,Z,t,y, and matrices A, S, G, W,H based on known n,Z,t,y 
#' data = generate_information(n= data0$n, Z= data0$Z, t= data0$t, y=data0$y)
#' # Perform AI-ReML algorithm for estimation of unknown variance parameters
#' # For instance, we choose an arbitrary input for AI-ReML algorithm
#' initial.par = theta
#' result <- AI_ReML(par = initial.par, l_ReML,  maxit = 1, maxtol = 1e-4, data = data, f_V, AI_DL, woodbury_inverse_V)
#' # print results from AI-ReML algorithm
#' print(result)
#' # Estimated variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
#' est.par= as.matrix(result$par, ncol=1)
#' ### estimation for coefficient of fixed effects xi
#' est_fixed_effects(est.par, data, f_V, woodbury_inverse_V)
#' # Prediction of effect sizes (alpha_p and eta_p) and genetic effects (g_i and g_i^*)
#' est_genetic_effects(est.par, data)
#' 
#'
#'
#' @references Johnson, DL and Thompson, Robin (1995) Restricted maximum likelihood estimation of variance components for univariate animal models using sparse matrix techniques and average information.  \emph{Journal of dairy science}, Vol. 78, No. 2, 449-456. 
#' @references Yang, Jian and Lee, S Hong and Goddard, Michael E and Visscher, Peter M (2011) {GCTA}: a tool for genome-wide complex trait analysis. \emph{The American Journal of Human Genetics}, Vol. 88, No. 1, 76-82.
#' @export



####################  AI-ReML algorithm for variance components #############
AI_ReML<- function(par = c(0.005, 0.002, 0.5, 0.2, 0.2) , l_ReML, maxit = 1000,  maxtol = 1e-4, data, f_V, AI_DL, woodbury_inverse_V)
{ 
  # par : initial input for unknown variance components (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e), default as c(0.005, 0.002, 0.5, 0.2, 0.2)
  # l_REML : function for computing restricted log likelihood
  # maxit : maximum of iterations, default as 1000 
  # maxtol : tolerance for convergence, default as 1e-4
  
  
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  # f_V: function for covariance matrix V 
  # AI_DL: function for computing AI and DL matrices
  # woodbury_inverse_V: function for inverse of covariance matrix V
  
  # Assuming 'par' is a vector and n0 is its length
  n0 = length(par)
  
  # initial setup
  old_par<- par
  
  
  
  # Function to update parameters
  update_parameters <- function(par) {
    AD <- AI_DL(par, data, f_V, woodbury_inverse_V)
    new_par <- par - solve(AD$AI) %*% AD$DL
    # Ensure parameters are non-negative
    new_par <- pmax(new_par, 5e-6)
    return(new_par)
  }
  
  
  # Initial log-likelihood
  old_ll <- l_ReML(old_par, data, f_V, woodbury_inverse_V) 
  
  
  # Initial update
  new_par <- update_parameters(old_par)
  new_ll <- l_ReML(new_par, data, f_V, woodbury_inverse_V)
  
  
  tol <- abs(new_ll - old_ll)
  it <- 1
  
  
  # Iterative update until convergence
  while (it < maxit && tol > maxtol) {
    old_par <- new_par
    old_ll <- new_ll
    
    new_par <- update_parameters(old_par)
    new_ll <- l_ReML(new_par, data, f_V, woodbury_inverse_V)
    
    tol <- abs(new_ll - old_ll)
    it <- it + 1
  }
  
  # Check convergence status
  convergence <- ifelse(tol < maxtol, 0, 1)
  
  
  return(list(par=new_par, convergence=convergence, tol=tol))
  
}

############################################################




##### Restricted log likelihood ######
l_ReML<-function(par, data, f_V, woodbury_inverse_V)
{
  
  # par: column vector of unknown variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  # woodbury_inverse_V: function for inverse of covariance matrix V
  
  
  
  G= data$G
  n= data$n
  t= data$t
  A= data$A
  S= data$S
  W= data$W
  y= data$y
  
  # covariance matrix V
  V =  f_V(par, data)
  
  
  ## Case1: When P is not very large, the inverse of covariance matrix V with the Woodbury Matrix Identity Formula:
  # Since P (= 253 or 200 or 100) not very large, then here we will apply the woodbury matrix identity formula as below.
  V1<- woodbury_inverse_V(par=par, data)
  
  # a1 = log(det(V))
  # a1 <- determinant(V, logarithm =TRUE)$modulus[1]
  a1 = 2*sum(log(diag(chol(V))))
  
  V2 = t(S)%*%V1%*%S
  
  # cheolsky decomosition for V2
  L <- chol(V2)
  
  # Log determinant of V2
  # a2<- determinant(V2, logarithm =TRUE)$modulus[1]
  a2 <- 2 * sum(log(diag(L)))

  # inverse of V2: 
  inverse_V2 = chol2inv(L)
  
  V3 = V1%*%S
    
  # R1<- V1 - V3%*% solve(V2)%*%t(V3) 
  R1 = V1 - V3 %*%  inverse_V2 %*% t(V3)
  
  # The Restricted log likelihood function: 
  l<- - (0.5)*( t(y) %*% R1 %*% y + a1 + a2) 
  
  l<- as.numeric(l)
  
  
  return(l)
  
  
  
}









########  Covariance Matrix V ############ 
f_V<-function(par, data)
{
  
  
  # par: column vector of unknown variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  
  G= data$G
  n= data$n
  t= data$t
  
  # total number of subjects: 
  N = length(n)
  
  # Number of unknown variance parameters 
  n0 = length(par)
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  
  
  
  V <- matrix(0, nrow= sum(n), ncol=sum(n))
  
  
  Sigma<- matrix(c(par[1], 0, 0, par[2]), nrow=2,ncol=2)
  
  
  B<- matrix(0, nrow=4, ncol=4)
  B[3,3]<- par[3]
  B[4,4]<- par[4] 
  
  
  
  
  # ith row and jth column block: 
  for (i in 1:N)
  {
    index1_i<- cumulative_sum(i,n) 
    index2_i<- cumulative_sum(i+1,n)
    
    index_set_i = (index1_i+1) : index2_i
    
    for (j in 1:N)
    {
      
      index1_j<- cumulative_sum(j,n) 
      index2_j<- cumulative_sum(j+1,n)
      
      index_set_j =  (index1_j+1) : index2_j
      
      if ( j !=i)
      {
        # Off-diagonal blocks
        A1 <- cbind(1, t[index_set_i])
        A2 <- cbind(1, t[index_set_j])
        
        V[ index_set_i,  index_set_j ] <- (G[i,j]) *A1 %*% Sigma %*% t(A2) 
      }
      else
      {
        # ith diagonal block V_i: (i>1)
        I<- diag(1,nrow=n[i],ncol=n[i])
        
        B[1:2,1:2]<-  (G[i,i])* Sigma
        
        R <- cbind(1, t[index_set_i], 1, t[index_set_i])
        
        V[ index_set_i,  index_set_i] <- par[5]* I +  (R %*% B %*% t(R))
      }
      
    }
  }  
  
  
  return(V)
  
}




#######################################################################################
## Function for first derivatives and average information related to second derivatives 
# AI_DL: function for computing AI and DL matrices
AI_DL<- function(par, data, f_V,  woodbury_inverse_V)
{
  # par: column vector of unknown variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  # f_V: function for covariance matrix V 
  # woodbury_inverse_V: function for inverse of covariance matrix V
  
  G=data$G
  t=data$t
  n=data$n
  y=data$y
  H=data$H
  A=data$A
  S=data$S
  
  # covariance matrix 
  V = f_V(par, data)
  
  ## Case1: When P is not very large, the inverse of covariance matrix V with the Woodbury Matrix Identity Formula:
  # Since P (= 253 or 200 or 100 ) not very large, then here we will apply the woodbury matrix identity formula as below.
  V1<- woodbury_inverse_V(par, data)

  V2 = V1%*%S
  
  R1<- V1 - V2%*% solve(t(S)%*%V2)%*%t(V2) 
  
  n0= length(par)
  
  AI_matrix<- matrix(0, n0, n0)
  DL_matrix<- matrix(0, n0, 1) 
  
  
  for (i in 1:n0)
  { 
    DL_matrix[i,1]<-  as.numeric(-0.5 *(sum(diag(R1%*%H[((i-1)*sum(n)+1) : (i*(sum(n))), ] )) - t(y) %*%R1 %*% H[((i-1)*sum(n)+1) : (i*(sum(n))), ] %*%R1 %*%y ))
  }
  
  for (i in 1:n0)
  {
    for (j in 1:n0)
    {
      AI_matrix[i,j]<-  as.numeric(-0.5 * t(y) %*%R1 %*% H[((i-1)*sum(n)+1) : (i*(sum(n))), ] %*%R1 %*% H[((j-1)*sum(n)+1) : (j*(sum(n))), ]%*%R1%*%y)
    }
  }
  
  return(list(AI= AI_matrix, DL=DL_matrix))
  
}














####### inverse function for the inverse of V via Woodbury Matrix Formula ####
woodbury_inverse_V<- function(par, data)
{
  # par: column vector of unknown variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  
  W= data$W
  A = data$A
  n= data$n
  Z= data$Z
  
  P = dim(Z)[2]
  N = length(n)
  
  # the inverse of matrix Q : Q1
  Q1<- matrix(0, nrow= 2*(P+N), ncol= 2*(P+N))
  
  # the inverse of matrix Sigma: Sigma1
  Sigma1<- matrix(c(par[2], 0, 0, par[1])/(par[1]*par[2]), nrow=2,ncol=2)
  
  for(p in 1:P)
  {
    Q1[(2*p-1):(2*p), (2*p-1):(2*p)]<- Sigma1
  }
  for (i in 1:N)
  {
    Q1[(2*P + i), (2*P + i)]<- 1/par[3]
    Q1[(2*P+N + i), (2*P+ N + i)]<- 1/par[4]
  }  
  
  
  V1<- diag(1/par[5],sum(n), sum(n))-  (W %*% (chol2inv(chol(Q1 + (t(W)%*%W)/par[5])) ) %*%t(W))/((par[5])^2)
  
  return(V1)
}
##################################################



### Function about estimation for coefficient of fixed effects xi
est_fixed_effects <- function(est.par, data, f_V, woodbury_inverse_V)
{
  
  # est.par : estimated variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  # f_V: function for covariance matrix V 
  # woodbury_inverse_V: function for inverse of covariance matrix V 
  
  y= data$y
  S=data$S
  
  
  
  # Estimated covariance matrix V
  V = f_V(est.par, data)
  # The inverse of V
  V1 = woodbury_inverse_V(est.par,data)
  # Estimated variance-covariance matrix of estimated xi
  
  var.est.xi <- solve(t(S)%*%V1%*%S)
  # Estimated xi
  est.xi  <- var.est.xi %*%t(S)%*%V1%*%y 
  
  return(list(est.xi= est.xi, var.est.xi = var.est.xi ))
  
  
}


# Function about prediction of effect sizes (alpha_p and eta_p) and genetic effects (g_i and g_i^*)
est_genetic_effects<-function(est.par, data)
{
  
  
  # est.par: estimated variance components = (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  
  
  ####### data : generated data with information listed as below ###################
  # G : genetic relationship matrix 
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  # H : composite matrix for AI matrix and DL matrix
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  ###################################################################
  
  W= data$W
  n= data$n
  Z= data$Z
  y=data$y
  
  P = dim(Z)[2]
  N = length(n)  
  
  
  
  # The inverse of matrix Q : Q1
  Q1<- matrix(0, nrow= 2*(P+N), ncol= 2*(P+N))
  
  # The inverse of matrix Sigma: Sigma1 where Sigma is the covariance of (alpha_p, eta_p)
  Sigma1 <- matrix(c(est.par[2], 0, 0, est.par[1])/(est.par[1]*est.par[2]), nrow=2,ncol=2)
  
  for(p in 1:P)
  {
    Q1[(2*p-1):(2*p), (2*p-1):(2*p)]<- Sigma1
  }
  
  for (i in 1:N)
  {
    Q1[(2*P + i), (2*P + i)]<- 1/est.par[3]
    Q1[(2*P+N + i), (2*P+ N + i)]<- 1/est.par[4]
  }  
  
  Q2 <- (chol2inv(chol( (Q1 + t(W)%*%W)/est.par[5] )))
  est.v<- Q2 %*%t(W)%*%y/est.par[5]
  
  est.alpha<- rep(0, P)
  est.eta <- rep(0,P)
  
  
  for (p in 1:P)
  {
    # estimated alpha's : 
    est.alpha[p]<- est.v[2*p-1]
    # estimated eta's : 
    est.eta[p]<- est.v[2*p]
  }
  
  
  est.g1<- Z%*% est.alpha
  est.g2<- Z%*% est.eta
  
  return(list(est.alpha= est.alpha, est.eta= est.eta, est.g.baseline = est.g1, est.g.slope=est.g2))
  
}






###### For instance:  
###### For a given real data set with information of n, Z, t, and y:#################### 
# Then the generated matrices can be obtained by generate_information() function 
generate_information <- function(n, Z, t, y)
{
  
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # y : the  sum(n) x 1 column vector of response 
  
  P = dim(Z)[2]
  N = length(n)
  
  # Create the covariate matrix A for fixed effects beta_0 and beta_1
  A <- matrix(1, nrow = sum(n), ncol = 2)
  A[, 2] <- t
  
  # Generate matrix G: Genetic relationship matrix
  #### Obtain the genotype information matrix G with dimension N x N ######## 
  f_G <- function(Z)
  {  
    N = dim(Z)[1]
    
    G<- matrix(0,nrow=N, ncol=N)
    
    for(i in 1:N)
    {
      for(j in 1:N)
      {
        G[i,j]<- sum(Z[i,]*Z[j,])
      }
    } 
    return(G)
  }
  ###########################################################################
  G <- f_G(Z)
  ###########################################################################
  
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  
  
  ## covariate matrix of xi = (beta_0, beta_1, mu_alpha, mu_eta)
  f_S <- function(Z,t,n)
  {
    
    # covariate matrix C: C_i = {sum_{p=1}^{P} z_{ip}}* A_i
    C<- matrix(0, nrow=sum(n), ncol=2)
    for (i in 1:N) {
      start_idx <- cumulative_sum(i, n) + 1
      end_idx <- cumulative_sum(i, n) + n[i]
      index_set =  start_idx: end_idx
      C[index_set, ] <- (sum(Z[i, ])) * A[index_set, ]
    }
    
    #  covariate matrix S = (A, C) for fixed effects: beta_0, beta_1, mu_alpha and mu_eta
    S<- matrix(0, nrow=sum(n), ncol=4)
    S[,1:2] <- A 
    S[,3:4] <- C
    
    return(S)
  }
  ###############################################################
  S <- f_S(Z, t, n)
  ###############################################################
  
  
  f_W <- function(Z, t, n)
  {  
    
    
    #### Covariate matrix of random effects v ##########
    W <-  matrix (0, nrow= sum(n) , ncol= 2*(P+N) )
    
    
    # for ith row block: (i>1)
    for (i in 1:N)
    {
      start_idx <- cumulative_sum(i, n) + 1
      end_idx <- cumulative_sum(i, n) + n[i]
      index_set =  start_idx: end_idx
      
      for(p in 1:P)
      {
        W[index_set,  (2*p-1): (2*p)]<- Z[i,p]*A[index_set, ]
      }
    }
    
    for( i in 1: N)
    {
      start_idx <- cumulative_sum(i, n) + 1
      end_idx <- cumulative_sum(i, n) + n[i]
      index_set =  start_idx: end_idx
      W[index_set, (2*P+i)] <- 1
      W[index_set, (2*P+N+i)]<- t[index_set]
    }
    
    return(W)
  }
  #####################################################
  W <- f_W(Z, t, n)
  #####################################################
  
  # Generate H matrix (assumed to handle variance structure)
  f_H<- function(Z, t, n)
  {
    
    # the number of unknown variance parameters 
    n0 <- 5
    
    #########################################
    H<- matrix(0,nrow=n0*sum(n), ncol=sum(n))
    H1<- matrix(0,nrow=sum(n), ncol=sum(n))
    H2<- matrix(0,nrow=sum(n), ncol=sum(n))
    H3<- matrix(0,nrow=sum(n), ncol=sum(n))
    H4<- matrix(0,nrow=sum(n), ncol=sum(n))
    H5<- diag(1,nrow=sum(n), ncol=sum(n))
    
    
    
    ##### H1, H2, H3, H4 ######
    # ith row and jth column block: 
    for (i in 1:N)
    {
      index1_i<- cumulative_sum(i,n) 
      index2_i<- cumulative_sum(i+1,n)
      index_set_i = (index1_i +1) : index2_i
      
      
      J1<- matrix(1,nrow=n[i],ncol=n[i])
      B1 <- t[index_set_i]
      
      # For H3 and H4: 
      # diagonal block: 
      H3[ index_set_i,  index_set_i] <- J1
      H4[ index_set_i,  index_set_i] <- B1%*%t(B1)
      
      A1<- matrix(1,nrow=n[i], ncol=1)
      
      for (j in 1:N)
      {
        
        index1_j<- cumulative_sum(j,n) 
        index2_j<- cumulative_sum(j+1,n)
        index_set_j = (index1_j +1) : index2_j
        
        # ith row and jth column block: 
        A2<- matrix(1,nrow=n[j], ncol=1)
        B2<- t[index_set_j]
        
        H1[index_set_i, index_set_j]<- (G[i,j]) *A1 %*% t(A2) 
        H2[index_set_i, index_set_j]<- (G[i,j]) *B1 %*% t(B2) 
      }
    }
    
    
    H[1: sum(n), ] <- H1
    H[(sum(n) + 1): (2*sum(n)), ]<- H2
    H[ (2*sum(n) + 1): (3*sum(n)), ]<- H3
    H[ (3*sum(n) + 1): (4*sum(n)), ]<- H4
    H[ (4*sum(n) + 1): (5*sum(n)), ]<- H5
    
    return(H) 
  }
  #####################################################
  H <- f_H(Z, t, n)
  #####################################################
  
  

 
  
  # n : the N x 1 column vector represents total number of measurements for each subjects
  # Z : standardized genotypic values matrix 
  # t : the sum(n) x 1 column vector of time effect
  # y : the  sum(n) x 1 column vector of response 
  # A : covariate matrix for fixed effects beta_0 and beta_1
  # G : genetic relationship matrix 
  # S : covariate matrix for fixed effects xi = (beta_0, beta1, mu_alpha, mu_eta)
  # W : covariate matrix for random effects
  # H : composite matrix for AI matrix and DL matrix
  
  return(list(
    n = n,
    Z = Z,
    t = t,
    y = y,
    A = A,
    S = S,
    G = G,
    W = W,
    H = H
  ))
}



# Function for generating n, Z, t, y: 
generate_data <- function(P, N, J, theta, xi) {
  # Required library for mvrnorm function
  library(MASS)
  
  # P: number of SNPs
  # N: number of subjects
  # J: maximum of measurements
  # theta: column vector of variance components (sigma^2_alpha, sigma^2_eta, sigma^2_b0, sigma^2_b1, sigma^2_e)
  # xi: coefficients of fixed effects beta_0, beta_1, mu_alpha, mu_eta
  # generate_information: a function to obtain matrices information based on known n, Z, t, y 
  
  ## Initialize fixed effects: entry age and number of measurements for ith individual
  set.seed(1)
  # Assign different probabilities for measurements from 1 to J for N individuals
  n <- sample(1:J, N, replace = TRUE, prob = c(0.01, 0.02, 0.02, 0.1, 0.15, 0.7))
  # Range of entry age from 55 to 74 years old to mimic the PLCO screening trial prostate cancer data set
  age <- sample(55:74, N, replace = TRUE)
  
  ##  Generate normalized genotypic values
  # Settings for minor allele frequency for all P SNPs 
  f0 <- runif(P, min = 0.01, max = 0.5)
  # Standardized genotypic values for pth SNP of ith individual
  Z <- matrix(0, nrow = N, ncol = P)
  for (p in 1:P) 
  {
    # Count number for pth SNP of ith individual
    x <- sample(0:2, N, replace = TRUE, prob = c((1 - f0[p])^2, 2 * f0[p] * (1 - f0[p]), f0[p]^2))
    # formula to calculate normalized genotypic value z_{ip}
    Z[, p] <- (x - 2 * f0[p]) / sqrt(2 * f0[p] * (1 - f0[p]))
  }
  
  ## Generate time effects
  t0 <- matrix(0, nrow = N, ncol = J)
  for (j in 1:J) 
  {
    # screening days for jth measurement of ith individual
    t0[, j] <- sample(((j - 1) * 365 + 1):(j * 365), N, replace = TRUE) / 365
  }
  
  
  
  ## Helper function to calculate the cumulative sum of n up to the (i-1)th element
  cumulative_sum<- function(i,n)
  {
    ifelse(i==1, 0, sum(n[1:(i-1)]) )
  }
  
  # Definition of time effect: estimated age minus 55
  t <- rep(0, sum(n))
  for (i in 1:N) {
    
    start_idx <- cumulative_sum(i, n) + 1
    end_idx <- cumulative_sum(i, n) + n[i]
    # Time effect: estimated age minus 55
    t[start_idx:end_idx] <- age[i] + sort(sample(t0[i, ], n[i])) - 55
  }
  
  
  # Simulate settings of genetic effects, random effects, and errors
  set.seed(1)
  effect_sizes <- MASS::mvrnorm(P, mu = xi[3:4], Sigma = matrix(c(theta[1], 0, 0, theta[2]), nrow = 2))
  colnames(effect_sizes) <- c("alpha", "eta")
  alpha <- effect_sizes[, 1]
  eta <- effect_sizes[, 2]
  b0 <- rnorm(N, 0, sqrt(theta[3]))
  b1 <- rnorm(N, 0, sqrt(theta[4]))
  e <- rnorm(sum(n), 0, sqrt(theta[5]))
  
  
  # Generate response variable y
  y <- matrix(0, nrow = sum(n), ncol = 1)
  g1 <- matrix(0, nrow = N, ncol = 1)
  g2 <- matrix(0, nrow = N, ncol = 1)
  
  for (i in 1:N) {
    # total genetic effects on baseline
    g1[i] <- Z[i, ] %*% alpha
    # total genetic effects on slope
    g2[i] <- Z[i, ] %*% eta
    for (j in 1:n[i]) {
      idx <- cumulative_sum(i, n) + j
      y[idx] <- xi[1] + xi[2] * t[idx] + g1[i] + g2[i] * t[idx] + b0[i] + b1[i] * t[idx] + e[idx]
    }
  }
  

  return(list(n=n, Z=Z, t=t, y=y))
  
}






