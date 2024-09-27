f<- function(a)
{
  
  # In the simulation studies, 
  # let M = 1000, set.seed(a), where a=1,2,...,M, with fixed fixed genotypic values z_{ip},
  # and random effects: alpha_p, eta_p, b0, b1, e_{ij}. 
  
  
  library("dplyr")
  library("stats")
  library("base")
  library("MASS")
  
  # Set a specific seed for reproducibility
  set.seed(1)
  # Data settings: 
  # Number of SNPs:
  P<- 100
  # Sample size: 
  N<- 2000
  # maximum of repeated measurements per person: 
  J<- 6
  
  # each person have at least one measurement:
  # individuals with different number of measurements: 
  n<- sample(1:J,N,replace = TRUE, prob= c(0.01,0.02,0.02,0.1,0.15,0.7))
  
  # entry age for each participant: 
  age<- sample(55:74, N, replace = TRUE)
  
  # True value of parameters: theta= (sigma^2_alpha, sigma^2_eta,sigma^2_b0, sigma^2_b1,  sigma^2_e)
  # Assume that There is no correlation between alpha_p and eta_p: cov(alpha_p, eta_p)=0
  # true values of variance components in the proposed model
  theta<- c(0.005, 0.002 , 0.5, 0.2, 0.2)
  
  # xi= (beta_0, beta_1, mu_alpha, mu_eta)
  # xi<- c(0.5, 0.03, 0.01, 0.005)
  xi<- c(0.5, 0.05, 0, 0)
  
  initial.par<- c(0.005, 0.002, 0.5, 0.2, 0.2)
  # initial.par<- c(0.5,0.15,0.4,0.01,0.1)
  # initial input for variance components: 
  # initial.par<- c(0.9,0.9,0.11,0.21,0.9)
  
  
  # genetic information: 
  Z<- matrix(0, nrow= N, ncol =P) 
  # minor allele frequency: 
  f0<- runif(P, min = 0.01, max = 0.5)
  # normalized genetic information: 
  for (p in 1:P)
  {
    # the number of minor allele for pth SNP:
    x<-  sample(0:2, N, replace= TRUE, prob= c((1-f0[p])^{2},  2*f0[p]*(1-f0[p]), f0[p]^{2}) ) 
    
    ## Here we will normalize the SNPs information as below: 
    Z[,p]<-   (x- 2*f0[p])/sqrt(2*f0[p]*(1-f0[p]))
  }
  
  
  
  
  # Estimated variance parameters:
  est.par<- matrix(0, nrow=5, ncol=1)
  est.xi<- matrix(0,nrow=4,ncol=1)
  
  
  # screening days: 
  t0<- matrix(0, nrow=N,ncol=J)
  for (j in 1:J)
  {
    t0[,j]<- sample( ((j-1)*365+1):(j*365) , N, replace = TRUE)/365
  }
  
  # time effect: 
  t<-matrix(0,nrow=sum(n), ncol=1)
  
  # For 1st person: 
  t[1:n[1]] <- age[1] + sort(sample(t0[1,],n[1])) -55
  
  
  # For ith person, i=2,..., N 
  for (i in 2:N)
  {
    # age randomization till jth measurements:
    t[(sum(n[1:(i-1)]) + 1) : sum(n[1:i])]<-  age[i]  +  sort(sample(t0[i,],n[i])) -55
    
  }
  
  
  
  
  
  
  
  set.seed(a)
  # In the simulation studies, let M = 1000, set.seed(a), where a=1,2,...,M, with fixed z_ip, random effects: alpha,eta,b0,b1,e_ij 
  
  # The true value of effect sizes:
  effect_sizes<- MASS::mvrnorm(P, mu=xi[3:4], Sigma= matrix(c(theta[1], 0, 0, theta[2]), nrow=2, ncol=2))
  colnames(effect_sizes)<- c("alpha", "eta")
  
  alpha<- effect_sizes[,1]
  eta<- effect_sizes[,2]
  
  
  b0<- rnorm(N, 0, sqrt(theta[3]))
  b1<- rnorm(N, 0, sqrt(theta[4]))
  e<- matrix(rnorm(sum(n),0, sd= sqrt(theta[5])), nrow= sum(n), ncol=1)
  
  
  # The response: 
  y<-matrix(0, nrow= sum(n) , ncol=1)
  
  # genetic effects: 
  g1<- matrix(0, nrow=N, ncol=1)
  g2<-  matrix(0, nrow=N, ncol=1)
  
  # For 1st person: 
  for (j in 1:n[1])
  {
    
    # PRSs:
    g1[1]<- Z[1,]%*% alpha
    g2[1]<- Z[1,]%*% eta
    y[j]<- xi[1] + xi[2]*t[j] + g1[1] + g2[1]*t[j] + b0[1]+ b1[1]*t[j] +  e[j]
  }
  
  # for ith person, i=2,..., N  
  for (i in 2:N)
  {
    for (j in 1:n[i])
    {
      
      # PRSs:
      g1[i]<- Z[i,]%*% alpha
      g2[i]<- Z[i,]%*% eta
      
      
      y[sum(n[1:(i-1)]) + j]<- xi[1] + xi[2]*t[sum(n[1:(i-1)]) + j] + g1[i] + g2[i]*t[sum(n[1:(i-1)]) + j] + b0[i]+ b1[i]*t[sum(n[1:(i-1)])+ j] +  e[sum(n[1:(i-1)])+ j]
      
    }
  }
  
  
  
  #### Obtain the genotype information matrix G with dimension N x N ######## 
  G<- matrix(0,nrow=N, ncol=N)
  
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      # Z[i,p]<- (x[i,p] - AF[p]) /sqrt(2*AF[p]*(1-AF[p]))
      # Z[j,p]<- (x[j,p] - AF[p])/ sqrt(2*AF[p]*(1-AF[p]))
      # Then, 
      # if AF[p]> 0.5, Z[i,p]<- - Z[i,p], Z[j,p]<- -Z[j,p] 
      #               => Z[i,p]*Z[j,p] will be the same 
      # Thus, we needn't change the sign of Z[i,p] when AF>0.05 
      # in fact, all we need is the product of Z[i,p] and Z[j,p]!!!!! 
      
      
      # Hence, each entry of G =  sum of Z[i,p]*Z[j,p] for all p from 1 to P
      G[i,j]<- sum(Z[i,]*Z[j,])
    }
  } 
  ###########################################################################
  
  
  
  
  
  
  # covariate matrix A: 
  A <- matrix(1, nrow= sum(n), ncol=2)
  A[,2]<- t
  
  # covariate matrix C: C_i = {sum_{p=1}^{P} z_{ip}}* A_i
  C<- matrix(0, nrow=sum(n), ncol=2)
  C[1:n[1],]<- (sum(Z[1,])) * A[1:n[1],]
  for(i in 2:N) 
  {
    C[(sum(n[1:(i-1)]) + 1) : sum(n[1:i]), ]<- (sum(Z[i,])) * A[(sum(n[1:(i-1)]) + 1) : sum(n[1:i]),]
  }
  
  ##  covariate matrix S = (A, C)
  S<- matrix(0, nrow=sum(n), ncol=4)
  S[,1:2]<-A 
  S[,3:4]<-C
  ####################################################
  
  
  
  
  
  
  #### Covariate matrix of random effects v ##########
  W <-  matrix (0, nrow= sum(n) , ncol= 2*(P+N) )
  # for 1st row with block:
  for ( p in 1 :P )
  {
    W[1:n[1],  (2*p-1): (2*p) ] <-  Z[1,p]*A[1:n[1], ]
  }
  
  
  # for ith row block: (i>1)
  for (i in 2:N)
  {
    for(p in 1:P)
    {
      W[(sum(n[1:(i-1)])+1) : sum(n[1:i]),  (2*p-1): (2*p)]<- Z[i,p]*A[(sum(n[1:(i-1)])+1) : sum(n[1:i]), ]
    }
  }
  
  
  
  W[1: n[1],  2*P+1 ] <- 1
  W[1: n[1],  2*P+1+N ] <- t[1:n[1]]
  
  
  for( i in 2: N)
  {
    W[(sum(n[1:(i-1)])+1): sum(n[1:i] ), (2*P+i)] <- 1
    W[(sum(n[1:(i-1)])+1): sum(n[1:i]), (2*P+N+i)]<- t[(sum(n[1:(i-1)])+1):sum(n[1:i])]
  }
  #####################################################
  
  
  
  
  
  # the number of unknown variance parameters 
  n0<- 5
  #########################################
  H<- matrix(0,nrow=n0*sum(n), ncol=sum(n))
  H1<- matrix(0,nrow=sum(n), ncol=sum(n))
  H2<- matrix(0,nrow=sum(n), ncol=sum(n))
  H3<- matrix(0,nrow=sum(n), ncol=sum(n))
  H4<- matrix(0,nrow=sum(n), ncol=sum(n))
  H5<- diag(1,nrow=sum(n), ncol=sum(n))
  #H6<- matrix(0,nrow=sum(n), ncol=sum(n))
  
  
  ##### H1, H2, H6 #### 
  # 1st row and jth column block: (j>1)
  for(j in 2:N) 
  {
    A1<- matrix(1,nrow=n[1], ncol=1)
    A2<- matrix(1,nrow=n[j], ncol=1)
    B1<- matrix(0,nrow=n[1], ncol=1)
    B2<- matrix(0,nrow=n[j], ncol=1)  
    B1<- t[1:n[1]]
    B2<- t[(sum(n[1:(j-1)])+1) :sum(n[1:j])]
    
    H1[ 1:n[1], (sum(n[1:(j-1)])+1) :sum(n[1:j])]<- G[1,j]* A1%*%t(A2)
    H2[ 1:n[1], (sum(n[1:(j-1)])+1) :sum(n[1:j]) ]<- G[1,j]* B1%*%t(B2)
    #H6[ 1:n[1], (sum(n[1:(j-1)])+1) :sum(n[1:j]) ]<- G[1,j]* ( B1%*%t(A2) + A1%*%t(B2) )
  }
  
  
  
  # ith row and 1st column block: (i>1)
  for(i in 2:N) 
  {
    
    A1<- matrix(1,nrow=n[i], ncol=1)
    A2<- matrix(1,nrow=n[1], ncol=1)
    B1<- matrix(0,nrow=n[i], ncol=1)
    B2<- matrix(0,nrow=n[1], ncol=1)  
    B1<- t[(sum(n[1:(i-1)])+1) :sum(n[1:i])]
    B2<- t[1:n[1]]
    
    H1[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , 1:n[1]]<- G[i,1]*A1%*%t(A2)
    
    H2[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , 1:n[1]]<- G[i,1]*B1%*%t(B2)
    
    # H6[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , 1:n[1]]<- G[i,1]*( B1%*%t(A2) + A1%*%t(B2) )
  }
  
  
  
  # ith row and jth column block: (i>1, j>1)
  for (i in 2:N)
  {
    for (j in 2:N)
    {
      # ith row and jth column block: (i>1, j>1)
      A1<- matrix(1,nrow=n[i], ncol=1)
      A2<- matrix(1,nrow=n[j], ncol=1)
      B1<- matrix(0,nrow=n[i], ncol=1)
      B2<- matrix(0,nrow=n[j], ncol=1)  
      B1<- t[ (sum(n[1:(i-1)])+1) :sum(n[1:i])  ]
      B2<- t[ (sum(n[1:(j-1)])+1) :sum(n[1:j])  ]
      
      H1[(sum(n[1:(i-1)])+1) :sum(n[1:i]), (sum(n[1:(j-1)])+1) :sum(n[1:j])]<- G[i,j] *A1 %*% t(A2) 
      H2[(sum(n[1:(i-1)])+1) :sum(n[1:i]), (sum(n[1:(j-1)])+1) :sum(n[1:j]) ]<- G[i,j] *B1 %*% t(B2) 
      #H6[(sum(n[1:(i-1)])+1) :sum(n[1:i]), (sum(n[1:(j-1)])+1) :sum(n[1:j])]<- G[i,j] *( B1%*%t(A2) + A1%*%t(B2))
      
    }
  }
  
  # For H3 and H4: 
  # diagonal block: (i>1)
  for (i in 2:N)
  {
    J1<- matrix(1,nrow=n[i],ncol=n[i])
    B1<- matrix(0, nrow=n[i], ncol=1)
    B1 <- t[ (sum(n[1:(i-1)])+1) : sum(n[1:i])]
    H3[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , (sum(n[1:(i-1)])+1) :sum(n[1:i])] <-  J1
    H4[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , (sum(n[1:(i-1)])+1) :sum(n[1:i])] <-  B1%*%t(B1)
  }    
  
  ## For H3 and H4: 
  # 1st block 
  J1<- matrix(1,nrow=n[1],ncol=n[1])
  B1<- matrix(0, nrow=n[1], ncol=1)
  B1<- t[1:n[1]]
  H3[ 1:n[1] , 1:n[1]] <-  J1
  H4[ 1:n[1] , 1:n[1]] <-  B1%*%t(B1)
  
  
  H[1: sum(n), ] <- H1
  H[(sum(n) + 1): (2*sum(n)), ]<- H2
  H[ (2*sum(n) + 1): (3*sum(n)), ]<- H3
  H[ (3*sum(n) + 1): (4*sum(n)), ]<- H4
  H[ (4*sum(n) + 1): (5*sum(n)), ]<- H5
  #H[ (5*sum(n) + 1): (6*sum(n)), ]<- H6
  #########################################
  
  
  
  
  
  ######## Compute Covariance Matrix V ############ When P is very large or not  
  f_V<-function(G,t,par)
  {
    
    V<- matrix(0, nrow= sum(n), ncol=sum(n))
    
    
    Sigma<- matrix(c(par[1], 0, 0,par[2]), nrow=2,ncol=2)
    
    B<- matrix(0, nrow=4, ncol=4)
    B[3,3]<- par[3]
    B[4,4]<- par[4] 
    
    
    # 1st row and jth column block: (j>1)
    for(j in 2:N) 
    {
      A1<- matrix(1,nrow=n[1], ncol=2)
      A2<- matrix(1,nrow=n[j], ncol=2)
      A1[,2]<-  t[1:n[1]]
      A2[,2]<-  t[(sum(n[1:(j-1)])+1) :sum(n[1:j])]
      V[ 1:n[1] , (sum(n[1:(j-1)])+1) :sum(n[1:j]) ]<- G[1,j]* A1%*%Sigma%*%t(A2)
    }
    
    
    
    # ith row and 1st column block: (i>1)
    for(i in 2:N) 
    {
      
      A1<- matrix(1,nrow=n[i], ncol=2)
      A2<- matrix(1,nrow=n[1], ncol=2)
      A1[,2]<- t[ (sum(n[1:(i-1)])+1) :sum(n[1:i])  ]
      A2[,2]<-  t[1:n[1]]
      V[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , 1:n[1]]<- G[i,1]*A1%*%Sigma%*%t(A2)
      
    }
    
    
    
    # ith row and jth column block: (i>1, j>1, i!=j)
    for (i in 2:N)
    {
      for (j in 2:N)
      {
        if ( j !=i)
        {
          # ith row and jth column block: (i>1, j>1, i!=j)
          A1<- matrix(1,nrow=n[i], ncol=2)
          A2<- matrix(1,nrow=n[j], ncol=2)
          A1[,2]<- t[ (sum(n[1:(i-1)])+1) :sum(n[1:i])  ]
          A2[,2]<- t[ (sum(n[1:(j-1)])+1) :sum(n[1:j])  ]
          V[(sum(n[1:(i-1)])+1) :sum(n[1:i]), (sum(n[1:(j-1)])+1) :sum(n[1:j]) ] <- G[i,j] *A1 %*% Sigma %*% t(A2) 
        }
        else
        {
          # ith diagonal block V_i: (i>1) 
          I<- diag(1,nrow=n[i],ncol=n[i])
          B[1:2,1:2]<-  G[i,i]* Sigma
          R<- matrix(1, nrow=n[i], ncol=4)
          R[,2] <- t[ (sum(n[1:(i-1)])+1) : sum(n[1:i])]
          R[,4]<-  t[ (sum(n[1:(i-1)])+1) : sum(n[1:i])]
          V[(sum(n[1:(i-1)])+1) :sum(n[1:i]) , (sum(n[1:(i-1)])+1) :sum(n[1:i]) ] <-  par[5]* I +  (R %*% B %*% t(R))
        }
        
      }
    }  
    
    # 1st block: V1: 
    I<- diag(1,nrow=n[1],ncol=n[1])
    B[1:2,1:2]<-  G[1,1]* Sigma
    R<- matrix(1, nrow=n[1], ncol=4)
    R[,2]<- t[1:n[1]]
    R[,4]<- t[1:n[1]]
    V[ 1:n[1] , 1:n[1]] <-  par[5]* I +  (R %*% B %*% t(R))
    
    return(V)
  }
  
  
  
  ####### V= sigma^2_e I + WQW^{T} ###### 
  ####### Method 2 to calculate covariance matrix V ### When P is not very large
  f_V_1<- function(par, W)
  {
    # covariance matrix: Q 
    Q<- matrix(0, nrow= 2*(P+N), ncol= 2*(P+N))
    
    # Sigma:
    Sigma<- matrix(c(par[1], 0, 0, par[2]), nrow=2, ncol=2)
    
    for(p in 1:P)
    {
      Q[(2*p-1):(2*p), (2*p-1):(2*p)]<- Sigma
    }
    
    for (i in 1:N)
    {
      Q[(2*P + i), (2*P + i)]<- par[3]
      Q[(2*P+N + i), (2*P+ N + i)]<- par[4]
    }  
    
    # covariance matrix V: 
    V<- diag(par[5],sum(n), sum(n))  + W%*%Q%*%t(W)
    
    return(V)
  }
  
  
  ####### inverse function for the inverse of V via Woodbury Matrix Formula ####
  woodbury_inverse<- function(par, W)
  {
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
  
  
  ##### When variance components are unknown ######
  l_ReML<-function(par)
  {
    V<- f_V(G=G, t=t,par=par)
    
    ## Case1: When P is not very large, the inverse of covariance matrix V with the Woodbury Matrix Identity Formula:
    # Since P = 253 not very large, then here we will apply the woodbury matrix identity formula as below.
    V1<- woodbury_inverse(par=par, W=W)
    
    
    ##Case2: When P is very large like 7 million, we will directly apply the below methods to compute
    # Due to the dimension of V independent with P
    # Method1: solve() function  
    #V1<- solve(V)
    
    # Method2: Apply decomposition for the inverse of V: 
    #V1<- chol2inv(chol(V))
    
    
    
    # a1 = log(det(V))
    a1<- determinant(V, logarithm =TRUE)$modulus[1]

    V2 = t(S)%*%V1%*%S
    
    V3 = V1%*%S
    
    R1<- V1 - V3%*% solve(V2)%*%t(V3) 
    
    a2<- determinant(V2, logarithm =TRUE)$modulus[1]
    
    
    # The Reml function: 
    l<- - (0.5)*( t(y) %*% R1 %*% y + a1 + a2) 
    
    l<- as.numeric(l)
    
    return(l)
  }
  
  
  
  
  
  #############################################################
  AI_DL<- function(par)
  {
    
    V<- f_V(G=G,t=t,par=par)
    
    ## Case1: When P is not very large, the inverse of covariance matrix V with the Woodbury Matrix Identity Formula:
    # Since P = 253 not very large, then here we will apply the woodbury matrix identity formula as below.
     # V1<- woodbury_inverse(par=par, W=W)
    
    ##Case2: When P is very large like 7 million, we will directly apply the below methods to compute
    # Due to the dimension of V independent with P
    # Method1: solve() function  
    V1<- solve(V)
    
    # Method2: Apply decomposition for the inverse of V: 
    #V1<- chol2inv(chol(V))
    
    V2= V1%*%S

    
    R1<- V1 - V2%*% solve(t(S)%*%V2)%*%t(V2)
    
    
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
  
  
  
  #################### adapted  AI-ReML algorithm #############
  AI_ReML<- function(par, f,maxit, maxtol)
  { 
    
    old_par<- par
    
    ## Newton's method to update unknown variance parameters: 
    AD<- AI_DL(par=old_par)
    new_par<- old_par  - (solve(AD$AI)) %*% (AD$DL)
    
    for(i in 1:n0)
    {
      if(new_par[i]< 0 || new_par[i]==0) {new_par[i]<- 5e-6}
    } 
    tol<- abs(f(new_par) - f(old_par))
    it<-1
    
    while( it < maxit && tol > maxtol ) 
    {
      
      old_par<- new_par
      ## Newton's method to updated unknown variance parameters: 
      AD<- AI_DL(par=old_par)
      new_par<- old_par  - (solve(AD$AI)) %*% (AD$DL)
      for(i in 1:n0)
      {
        if(new_par[i]< 0 || new_par[i]==0) {new_par[i]<- 5e-6}
      } 
      tol<- abs(f(new_par) - f(old_par)) 
      it<-it+1
      
    }
    
    if(tol< maxtol){convergence<-0}else{convergence <-1}  
    
    return(list(par=new_par, convergence=convergence, tol=tol))
    
  }
  
  ############################################################
  
  
  
  result<- AI_ReML(par=initial.par,    
                   f=l_ReML, #the function to maximize
                   maxit=1000, 
                   maxtol=1e-4 #  tolerance over-one step
  )
  
  
  # initial.par<- c(0.9,0.9,0.11,0.21,0.9)
  # input.par<- log(initial.par)
  # result<- optim(input.par,    
  #               neg_l_REML, # the function to minimize
  #               control=c(
  #                 maxit=500 , 
  #                  reltol=1e-8 # response tolerance over-one step
  #                ))
  # est.par<-exp(result$par)
  # convergence<- result$convergence
  
  
  
  est.par<- result$par
  convergence<- result$convergence
  df<- rep(0,14)
  df[1:5]<- est.par
  df[6]<- convergence
  
  
  ### the process of estimation for xi:
  V<- f_V(G=G,t=t, par=est.par)
  # the inverse of V: 
  V1<- woodbury_inverse(par=est.par,W=W)
  est.var <- solve(t(S)%*%V1%*%S)
  est.xi  <- est.var%*%t(S)%*%V1%*%y 
  
  
  ## estimated fixed effects: 
  df[7:10]<- est.xi
  df[11]<- est.var[1,1]
  df[12]<- est.var[2,2]
  df[13]<- est.var[3,3]
  df[14]<- est.var[4,4]
  df<- data.frame(df)
  write.csv(df, paste0(a, ".csv"))
  
}
