
library(MASS)
# simulation correlated data generation function -----
sim_data_cor <- function(n, p = 1000, r = 0.3, rou = 0.25,
                         phi = c(0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3),
                         alpha1.8 = c(0.4,0.4,0.5,0.5,0.5,0.5,0,0),
                         beta1.8 = c(0.4,0.5,0.5,0.6,0,0,0.5,0.5)) {
  # -------- Set the dimension of mediators M and indirect effect -------- #
  # p: number of potential mediators
  alpha <- rep(0,p)
  beta <- rep(0,p)
  
  alpha[1:8] <- alpha1.8
  beta[1:8] <- beta1.8
  # The indirect effect of the first eight mediators
  # 0.16 0.20 0.25 0.30 0.00 0.00 0.00 0.00
  
  # -------- Generate exposure X -------- #
  Sigma <- matrix(c(1,r,r,r,r,1,r,r,r,r,1,r,r,r,r,1),4,4)
  C1 <- mvrnorm(n, rep(0, 4), Sigma)              # 4 continuous confounder
  C2 <- matrix(rbinom(4*n, 1, 0.3), n, 4)         # 4 binary confounder
  COV <- cbind(C1, C2)                            # 8 confounders
  colnames(COV) <- paste0("Cov", 1:ncol(COV))
  
  phi_x <- phi  # phi: confounders --> exposure
  CX_prod <- COV %*% phi_x           # vector: product of confounders(C->X)
  
  pr <- 1/(1+exp( -CX_prod ))        # appropriate pr value
  X <- rbinom(n, 1, pr)              # Exposure X
  
  # -------- Generate mediators M -------- #
  phi_m <- phi  # phi: confounders --> mediators
  
  mu <- matrix(0,p,1)
  sigma_e <- matrix(0,p,p)  # correlation matrix
  for (i in 1:p) {
    for (j in 1:p) {
      sigma_e[i,j]=(rou^(abs(i-j)));
    }
  }
  e <- mvrnorm(n, mu, sigma_e)  # the error terms
  CM_prod <-matrix(COV %*% phi_m, n, p) 
  M <- matrix(0,n,p)
  M <- t(t(X)) %*% alpha + e + CM_prod  # mediators
  colnames(M) <- paste0("M", 1:ncol(M))
  
  # -------- Generate outcome Y -------- #
  phi_y <- phi  # phi: confounders --> outcome
  CY_prod <- COV %*% phi_y  # vector: product of confounders(C->Y)
  
  linpred <- 0.5*X + M %*% beta + CY_prod             # Linear combination of Y
  
  e1 <- rnorm(n)                     # The error term of Y
  Y <- as.numeric( 0.5 + linpred + e1 )  
  
  # -------- The data list consists of all variables -------- #
  out <- list(Y=Y, X=X, COV=COV, M=M)
  return(out)
}

# For example:
set.seed(111)
dat <- sim_data_cor(n=300, rou = 0.75)

