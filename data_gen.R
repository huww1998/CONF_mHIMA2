
library(MASS)
# simulation data generation function -----
sim_data <- function(n, p=10000, r=0.3, 
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
  
  phi_x <- c(0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3)  # phi: confounders --> exposure
  CX_prod <- COV %*% phi_x           # vector: product of confounders(C->X)
  
  pr <- 1/(1+exp( -CX_prod ))        # appropriate pr value
  X <- rbinom(n, 1, pr)              # Exposure X
  
  # -------- Generate mediators M -------- #
  phi_m <- c(0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3)  # phi: confounders --> mediators
  CM_prod <- COV %*% phi_m           # vector: product of confounders(C->M)
  
  ck <- t(runif(p, 0, 1))            # constant term
  M <- matrix(0, n, p)               # mediators M 
  for(i in 1:n){
    e <- rnorm(p, 0, 1.2)            # The error term of M
    M[i,] <- ck + X[i]*alpha + e + CM_prod[i]
  }
  colnames(M) <- paste0("M", 1:ncol(M))
  
  # -------- Generate outcome Y -------- #
  phi_y <- c(0.2, 0.2, 0.3, 0.3, 0.2, 0.2, 0.3, 0.3)  # phi: confounders --> outcome
  CY_prod <- COV %*% phi_y  # vector: product of confounders(C->Y)
  
  linpred <- 0.5*X + M %*% beta + CY_prod             # Linear combination of Y
  
  e1 <- rnorm(n)                     # The error term of Y
  Y <- as.numeric( 0.5 + linpred + e1 )  
  
  # -------- The data list consists of all variables -------- #
  out <- list(Y=Y, X=X, COV=COV, M=M)
  return(out)
}

# For example:
# set.seed(111)
# dat <- sim_data(n=300, p=10000)
# X <- dat$X
# Y <- dat$Y
# M <- dat$M
# COV <- dat$COV

