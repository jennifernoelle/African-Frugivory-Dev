## Simulate from the autoregressive Von Mises Fisher process

# Load packages
library(rstiefel)

## Set up simulation parameters
P <- 10 # number of nodes
R <- 5  # dim of latent eigenvalue decomposition
T.0 <- 20 # number of periods

M <- matrix(rnorm(P*R), R, P) # Modal directions matrix for initialization
rho <- 0.9
delta <- 0.9
phi <- 1
sd.e <- 0.1

## Simulate {U_t}: time series of RxP orthonormal matrices
U <- array(NA, dim = c(R,P,T.0))
U[,,1] <- U1 <- rmf.matrix(M)

for(t in 2:T.0){
  U[,,t] <- rmf.matrix(rho*U[,,t-1])
}

## Simulate {Lambda_t}: time sereis of R-dim latent
L <- matrix(NA, nrow = R, ncol = T)
L1 <- L[,1] <- rnorm(R, 0, 1)
for(t in 2:T.0){
  L[,t] <- rnorm(1, mean = delta*L[, t-1], sd = phi)
}

## Simulate {Z_t}: time series of latent connection propensities
# Z is symmetric with diagonal zero
Z <- array(NA, dim = c(P,P, T.0))
for(t in 1:T.0){
  E_t <- matrix(NA, P, P)
  e_t <- rnorm(P*(P-1)/2, 0, sd.e)
  E_t[upper.tri(E_t)] <- e_t
  E_t[lower.tri(E_t)] <- t(E_t)[lower.tri(E_t)]
  Z[,,t] <- t(U[,,t])%*% diag(L[,t])%*% U[,,t] + E_t
}

