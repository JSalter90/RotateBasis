source('C:/Users/James/Dropbox/RotateBasis/R/rotation_functions.R')
library(GenSA)
library(lhs)
library(tensor)
library(far)

# Generate data from fn
n <- 60 # ensemble size
sample <- as.data.frame(2*maximinLHS(n,6) - 1)
colnames(sample) <- c("x1","x2","x3","x4","x5","x6")
output <- array(c(rep(0,100*n)), dim=c(10,10,n))
for (i in 1:n){
  output[,,i] <- fn(as.numeric(sample[i,]))
}
dim(output) <- c(100, n) # vectorising spatial output, so each column of the data matrix is a single realisation of the function

# Define observations
obs <- c(fn(c(0.7,0.01,0.01,0.25,0.8,-0.9)))

# Consider 2 different examples: 1) W = identity 2) structured W
DataBasis1 <- MakeDataBasis(data = output, RemoveMean = TRUE)
dim(DataBasis1$tBasis)
sum(DataBasis1$tBasis[,1] * DataBasis1$tBasis[,2]) # orthogonal basis vectors in L2

tmp <- t(runif(100,0.01,5))
#W <- t(tmp) %*% tmp + runif(100,0.01,0.1)*diag(100) # some positive definite 100x100 weight matrix
W <- diag(100)
diag(W) <- runif(100,0.1,10)
Winv <- GetInverse(W) # invert via chol
DataBasis2 <- MakeDataBasis(data = output, weightinv = Winv, RemoveMean = TRUE)
DataBasis2$tBasis[,1] %*% Winv %*% DataBasis2$tBasis[,2] # orthogonal basis vectors in W

# Centre observations by the ensemble mean
obsc <- obs - DataBasis1$EnsembleMean

# Look at the VarMSEplots
v1 <- VarMSEplot(DataBasis = DataBasis1, obs = obsc)
v2 <- VarMSEplot(DataBasis = DataBasis2, obs = obsc, weightinv = Winv)

# Check variances correct
svd_d <- wsvd(t(DataBasis2$CentredField), weightinv = Winv)$d
sum(svd_d[1:2]^2 / sum(svd_d^2)) # matches up to that given by VarExplained

# Project the (centred) ensemble onto the first q vectors of the SVD basis
q1 <- ExplainT(DataBasis1, vtot = 0.95)
q2 <- ExplainT(DataBasis2, vtot = 0.95, weightinv = Winv)

Coeffs1 <- CalcScores(data = DataBasis1$CentredField, basis = DataBasis1$tBasis[,1:q1])
Coeffs2 <- CalcScores(data = DataBasis2$CentredField, basis = DataBasis2$tBasis[,1:q2])
plot(Coeffs1[,4], Coeffs2[,4], pch = 18)

# Rotate!
RotatedBasis1 <- RotateBasis(DataBasis = DataBasis1, obs = obsc, kmax = 3, v = c(0.4,0.1,0.1), MaxTime = 5)
RotatedBasis2 <- RotateBasis(DataBasis = DataBasis2, obs = obsc, weightinv = Winv, kmax = 3, v = c(0.4,0.1,0.1), MaxTime = 5)
# Are resulting vectors orthogonal in W?
sum(RotatedBasis1$tBasis[,1] * RotatedBasis1$tBasis[,7])
RotatedBasis2$tBasis[,2] %*% Winv %*% RotatedBasis2$tBasis[,4]
# Do we minimise R_W?
par(mfrow=c(1,2), mar = c(4,2,2,2))
VarMSEplot(RecVarData = v1, obs = obsc, ylim = c(0,25))
VarMSEplot(DataBasis = RotatedBasis1, obs = obsc, ylim = c(0,25))

VarMSEplot(RecVarData = v2, obs = obsc, weightinv = Winv, ylim = c(0,5.5))
VarMSEplot(DataBasis = RotatedBasis2, obs = obsc, weightinv = Winv, ylim = c(0,5.5))

# Instead define some pattern, find residual basis - check orthogonal in correct norm




