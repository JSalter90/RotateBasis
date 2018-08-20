
#### Inputs required ####
# ensemble (n runs, l-dimensional output field), lxn matrix (each run of the model is a column)
# observations (vector of length l)


#### Formatting data ####
# MakeDataBasis formats the ensemble so that later functions can be applied
## $CentredField - the ensemble that the basis will be calculated from. 
##                 If RemoveMean = FALSE, then this is just the ensemble
##                 If RemoveMean = TRUE, the ensemble mean is calculated, and removed
## $EnsembleMean - the ensemble mean, a vector of length l
## $tBasis - the SVD (or weighted SVD basis), calculated across CentredField
## $Q, $Lambda - from the eigendecomposition of W^{-1}. These are used if calculating the
##               weighted SVD basis, and will be useful later when calculating the residual
##               basis at each iteration of the rotation algorithm



#### Setting W, W^{-1} ####
# The code defaults to W = identity matrix (i.e. uncorrelated, equal errors)
# We think of W as the sum of the observation error and discrepancy matrices, as this
# gives a parallel between the reconstructon error and history matching, and we can use
# the history match bound to assess whether the basis representation of the observations
# will be ruled out.
# If W is not the identity matrix, calculate the inverse via GetInverse, as this assigns 
# attributes to W^{-1}:
## $diagonal - is W^{-1} a diagonal matrix?
## $identity - is W^{-1} the identity matrix?
# These allow computational improvements to be made in other functions


