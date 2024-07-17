
##########################################################
#
# Supplementary Script S0: A framework for simulating GEI
#
# Simulating an example TPE and sampling a MET dataset
#
# Script authors: J. Bancic and D.J. Tolhurst
#
##########################################################

# This script simulates a small example TPE, from which a MET
# dataset is sampled. A randomised complete block design is
# used for each environment with two replicate blocks.

set.seed(123)
p  <- 100 # no. environments in TPE
pm <- 10  # no. environments in MET
v  <- 400 # no. genotypes per environment
r  <- 2   # no. replicate blocks per environment
mu <- 4   # trait mean
k  <- 7   # no. multiplicative terms

#> Simulate GE effects
# 1. Construct Ge matrix, representative of a TPE
Ce <- diag(1/2, p)
Ce[upper.tri(Ce)] <- runif(p*(p-1)/2); Ce <- Ce + t(Ce)
De <- diag(rgamma(p, shape = 1.5, scale = 1))
Ge <- De^(1/2) %*% Ce %*% De^(1/2)

# 2. Decompose Ge, taking the first k terms
eigen_decom <- eigen(Ge)
Uk <- eigen_decom$vectors[, 1:k]
Lk <- diag(eigen_decom$values[1:k])

# 3. Obtain environmental covariates and genotype slopes
Sk <- Uk
fk <- scale(matrix(rnorm(v*k), ncol = k)) %*% Lk^(1/2)

# 4. Construct GE effects in the TPE
u <- fk %*% t(Sk)

#> Construct MET dataset
# 1. Simulate environmental main effects
X <- kronecker(diag(pm), rep(1, v*r))
tau <- c(scale(rnorm(pm)))

# 2. Sample GE effects from TPE
um <- u[, sample(1:p, pm)]
# Check MET-TPE alignment
cor(rowMeans(um), rowMeans(u))
# Randomise genotypes to plots in blocks
env <- factor(rep(1:pm, each = v*r))
block <- factor(rep(1:r, each = v))
id <- unlist(lapply(seq_len(pm*r), function(x) factor(sample(1:v,v))))
Z <- model.matrix(~ id:env - 1)

# 3. Simulate plot errors
e <- c(scale(matrix(rnorm(v*pm*r, sd = 2), ncol = pm), scale = FALSE))

# 4. Generate phenotypes
y <- mu + X %*% tau + Z %*% c(um) + e
# Construct MET dataframe
MET.df <- data.frame(env = env,
                     block = block,
                     id  = id,
                     y = y,
                     u = Z %*% c(um))
# Check main effect accuracy in TPE and MET, and accuracy of GE effects in MET
cor(with(MET.df, tapply(y, id, mean)), rowMeans(u))
cor(with(MET.df, tapply(y, id, mean)), rowMeans(um))
summary(diag(cor(with(MET.df, tapply(y, list(id, env), mean)), um)))

# Check genetic variances
round(sigG <- apply(um,2,var),2)
mean(sigG)

# Check error variances
round(sigE <- apply(matrix(e, ncol = pm),2,var),2)
mean(sigE)

# Check heritabilities
round(sigG/(sigG + sigE),2)
mean(sigG/(sigG + sigE))

# end of script




