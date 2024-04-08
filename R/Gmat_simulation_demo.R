
##########################################################
#
# Simulating a between-environment genetic variance matrix
#
# Script authors: D.J. Tolhurst and J. Bancic
#
##########################################################

# This script demonstrates how to simulate a between-environment genetic variance
# matrix, Ge, and how to tune this matrix using measures of variance explained.
# The methods are demonstrated using FieldSimR's functions.

# The simulation involves three steps:
# 1. Simulating a diagonal genetic variance matrix using rand_diag_mat()
# 2. Simulating a between-environment genetic correlation matrix using rand_diag_mat()
# 3. Constructing the between-environment genetic variance matrix by combining the matrices in 1. and 2.
# 4. Obtaining measures of variance explained using measure_interactions()

# Load FieldSimR
library(FieldSimR)

# Simulation parameters
p <- 1000      # Number of environments in the TPE
min_var <- 0   # minimum genetic variance within environments
max_var <- 0.2 # maximum genetic variance within environments
rho <- 0       # baseline genetic correlation between environments in the TPE
k <- 7         # rank of the between-environment genetic correlation matrix
gamma <- -0.7  # skewness of the genetic correlation distribution


# The simulation
# 1. Use FieldSimR's function rand_diag_mat() to simulate a diagonal genetic
# variance matrix for all environments in the TPE
De <- rand_diag_mat(n = p,
                    min.var = min_var,
                    max.var = max_var)
De <- skew_diag_mat(n = p,
                    min.var = min_var,
                    max.var = max_var)
# Plot the genetic variances using FieldSimR's function plot_hist() <-- TODO
plot_hist(df = diag(De), density = TRUE)
mean(diag(De)) # mean genetic variance within environments


# 2. Use FieldSimR's function rr_cor_mat() to simulate a between-environment genetic
# correlation matrix for all environments in the TPE
Ce <- rr_cor_mat(n = p,
                 base.cor = rho,
                 rank = k,
                 skew = gamma)
# Plot the correlation matrix using FieldSimR's functions plot_hist() and plot_matrix()
hist(Ce[upper.tri(Ce)])
mean(Ce[upper.tri(Ce)])
plot_matrix(Ce[1:11,1:11],
            order = TRUE,
            labels = FALSE)


# 3. Construct the between-environment genetic variance matrix
Ge <- sqrt(De) %*% Ce %*% sqrt(De)


# 4. Obtain measures of variance explained using FieldSimR's function measure_interactions()
# to tune the matrices in Steps 1 and 2
measure_interactions(Ge)


# end of script


