##########################################################
#
# Simulating a between-environment genetic variance matrix
#
# Script authors: D.J. Tolhurst and J. Bancic
#
##########################################################

# This script demonstrates how to simulate a reduced rank between-environment 
# genetic variance matrix, Ge, that captures GEI. The genetic variance matrix 
# will be summarised using measures of variance explained, and will form the
# basis of a reduced rank multiplicative model in the next script
# The simulation is demonstrated with FieldSimR functions.

# The simulation involves four steps:
# 1. Simulating a diagonal genetic variance matrix using rand_diag_mat()
# 2. Simulating a between-environment genetic correlation matrix using rr_cor_mat()
# 3. Constructing the between-environment genetic variance matrix by combining the 
#    matrices in Steps 1. and 2.
# 4. Obtaining measures of variance explained using measure_interactions()

# Load FieldSimR
library(FieldSimR)

# Simulation parameters
p <- 1000      # Number of environments in the TPE
# Genetic variance parameters
alpha <- 1.5  # shape of the inverse gamma distribution
beta <- 1     # scale of the inverse gamma distribution 
# Genetic correlation parameters
rho <- 0       # baseline genetic correlation between environments in the TPE
k <- 7         # rank of the between-environment genetic correlation matrix
gamma <- -0.7  # skewness of the genetic correlation distribution


# The simulation
# 1. Simulate a diagonal genetic variance matrix for all environments in the TPE 
# using FieldSimR's function rand_diag_mat().
De <- skew_diag_mat(n = p,
                    shape = alpha,
                    scale = beta)
# Plot the genetic variances using FieldSimR's function plot_hist()
plot_hist(df = diag(De))
mean(diag(De)) # mean genetic variance within environments


# 2. Simulate a between-environment genetic correlation matrix for all environments 
# in the TPE using FieldSimR's function rr_cor_mat().
Ce <- rr_cor_mat(n = p,
                 base.cor = rho,
                 rank = k,
                 skew = gamma)
# Plot the correlation matrix using FieldSimR's functions plot_hist() and plot_matrix()
plot_hist(df = Ce[upper.tri(Ce)])
mean(Ce[upper.tri(Ce)])
plot_matrix(Ce[1:200,1:200],
            order = TRUE,
            labels = FALSE)


# 3. Construct the between-environment genetic variance matrix.
Ge <- sqrt(De) %*% Ce %*% sqrt(De)
# Plot the variance matrix using FieldSimR's functions plot_hist() and plot_matrix()
plot_hist(df = Ge[upper.tri(Ge)])
mean(Ge[upper.tri(Ge)])
plot_matrix(Ge[1:200,1:200],
            order = TRUE,
            labels = FALSE)


# 4. Obtain measures of variance explained to tune the matrices in Steps 1 and 2 using 
# FieldSimR's function measure_interactions().
measure_interactions(mat = Ge)


# 5. Potential to simulate some Ge with groupings
# - ar1 base function
# - gxy and gxl using an extended version of rr_cor_mat

# Generate MET samples
# n_samples = 20
# p_year = 20
# sample_covs = lapply(X = 1:n_samples, 
#                      function(X) sample(1:p, size = p_year, replace = F))


## Save RData
rm(alpha,beta,gamma,rho,n_samples,p_year,p)
save.image("simulatedGe.RData")

# end of script



