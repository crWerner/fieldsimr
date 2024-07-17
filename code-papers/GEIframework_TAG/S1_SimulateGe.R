
##########################################################
#
# Supplementary Script S1: A framework for simulating GEI
#
# Simulating a between-environment genetic variance matrix
#
# Script authors: D.J. Tolhurst and J. Bancic
#
##########################################################

# This script demonstrates how to simulate a reduced rank between-environment
# genetic variance matrix, Ge, that captures GEI. The genetic variance matrix
# will be summarised using measures of variance explained, and will form the
# basis of a reduced rank multiplicative model in the next script.
# The simulation is demonstrated with FieldSimR functions.

# The simulation involves four steps:
# 1. Simulating a diagonal genetic variance matrix using rand_diag_mat()
# 2. Simulating a between-environment genetic correlation matrix using struc_cor_mat()
# 3. Constructing the between-environment genetic variance matrix by combining the
#    matrices in Steps 1. and 2.
# 4. Obtaining measures of variance explained using measure_interactions()

# It is also demonstrated how to simulate the between-environment genetic covariance matrix
# with specific structure, including genotype by year and genotype location interaction,
# as well as multiple TPE and multiple phenotypic traits

# Install FieldSimR from github using devtools
# devtools::install_github("crWerner/fieldsimr")

# Load FieldSimR
library(FieldSimR)

# Simulation parameters
p <- 1000      # number of environments in the TPE
# Genetic variance parameters
alpha <- 1.5   # shape of the gamma distribution
beta  <- 1     # scale of the gamma distribution
# Genetic correlation parameters
rho <- 0.2     # baseline genetic correlation between environments in the TPE
k <- 7         # rank of the between-environment genetic correlation matrix
gamma <- -0.5  # skewness of the genetic correlation distribution


# The simulation
# 1. Simulate a diagonal genetic variance matrix, De, for all environments in the TPE
# using FieldSimR's function rand_diag_mat().
De <- skew_diag_mat(n = p,
                    shape = alpha,
                    scale = beta)
# Plot the genetic variances using FieldSimR's function plot_hist()
plot_hist(df = diag(De))
mean(diag(De)) # mean genetic variance within environments


# 2. Simulate a between-environment genetic correlation matrix, Ce, for all environments
# in the TPE using FieldSimR's function struc_cor_mat().
Ce <- struc_cor_mat(n = p,
                    base.cor = rho,
                    rank = k,
                    skew = gamma)
# Plot the correlation matrix using FieldSimR's functions plot_hist() and plot_matrix()
plot_hist(df = Ce[upper.tri(Ce)])
mean(Ce[upper.tri(Ce)])
plot_matrix(Ce,
            order = TRUE,
            labels = FALSE)


# 3. Construct the between-environment genetic variance matrix, Ge.
Ge <- sqrt(De) %*% Ce %*% sqrt(De)
# Plot the variance matrix using FieldSimR's functions plot_hist() and plot_matrix()
plot_hist(df = Ge[upper.tri(Ge)])
mean(Ge[upper.tri(Ge)]) # mean genetic correlation between environments


# 4. Obtain measures of variance explained to tune the matrices in Steps 1 and 2 using
# FieldSimR's function measure_interactions().
measure_variances(mat = Ge)
# The mean genetic variance is partitioned into main effect and interaction variance.
# The interaction variance is further partitioned into heterogeneity of genetic variance
# and lack of correlation.
# Also presented are new measures which partition the mean genetic variance into
# non-crossover and crossover variance, with the former representing the variation
# in each environment attributed to perfect positive correlation with the genotype
# main effects.


# The below is not strictly required to run Supplementary Script S2.
# It is also possible to simulate Ce with various structure, including
# (a) genotype by year with decaying correlations over time
# (b) genotype by year and genotype by location interaction
# (c) multiple groups representing multiple TPE or traits

# (a) Simulate Ce to capture genotype by year interaction with decaying correlations over
# time based on an autoregressive process of order one for the base function
pm <- 20       # number of environments
rho_ar1 <- 0.7 # baseline genetic correlation
ar1_base <- rho_ar1^abs(outer(X = 1:pm, Y = 1:pm, FUN = "-")) # base function
# Plot the base correlation matrix using FieldSimR's function plot_matrix()
plot_matrix(ar1_base)
# Add structured noise to the base correlation matrix using FieldSimR's function struc_cor_mat()
epsilon <- (1-rho_ar1)/(1+rho_ar1) # magnitude of noise, ensures Ce is positive (semi)-definite
Ce_a <- struc_cor_mat(base.mat = ar1_base,
                      range = epsilon,
                      rank = 3,
                      skew = 0)
# Plot Ce using FieldSimR's functions plot_hist and plot_matrix()
plot_hist(df = Ce_a[upper.tri(Ce_a)])
mean(Ce_a[upper.tri(Ce_a)])
plot_matrix(Ce_a)
# The above is analogous to simulating Ce to capture genotype by location interaction with
# decaying correlations over space

# (b) Simulate Ce to capture genotype by year and genotype by location interaction
# based on a variance component model for the base function
ym <- 4         # number of years
lm <- 5         # number of locations
pm <- ym*lm     # number of environments
prop_g <- 0.2   # proportion of genotype main effect variance
prop_gxy <- 0.3 # proportion of genotype by year variance
prop_gxl <- 0.2 # proportion of genotype by location variance
prop_gxe <- 0.3 # proportion of genotype by location variance
base_g <- matrix(1, ncol = pm, nrow = pm) * prop_g  # Base function for main effects
base_gxy <- kronecker(diag(ym), matrix(1, ncol = lm, nrow = lm)) * prop_gxy # base function for genotype by year interaction
base_gxl <- kronecker(matrix(1, ncol = ym, nrow = ym), diag(lm)) * prop_gxl # base function for genotype by location interaction
base_gxe <- diag(pm) * prop_gxe                     # base function for genotype by location interaction
base_var <- base_g + base_gxy + base_gxl + base_gxe # overall base function
# Plot the base correlation matrix using FieldSimR's function plot_matrix()
plot_matrix(base_var)
# Add structured noise to the base correlation matrix using FieldSimR's function struc_cor_mat()
epsilon <- 0.2 # magnitude of noise
Ce_b <- struc_cor_mat(base.mat = base_var,
                      range = epsilon,
                      rank = 3,
                      skew = 0)
# Plot Ce using FieldSimR's functions plot_hist and plot_matrix()
plot_hist(df = Ce_b[upper.tri(Ce_b)])
mean(Ce_b[upper.tri(Ce_b)])
plot_matrix(Ce_b,
            group.df = rep(1:ym, each = lm),
            order = TRUE)

# (c) Simulate Ce with multiple groups representing multiple TPE or traits based on a block diagonal matrix
p1 <- 5   # number of environments in group 1
p2 <- 15  # number of environments in group 2
pm <- p1 + p2 # total number of environments
rho <- 0   # baseline genetic correlation between groups
rho1 <- 0.3  # baseline genetic correlation within group 1
rho2 <- 0.5  # baseline genetic correlation within group 2
base <- matrix(1, ncol = pm, nrow = pm) * rho # base function between groups
base1 <- matrix(1, ncol = p1, nrow = p1) * (rho1-rho) # base function within group 1
base2 <- matrix(1, ncol = p2, nrow = p2) * (rho2-rho) # base function within group 2
base_group <- base + magic::adiag(base1, base2) # overall base function
diag(base_group) <- 1
# Plot the base correlation matrix using FieldSimR's function plot_matrix()
plot_matrix(base_group)
# Add structured noise to the base correlation matrix using FieldSimR's function struc_cor_mat()
epsilon <- 0.2 # magnitude of noise
Ce_c <- struc_cor_mat(base.mat = base_group,
                      range = epsilon,
                      rank = 3,
                      skew = 0)
# Plot Ce using FieldSimR's functions plot_hist and plot_matrix()
plot_hist(df = Ce_c[upper.tri(Ce_c)])
mean(Ce_c[upper.tri(Ce_c)])
plot_matrix(Ce_c,
            group.df = rep(1:2, times = c(p1, p2)),
            order = TRUE)
# The above can be tailored to simulating Ce for multiple TPE for multiple phenotypic traits
# by considering no.tpe x no.trait groups.



# end of script



