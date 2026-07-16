
#######################################
#
# Supplementary Script S1
#
# Working example: Application to the two
# hypothetical environments using Base R
# functions
#
# Script Author: D.J. Tolhurst
#
#######################################

# This example generates genotype effects for 20 genotypes
# in 2 environments, with genetic variances of 1 and 3, and
# genetic correlation of 0. It will be demonstrated how to
# disentangle non-crossover and crossover genotype by
# environment interaction (GEI).

# The example involves an initial simulation step, followed by
# four steps from the new framework:
#
#  0. Generate variance matrix between environments
#     and genotype effects for each environment
#  1. Obtain genotype main effects and genotype by environment
#     interaction effects
#  2. Calculate the regression coefficient between the main effects
#     and interaction effects for each environment
#  3. Construct non-crossover and crossover effects by transferring
#     all variation in the interaction effects arising from perfect
#     positive correlation with the main effects
#  4. Investigate non-crossover and crossover interaction
#     from the perspective of the genotypes and environments
#
# It will then be demonstrated how to investigate GEI from the
# perspective of the environments and genotypes, including
# measures of expected genetic gain and overall genotype
# performance and stability.


require(devtools)
require(FieldSimR)
require(ggplot2)


# Some functions
pop_var <- function(x, y = NULL) {
  n <- if (is.null(y)) NROW(x) else length(x)
  var(x, y) * (n - 1) / n
}

plot_reaction <- function(intercepts, covariates, slopes) {

  x <- drop(covariates)
  b <- drop(slopes)
  g <- drop(intercepts)

  rn <- g + outer(b, x)
  cols <- hcl.colors(length(g), "Spectral")

  matplot(x, t(rn),
          type = "o", pch = 16, lty = 1,
          col = cols,
          xlab = "Environmental covariate",
          ylab = "Genetic effect"); abline(v = 0, lty = 2)
}

plot_regression <- function(x, y, slope) {
  intercept <- mean(y) - slope * mean(x)
  plot(x, y, pch = 16,
       xlab = "x", ylab = "y")
  abline(a = intercept, b = slope, lwd = 2)
}

plot_table <- function(x) {
  heatmap(x, Rowv = NA, Colv = NA, scale = "none",
  col = colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlGn"))(100),
  xlab = "Environment", ylab = "Genotype")
}

plot_effects <- function(x) {
  matplot(t(x), col = hcl.colors(nrow(x), "Dark 3"),
          type = "o", pch = 16, lty = 1,
          xaxt = "n",
          xlab = "Environment", ylab = "Effect")
  axis(1, at = seq_len(ncol(x)), labels = colnames(x))
}


#######################################
# Step 0: Obtain variance matrix between environments and
# genotype effects within environments

# Number of genotypes and environments
nids <- 20
nenvs <- 2

# Construct variance matrix and simulate a genotype by environment table
(var_ge <- matrix(c(1, 0, 0, 3), ncol = nenvs))
# 1  0
# 0  3
set.seed(123) # reproduce results in this script
ge_mat <- scale(matrix(rnorm(nenvs*nids), ncol = nenvs))
ge_mat <- ge_mat %*% solve(chol(pop_var(ge_mat))) %*% chol(var_ge)
pop_var(ge_mat) # check variance matrix
# Note that we will use the population variance,
# since the effects are simulated (and hence known)

# Plot the genotype by environment table
plot_table(ge_mat)
plot_effects(ge_mat)


#######################################
# Step 1: Obtain genotype main effects and
# genotype by environment interaction effects

# Main effects
g_mat <- cbind(rowMeans(ge_mat))
# Main effect variance
pop_var(g_mat)
(var_g <-  mean(var_ge)) # equivalence using variance parameters
# 1
# Plot the main effects
plot_table(g_mat[, c(1,1)])
plot_effects(g_mat[, c(1,1)])

# Interaction effects
gxe_mat <- ge_mat - g_mat[, c(1,1)]
# Interaction variance matrix
pop_var(gxe_mat)
(var_gge <- rowMeans(var_ge) - var_g) # covariance between main effects and interaction effects
(var_gxe <- var_ge - matrix(c(var_g + 2*var_gge[1], var_g,
                              var_g, var_g + 2*var_gge[2]), ncol = 2)) # equivalence using variance parameters
#  1 -1
# -1  1
# Plot the interaction effects with plot_table()
plot_table(gxe_mat)
plot_effects(gxe_mat)

# Quick check:
range(ge_mat - (g_mat[, c(1,1)] + gxe_mat))
# 0


# ----------------------------------- #
# Note that we can construct (hypothetical) environmental
# covariates and genotype slopes for the interaction
# effects as:
cov_mat <- c(-1, 1)
slope_mat <- cbind(gxe_mat[,2])
# Quick check:
range(gxe_mat - (slope_mat %*% t(cov_mat)))
# 0
(var_slope <- pop_var(slope_mat))
# 1
cov_mat %*% var_slope %*% cov_mat # as above for var_gxe
#  1 -1
# -1  1
# Note that the intercepts in this regression are taken as the main effects.
# Obtain the variance matrix between the intercepts (main effects)
# and slopes as:
pop_var(cbind(g_mat, slope_mat))
(var_ms <- matrix(c(var_g, rowMeans(var_ge)[2] - var_g,
                    rowMeans(var_ge)[2] - var_g, var_gxe[2,2]), ncol = 2)) # equivalence using variance parameters
#   1  0.5
# 0.5    1
# we will use these covariates and slopes in parallel to the
# main effects and interaction effects above, to demonstrate
# the important analogy.
# Recreate the variance matrix between environments
cbind(1, cov_mat) %*% var_ms %*% t(cbind(1, cov_mat)) # as above for var_ge
# 1  0
# 0  3
# Plot the reaction norms
plot_reaction(intercepts = g_mat,
              slopes = slope_mat,
              covariates = cov_mat)


#######################################
# Step 2: Calculate the coefficients (slopes) from regressing
# the interaction effects for each environment onto the
# main effects

# Calculate the regression slopes
c(solve(pop_var(g_mat)) %*% pop_var(g_mat, gxe_mat))
(bj <- var_gge/var_g) # equivalence using variance parameters
# -0.5 0.5

# Plot the regressions on the main effects with plot_regression()
plot_regression(x = g_mat,
                y = gxe_mat[,1],
                slope = bj[1])
plot_regression(x = g_mat,
                y = gxe_mat[,2],
                slope = bj[2])
# with slopes given by bj

# Note that we can also obtain the regression of the genotype
# effects for each environment onto the main effects:
plot_regression(x = g_mat,
                y = ge_mat[,1],
                slope = 1+bj[1])
plot_regression(x = g_mat,
                y = ge_mat[,2],
                slope = 1+bj[2])
# with slopes given by
(1 + bj)
c(cov(g_mat, ge_mat))/c(var(g_mat)) # equivalence using genetic effects
rowMeans(var_ge)/var_g # equivalence using variance parameters
# 0.5 1.5


# ----------------------------------- #
# Similar to the above, we can calculate the coefficient (slope)
# from regressing the genotype slopes onto the main effects (intercepts)
# Calculate the regression slope
pop_var(g_mat, slope_mat)/pop_var(g_mat)
(br <- var_ms[1,2]/var_g) # equivalence using variance parameters
# Noting we can retrieve bj as:
(cov_mat * br)
# -0.5 0.5
# Plot the regression on the main effects with plot_regression()
plot_regression(x = g_mat,
                y = slope_mat,
                slope = br)
# with slopes given by br


#######################################
# Step 3: Construct non-crossover and crossover effects
# by transferring all variation in the interaction effects
# arising from perfect positive correlation with the main
# effects

# Construct non-crossover effects
ngxe_mat <- g_mat[, c(1,1)] + g_mat %*% bj
cor(ngxe_mat) # check for perfect positive correlation between environments
# Plot the non-crossover effects with plot_table()
plot_table(ngxe_mat)
plot_effects(ngxe_mat)

# Construct crossover effects
cgxe_mat <- gxe_mat - g_mat %*% bj
# 0
cor(cgxe_mat) # check for perfect negative correlation between environments
cor(cgxe_mat, g_mat) # check for independence between non-crossover and crossover effects
# Plot the crossover effects with plot_table()
plot_table(cgxe_mat)
plot_effects(cgxe_mat)


# ----------------------------------- #
# Note that we can also obtain the non-crossover and
# crossover effects from the regression on (hypothetical)
# environmental covariates by transferring all variation
# in the slopes arising from perfect positive correlation
# with the main effects (intercepts)

# Construct the non-crossover effects
ngxe_mat
g_mat_adj <- g_mat[, c(1,1)] + g_mat %*% (cov_mat * br)
# This is exactly equal to the non-crossover effects, but note the
# regression here with intercepts and slopes both given by the main effects
# and covariates given by bj = cov_mat * br
plot_reaction(covariates = bj, # note the covariate is centred at 0
              intercepts = g_mat,
              slopes = g_mat)

# Construct crossover effects
slope_mat_adj <- slope_mat - g_mat %*% br
plot_reaction(covariates = cov_mat,
              slopes = slope_mat_adj,
              intercepts = rep(0, nids))

# Recall, the reaction norm prior to disentangling the non-crossover and crossover effects:
plot_reaction(covariates = cov_mat,
              slopes = slope_mat,
              intercepts = g_mat)
# Therefore, the reactions have been completely disentangled
# into non-crossover and crossover interaction


#######################################
# Step 4: Investigate non-crossover and crossover interaction
# from the perspective of the genotypes and environments

# a) Non-crossover and crossover variances for each environment
# b) Genotype main effects and stability
# c) Response to selection and expected genetic gain
# d) Altering the selection focus


#######################################
# a) Non-crossover and crossover variances for each environment

# Calculate the non-crossover variances
diag(pop_var(ngxe_mat)) # for the two environments
(var_n <- (1+bj)^2*var_g) # equivalence using variance parameters
mean(var_n) # overall
# 1.25

# Calculate the crossover variances
diag(var_c <- pop_var(cgxe_mat))
diag(var_gxe) - bj^2*var_g # equivalence using variance parameters
mean(diag(var_c)) # overall
# 0.75

# Compare to the traditional main effect and interaction variances
# Print the main effect variance
var_g
# 1
# Print the interaction variance
diag(var_gxe) # equivalence using variance parameters
# 1 # overall


# ----------------------------------- #
# In addition to the above, we can inspect the
# non-crossover and crossover variances along the
# (hypothetical) environmental covariate

# Construct a hypothetical covariate spanning -2 to 2
cov_mat_all <- seq(-2,2,0.1)

# Calculate the non-crossover variance (recall bj = cov_mat * br)
diag(pop_var(ngxe_mat)) # for the two environments
(1 + cov_mat * br)^2*var_g # equivalence using variance parameters
# 0.25 2.25
plot(cov_mat_all, (1 + cov_mat_all * br)^2*var_g, type = "l") # for all values along the covariate
# Calculate the crossover variance
diag(pop_var(cgxe_mat)) # for the two environments
cov_mat^2*c(var_slope) - (cov_mat * br)^2*var_g # equivalence using variance parameters
# 0.75 0.75
plot(cov_mat_all, cov_mat_all^2*c(var_slope) - (cov_mat_all * br)^2*var_g, type = "l") # for all values along the covariate
# Calculate the total variance
diag(pop_var(ge_mat)) # for the two environments
(1 + cov_mat * br)^2*var_g + cov_mat^2*c(var_slope) - (cov_mat * br)^2*var_g # equivalence using variance parameters
plot(cov_mat_all, (1 + cov_mat_all * br)^2*var_g +
                   cov_mat_all^2 * c(var_slope) - (cov_mat_all * br)^2*var_g, type = "l") # for all values along the covariate
# Note that the covariance between these components is 0 for all values
# of the covariate.
tot_var <- c(var_g) + cov_mat_all^2*c(var_slope) + 2*cov_mat_all*c(pop_var(g_mat, slope_mat))
plot(cov_mat_all, tot_var, type = "l", ylim = c(-2, 7))
ncross_var <- (1+ cov_mat_all*br)^2 * var_g
lines(cov_mat_all, ncross_var, col = "darkgreen")
cross_var <- tot_var - ncross_var
lines(cov_mat_all, cross_var, col = "darkblue")

# Compare to the traditional main effect and interaction variances
pop_var(g_mat) # Calculate the main effect variance
var_g # equivalence using variance parameters
# 1
plot(cov_mat_all, rep(1, 41)*var_g, type = "l") # for all values along the covariate
# Calculate the interaction variance
diag(pop_var(gxe_mat)) # for the two environments
cov_mat^2*c(var_slope) # equivalence using variance parameters
# 1 1
plot(cov_mat_all, cov_mat_all^2 * c(var_slope), type = "l") # for all values along the covariate
# Note that the covariance between these components is not 0, except for the point
# where the intercepts and slopes are independent.
2 * cov_mat * var_ms[1,2] # for the two environments
plot(cov_mat_all, 2*cov_mat_all*var_ms[1,2], type = "l") # for all values along the covariate
# Calculate the total variance
c(1, 1)* var_g + cov_mat^2*c(var_slope) + 2*cov_mat %*% var_ms  # for the two environments
plot(cov_mat_all, rep(1, 41)*var_g +
                    cov_mat_all^2 * c(var_slope) +
                    2*cov_mat_all*var_ms[1,2], type = "l") # for all values along the covariate
# Note the same total variance as above, but notably different partitioning
plot(cov_mat_all, tot_var, type = "l", ylim = c(-2, 7))
lines(cov_mat_all, rep(var_g, length(cov_mat_all)), col = "darkgreen")
int_var <- cov_mat_all^2*c(var_slope)
lines(cov_mat_all, int_var, col = "darkblue")
cov_var <- 2*cov_mat_all*c(pop_var(g_mat, slope_mat))
lines(cov_mat_all, cov_var, col = "purple")


#######################################
# b) Genotype main effects and stability

# Obtain genotype main effects
g_mat
rowMeans(ngxe_mat) # equivalence using mean of non-crossover effects for each genotype
# Calculate genotype stability
stability_mat <- rowMeans(cgxe_mat^2) # variance of crossover effects for each genotype
plot(stability_mat, g_mat)
# Check mean of stability
mean(stability_mat)
mean(diag(var_c)) # equivalent to the crossover variance, by design
# 0.75

# Compare to traditional main effects and variances
stability_old_mat <- rowMeans(gxe_mat^2)
plot(stability_old_mat, g_mat)
plot(g_mat, g_mat); abline(a=0, b=1) # same main effects
plot(stability_mat, stability_old_mat); abline(a=0, b=1) # different stability
# Note the traditional stability is inflated compared to the new measures
# This is emphasiseed by the mean of traditional stability
mean(stability_old_mat)
mean(diag(var_gxe)) # equivalent to the interaction variance
# 1


#######################################
# c) Response to selection and adaptation potential

# Calculate the expected genetic gain at the selection focus
pop_var(g_mat, ge_mat)/c(pop_var(g_mat)) # for the two environments
(var_gge + var_g)/var_g # equivalence using variance parameters, i.e (1+bj)*sqrt(var_g)
sqrt(var_n) # equivalent to the non-crossover standard deviations
# 0.5  1.5
mean((var_gge + var_g)/var_g) # overall
sqrt(var_g) # equivalent to the main effect standard deviation
mean(sqrt(var_n)) # equivalent to the mean of the non-crossover standard deviations
# 1

# Calculate the expected genetic variability around the selection focus
# pop_var(g_mat, ge_mat)/c(pop_var(g_mat)) # contribution from the two environments
sqrt(diag(var_gxe) - bj^2*var_g) # equivalence using variance parameters, i.e
# (diag(var_gxe) - bj*var_g)^2

sqrt(diag(var_c)) # equivalent to the crossover standard deviations
# 0.866 0.866
# mean((var_gge + var_g)/var_g) # overall
sqrt(mean(diag(var_c))) # equivalent to the overall crossover standard deviation
# 0.866
mean(diag(var_gxe)) - mean(bj^2*var_g)


#######################################
# end of script

save.image()

