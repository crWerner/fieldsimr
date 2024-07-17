
##########################################################
#
# Supplementary Script S2: A framework for simulating GEI
#
# Simulating a multi-environment trial dataset
#
# Script authors: D.J. Tolhurst and J. Bancic
#
##########################################################

# This script demonstrates how to simulate a MET dataset that captures GEI and
# spatial variation within environments. The environments in the MET dataset
# will be sampled from a TPE, and will be summarised using
# measures of heritability, expected accuracy and MET-TPE alignment. The
# methods are demonstrated with FieldSimR and AlphaSimR, which are used
# to simulate GEI based on a reduced rank multiplicative model.

# The simulation involves six steps:
# 2. Simulating genetic values in AlphaSimR using FieldSimR's wrapper functions
#    multi_asr_input() and multi_asr_output()
#    - Simulating GE effects with GEI for all environments that constitute the TPE
# 3. Sampling a subset of environments for a MET from the full TPE and summarising the
#    MET-TPE alignment
# 4. Simulating plot errors with spatial variation using field_trial_error()
# 5. Constructing a MET dataset by combining the effects in 3. and 4. using make_phenotypes()
# 6. Obtaining measures of expected accuracy to summarise MET dataset

# Install FieldSimR from github using devtools
# devtools::install_github("crWerner/fieldsimr")

# Load FieldSimR and AlphaSimR
library(FieldSimR)
library(AlphaSimR)

# Simulation parameters
# Ge           # between-environment genetic variance matrix, obtained from Supplementary Script S1
load("Ge_matrices.RData") # alternatively, load a presimulated Ge with moderate GEI from the github
Ge <- Ge_mod
p  <- 1000    # number of environments in the TPE
v  <- 400     # number of genotypes in the breeding population
r  <- 2       # number of replicate blocks in each environment
# i.e. v*p = 400,000 genotype by environment combinations
# Genetic architecture and trait parameters
nchr <- 20   # number of chromosomes
nQTN <- 200  # number of QTN per chromosome
mu <- 4      # overall trait mean
H2 <- 0.3    # overall plot-level heritability
# Model parameters
k <- 7       # number of multiplicative terms (equivalent to rank of Ge)


# 1. Simulate genetic values in AlphaSimR using FieldSimR's wrapper functions
# multi_asr_input() and multi_asr_output()
# The function multi_asr_input() decomposes the between-environment genetic variance matrix
# to obtain the environmental covariates and stores the variance matrix of the genotype slopes
# to simulate use with AlphaSimR to generate appropriate population structure.
library(FieldSimR)
Ce <- cov2cor(Ge) # obtain between-environment genetic correlation matrix
de <- diag(Ge)    # obtain genetic variances
input_asr <- multi_asr_input(nenvs = p,
                             mean = mu,
                             var   = de,
                             corA  = Ce,
                             nterms = k)
covs_tpe <- input_asr$cov.mat

# Simulate a population of inbred genotypes
library(AlphaSimR)
pop <- quickHaplo(nInd = v,
                  nChr = nchr,
                  segSites = nQTN,
                  inbred = TRUE)
SP <- SimParam$new(pop)
SP$addTraitA(nQtlPerChr = nQTN,
             mean = input_asr$mean,
             var  = input_asr$var,
             corA = input_asr$corA)
pop <- newPop(pop)


# Obtain the true genetic values for all environments in the TPE using FieldSimR's
# wrapper function multi_asr_output()
gvs_tpe <- multi_asr_output(pop = pop,
                            nenvs = p,
                            nreps = r,
                            cov.mat = covs_tpe,
                            return.effects = TRUE)
head(gvs_tpe$gv.df) # genetic values for all environments in the TPE
head(gvs_tpe$Terms) # genotype slopes for all k terms
gvs_tpe <- gvs_tpe$gv.df # just take genetic values

# 3. Sample environments for the MET from the TPE using FieldSimR's function sample_met()
pm <- 20   # number of environment to be sampled in the MET dataset
v <- 400   # number of genotypes to be sampled in the MET dataset
r <- 2     # number of replicate blocks in each environment

envs_met <- sample_met(nenvs = p,
                       nsamples = 1000,
                       sample.size = pm,
                       replace = TRUE,
                       cov.mat = covs_tpe)
head(envs_met$sample)  # environments sampled in each MET
head(envs_met$cov.mat) # environmental covariates sampled in each MET


# Obtain true genetic values for the sample of genotypes and environments in the MET
gvs_met <- droplevels(gvs_tpe[gvs_tpe$env %in% envs_met$sample[[1]],]) # take first sample as an example
head(gvs_met) # genetic values for all genotypes and environments sampled in the MET

# Summarise the MET-TPE alignment based on the true mean genetic values in the MET and TPE
cor(with(gvs_met, tapply(gv.Trait1, id, mean)),
    with(gvs_tpe, tapply(gv.Trait1, id, mean)))
# note the expected MET-TPE alignment is
Ge_vars <- measure_variances(Ge)
sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm))


# 4. Simulate plot errors with spatial variation using FieldSimR's core function field_trial_error()
# with default parameters
# Generate errors variances to produce an overall heritability of 0.3
sigm2gj <- de[envs_met$sample[[1]]] # genetic variances in the MET
sigm2g <- mean(sigm2gj) # mean genetic variance in the MET
sigm2e <- sigm2g/H2 - sigm2g # mean error variance in the MET
sigm2ej <- diag(skew_diag_mat(n = pm, mean.var = sigm2e)) # simulate error variances in the MET
H2j <- c(sigm2gj/(sigm2gj+sigm2ej)) # heritabilities in the MET
# Display the heritabilities using FieldSimR's function plot_hist()
plot_hist(df = H2j)
mean(H2j)
# Simulate plot errors
error_df <- field_trial_error(nenvs = pm,
                              nblocks = r,
                              varR = sigm2ej,
                              spatial.model = "AR1",
                              ncols = 20,
                              nrows = 40,
                              return.effects = TRUE)
# Display the plot errors in Environment 1 using FieldSimR's function plot_effects()
plot_effects(error_df$error.df[error_df$error.df$env == 1,], effect = "e.Trait1")
plot_effects(error_df$Trait1[error_df$Trait1$env == 1,], effect = "e.spat")
plot_effects(error_df$Trait1[error_df$Trait1$env == 1,], effect = "e.rand")

# 5. Construct the MET dataset by combining the simulated genetic values from Step 3. with the
# plot errors from Step 4. using make_phenotypes(). Genotypes are allocated to plots according
# to a randomised complete block design.
met_df <- make_phenotypes(gv.df = gvs_met,
                          error.df = error_df,
                          randomise = TRUE,
                          return.effects = TRUE)
head(met_df$pheno.df) # phenotypes for all environments sampled in the MET
head(met_df$Trait1) # genetic values and errors for all environments sampled in the MET
met_df <- met_df$pheno.df # just take MET dataset
# Display the phenotypes in the first environment using FieldSimR's function plot_effects()
plot_effects(met_df[met_df$env == envs_met$sample[[1]][1],], effect = "y.Trait1")


# 6. Obtain measures of expected accuracy
# Main effect accuracy in the TPE (square-root of line-mean heritability across environments)
cor(with(met_df, tapply(y.Trait1, id, mean)), with(gvs_tpe, tapply(gv.Trait1, id, mean)))
# note the expected main effect accuracy in the TPE is
sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))

# Main effect accuracy in the MET
cor(with(met_df, tapply(y.Trait1, id, mean)), with(gvs_met, tapply(gv.Trait1, id, mean)))
# note the expected main effect accuracy in the MET is
sqrt((Ge_vars[1,2] + Ge_vars[2,2]/pm)/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))

# Accuracy of the GE effects in the MET (square-root of line-mean heritability within environments)
summary(diag(cor(with(met_df, tapply(y.Trait1, list(id, env), mean)),
         with(gvs_met, tapply(gv.Trait1, list(id, env), mean)))))
# note the expected accuracy of the GE effects in the MET is
sqrt((Ge_vars[1,2] + Ge_vars[2,2])/(Ge_vars[1,2] + Ge_vars[2,2] + sigm2e/r))


# end of script

