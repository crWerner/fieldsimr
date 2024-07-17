# ======================================================================
#
# Script: Simulate phenotypic line breeding program with moderate GEI
#
# Script authors: J. Bancic and  D.J. Tolhurst
#
# ======================================================================
#
# Description:
# This script demonstates how to simulate a phenotypic line breeding program
# with GEI using FieldSimR functionality. AlphaSimR is used to generate an
# additive genetic trait with population structure and FieldSimR is used to
# generate correlated plot errors based on a randomised complete block design.
# Each year, the simulation tracks the breeding progress by summarizing measures
# of genetic mean (average of genotype main effects in MET and TPE), variance
# (variance of genotype main effects in MET and TPE), accuracy of main genotype
# effects in MET and TPE, and MET-TPE alignment.
#
# The simulation involves six steps:
# 1. Import global parameters and a between-environment genetic variance
#     matrix, Ge, that constitutes the TPE.
# 2. Simulate founder parents with AlphaSimR and their genotype slopes with
#     FieldSimR's wrapper function multi_asr_input().
# 3. Sample a subset of environments for each year's MET from the full TPE.
# 4. Fill the breeding pipeline with unique genotypes using AlphaSimR
# 5. Perform yearly selections across breeding stages based on phenotypic.
#     values simulated with FieldSimR's wrapper functions multi_asr_output()
#     for simulating genetic values, field_trial_error() for simulating plot
#     errors with spatial variation and make_phenotypes() to combine genotype
#     values and errors.
# 6. Plot measures of genetic mean and variance (both observed in MET and
#     expected in TPE), accuracy of phenotypic selection (both observed
#     in MET and expected in TPE), and MET-TPE alignment.
# ======================================================================

# Load necessary libraries
rm(list = ls())  # Clean up the environment
library(AlphaSimR)
library(FieldSimR)
library(data.table)

# Scenario name
Scenario <- "Pheno"
GEI <- "Moderate"

# Import simulation parameters
source("GlobalParameters.R")

# Import simulated target population of environments (TPE)
Ge <- readRDS("Ge_ModerateGEI.rds") # Load presimulated Ge for moderate GEI from the github
Ce <- cov2cor(Ge)
De <- diag(diag(Ge))

# Run 20 years of breeding across multiple replicates
for (Rep in 1:nReps) {
  cat("\nWorking on replicate:", Rep, "\n")
  T
  # Simulate founders
  source("CreateFounders.R")

  # Import dataframe for storing variables
  source("StoreVariables.R")
  output$rep <- Rep

  # Sample MET envrionments from TPE for each simulation year
  # Note: Environments within each sample are sampled at random
  # and without ordering
  covs_tpe    <- input_asr$cov.mat
  samples_met <- sample_met(nenvs = nEnvTPE,
                            nsamples = nCycles,
                            sample.size = nEnvEYT,
                            replace = TRUE,
                            cov.mat = covs_tpe)$sample

  # Fill breeding pipeline
  source("FillPipeline.R")

  # Run simulation
  for (year in 1:nCycles) {
    time <- timestamp(prefix = "", suffix = "", quiet = TRUE)
    cat("Working on year: ", year, "    (", time, ")\n", sep = "")

    # Select new parents
    source("UpdateParents.R")

    # Advance the breeding year
    source("AdvanceYear.R")
  }

  # Save results for this replicate
  cat("Saving results \n")
  file_name <- paste0("Output_", Scenario, GEI, ".csv")
  write.table(output, file_name, sep = ",",
    col.names = !file.exists(file_name), row.names = FALSE, append = TRUE)
}

# Plot results
source("PlotResults.R")


# end of script

