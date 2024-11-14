
##########################################################
#
# Supplementary Script S3: A framework for simulating GEI
#
# Comparison of models fitted to a simulated MET dataset
#
# Script authors: J. Bancic and D.J. Tolhurst
#
##########################################################

# This script briefly demonstrates model comparison using ASReml-R

# It is important to note that the models here are intended
# for demonstration purposes only, and should be tailored to
# the research objectives.

# Load ASReml-R and FieldSimR
library(asreml)
library(FieldSimR)

# The following are generated in Supplementary Script S2
met_df   # MET data frame
gvs_met  # true genetic values in MET
gvs_tpe  # true genetic values in TPE
Ge       # true between-environment genetic variance matrix
load("MET_data.RData") # alternatively, load presimulated MET data from the github
pm <- nlevels(gvs_met$env)
v <- nlevels(gvs_met$id)
r <- nlevels(gvs_met$rep)

# Output data frame
modelNames <- c("Expected", "Main", "Comp", "MDiag", "Diag", "FA1")
output <- data.frame(
  GEI = rep("Moderate", length(modelNames)),
  Envs = rep(pm, length(modelNames)),
  Model = modelNames,
  LogLik = NA,
  AIC  = NA,
  r_g  = NA,
  r_m  = NA,
  r_mt  = NA,
  r_ge = NA,
  stringsAsFactors = FALSE
)
output

# Get true values
g_tpe_true <- with(gvs_tpe, tapply(gv.Trait1, id, mean))
g_met_true <- with(gvs_met, tapply(gv.Trait1, id, mean))
ge_met_true <- with(gvs_met, tapply(gv.Trait1, list(id, env), mean))
(Ge_vars <- measure_variances(Ge))
sigm2e

#-- Expected accuracies
# Main effect accuracy in TPE
output[output$Model == "Expected",]$r_g <-
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# Main effect accuracy in MET
output[output$Model == "Expected",]$r_m <-
  sqrt((Ge_vars[1,2] + Ge_vars[2,2]/pm)/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# MET-TPE alignment
output[output$Model != "Expected",]$r_mt <- cor(g_tpe_true, g_met_true) # observed
output[output$Model == "Expected",]$r_mt <-
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm))
# Accuracy for GE effects in MET
output[output$Model == "Expected",]$r_ge <-
  sqrt((Ge_vars[1,2] + Ge_vars[2,2])/(Ge_vars[1,2] + Ge_vars[2,2] + sigm2e/r))



#-- Model 1: Main effects only
asr_main <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df,
                   workspace = "1Gb")
asr_main <- update.asreml(asr_main)

# Predicted genotype main effects
g_met_pred <- asr_main$coefficients$random[
  grep(pattern = "^id", rownames(asr_main$coefficients$random)),]
# Repeated main effects across all environments to obtain genotype by environment effects
ge_met_pred <- matrix(g_met_pred, ncol = pm, nrow = v)
# Total genetic effects not calculated

# Store output
cat("LogLik: ", output[output$Model == "Main",]$LogLik <- asr_main$loglik,"\n")
cat("AIC: ",    output[output$Model == "Main",]$AIC <- summary.asreml(asr_main)$aic,"\n")
cat("r_g: ",  output[output$Model == "Main",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "Main",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "Main",]$r_ge <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")



#-- Model 2: Compound symmetry
asr_comp <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + env:id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df,
                   workspace = "1Gb")
asr_comp <- update.asreml(asr_comp)

# Predicted genotype main effects
g_met_pred <- asr_comp$coefficients$random[
  grep(pattern = "^id", rownames(asr_comp$coefficients$random)),]
# Predicted genotype by environment effects
ge_met_pred <- asr_comp$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr_comp$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge_met_pred

# Store output
cat("LogLik: ", output[output$Model == "Comp",]$LogLik <- asr_comp$loglik,"\n")
cat("AIC: ",    output[output$Model == "Comp",]$AIC <- summary.asreml(asr_comp)$aic,"\n")
cat("r_g: ",  output[output$Model == "Comp",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "Comp",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "Comp",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")



#-- Model 3: Main effects + diagonal
asr_mdiag <- asreml(y.Trait1 ~ 1 + env,
                    random   = ~ id + diag(env):id + diag(env):block,
                    residual = ~ dsum(~ ar1(col):ar1(row) | env),
                    data     = met_df,
                    workspace = "1Gb")
asr_mdiag <- update.asreml(asr_mdiag)

# Predicted genotype main effects
g_met_pred <- asr_mdiag$coefficients$random[
  grep(pattern = "^id", rownames(asr_mdiag$coefficients$random)),]
# Predicted genotype by environment effects
ge_met_pred <- asr_mdiag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr_mdiag$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge_met_pred

# Store output
cat("LogLik: ", output[output$Model == "MDiag",]$LogLik <- asr_mdiag$loglik,"\n")
cat("AIC: ",    output[output$Model == "MDiag",]$AIC <- summary.asreml(asr_mdiag)$aic,"\n")
cat("r_g: ",  output[output$Model == "MDiag",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "MDiag",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "MDiag",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")



#-- Model 4: Diagonal
asr_diag <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ diag(env):id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data = met_df,
                   workspace = "1Gb")
asr_diag <- update.asreml(asr_diag)

# Predicted genotype by environment effects
ge_met_pred <- asr_diag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr_diag$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Implicit genotype main effects (averages over environments)
g_met_pred <- rowMeans(ge_met_pred)
# Total genetic effects not calculated

# Store output
cat("LogLik: ", output[output$Model == "Diag",]$LogLik <- asr_diag$loglik,"\n")
cat("AIC: ",    output[output$Model == "Diag",]$AIC <- summary.asreml(asr_diag)$aic,"\n")
cat("r_g: ",  output[output$Model == "Diag",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "Diag",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "Diag",]$r_ge <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")



#-- Model 5: Factor analytic of order 1
asr_fa1 <- asreml(y.Trait1 ~ 1 + env,
                  random   = ~ rr(env, 1):id + diag(env):id + diag(env):block,
                  residual = ~ dsum(~ ar1(col):ar1(row) | env),
                  data = met_df,
                  workspace = "2Gb")
asr_fa1 <- update.asreml(asr_fa1)
asr_fa1 <- update.asreml(asr_fa1)

# Predicted regression genotype by environment effects
ge.reg_met_pred <- asr_fa1$coefficients$random[
  grep("^rr.*id", rownames(asr_fa1$coefficients$random)),]
ge.reg_met_pred <- ge.reg_met_pred[
  grep("Comp", names(ge.reg_met_pred), invert = TRUE)]
ge.reg_met_pred <- matrix(ge.reg_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted residual genotype by environment effects
ge.diag_met_pred <- asr_fa1$coefficients$random[
  grep("^env.*id", rownames(asr_fa1$coefficients$random)),]
ge.diag_met_pred <- matrix(ge.diag_met_pred, nrow = v, ncol = pm, byrow = F)
# Implicit genotype main effects (averages over environments based on latent regression)
g_met_pred <- rowMeans(ge.reg_met_pred)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge.diag_met_pred

# Store output
cat("LogLik: ", output[output$Model == "FA1",]$LogLik <- asr_fa1$loglik,"\n")
cat("AIC: ",    output[output$Model == "FA1",]$AIC <- summary.asreml(asr_fa1)$aic,"\n")
cat("r_g: ",  output[output$Model == "FA1",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "FA1",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "FA1",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")
# Note: The above can be extended to higher order factor analytic models by changing rr(env, 1):id accordingly



#-- View output example
output
#        GEI Envs    Model    LogLik      AIC       r_g       r_m      r_mt      r_ge
# 1 Moderate   20 Expected        NA       NA 0.8649752 0.9103657 0.9501404 0.6387166
# 2 Moderate   20     Main -18349.68 36861.36 0.8507815 0.8619960 0.9596626 0.5438892
# 3 Moderate   20     Comp -17762.12 35688.24 0.8817306 0.9090322 0.9596626 0.7065907
# 4 Moderate   20    MDiag -17197.51 34597.03 0.8333271 0.8303794 0.9596626 0.7234369
# 5 Moderate   20     Diag -17404.73 35009.45 0.8654571 0.9422233 0.9596626 0.6580564
# 6 Moderate   20      FA1 -16775.20 33790.40 0.8124046 0.8591280 0.9596626 0.7215955



# end of script



