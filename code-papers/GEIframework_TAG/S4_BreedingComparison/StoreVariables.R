# Store variables

# Set up dataframe to store simulation variables
output = data.frame(year = rep(1:nCycles, times = nStages),
                    rep  = numeric(nCycles*nStages),
                    scenario    = Scenario,
                    GEI         = GEI,
                    stage       = rep(c("Parents","HDRW","PYT","AYT","EYT"),
                                      each = nCycles),
                    mean_gv_tpe = numeric(nCycles*nStages),
                    mean_gv_met = numeric(nCycles*nStages),
                    var_gv_tpe = numeric(nCycles*nStages),
                    var_gv_met = numeric(nCycles*nStages),
                    acc_tpe = numeric(nCycles*nStages),
                    acc_met = numeric(nCycles*nStages),
                    met_tpe = numeric(nCycles*nStages))


# Function to automate saving output
save_output <- function(df, pop, stage, year, nenvs = NULL, gv_tpe, pheno_met = NULL) {
  # Note: if nenvs not specified, then MET values returned as NA
  # Note: if pheno_met not specified, then phenotypic values returned as NA
  require(data.table)

  # Genotype main effects
  # gv_means <- with(gv_tpe, tapply(gv.Trait1, list(id,env), mean)) # slow
  setDT(gv_tpe)
  gv_means <- gv_tpe[, .(mean_gv = mean(gv.Trait1)), by = .(env, id)]
  gv_means <- matrix(gv_means$mean_gv, nrow = pop@nInd, ncol = uniqueN(gv_tpe$env))
  mean_gv_tpe <- rowMeans(gv_means)
  mean_gv_met <- if (is.null(nenvs)) NA
                 else if (nenvs == 1) gv_means[, samples_met[[year]][1]]
                 else rowMeans(gv_means[, samples_met[[year]][1:nenvs]])

  # phenotypic main effects when available
  if (is.null(pheno_met)) {
    mean_pheno_met <- NA
  } else {
    pheno_means <- with(pheno_met, tapply(y.Trait1, list(id, env), mean))
    mean_pheno_met <- if (is.null(nenvs)) NA
                      else if (nenvs == 1) pheno_means[,1]
                      else rowMeans(pheno_means[,1:nenvs])
  }

  index <- which(df$year == year & df$stage == stage)
  df[index, c("mean_gv_tpe", "mean_gv_met")] <- c(mean(mean_gv_tpe), mean(mean_gv_met))
  df[index, c("var_gv_tpe", "var_gv_met")]   <- c(var(mean_gv_tpe), var(mean_gv_met))
  df[index, c("acc_tpe", "acc_met", "met_tpe")] <-
    c(if (is.null(pheno_met)) NA else cor(mean_gv_met, mean_pheno_met),
      if (is.null(pheno_met)) NA else cor(mean_gv_tpe, mean_pheno_met),
      if (is.null(nenvs))     NA else cor(mean_gv_met, mean_gv_tpe))

  return(df)
}
