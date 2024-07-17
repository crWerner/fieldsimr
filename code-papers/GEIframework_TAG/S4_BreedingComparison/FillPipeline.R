# Populate breeding stages with unique genotypes

for(year in 1:7){
  cat("FillPipeline year:", year, "of 7\n")
  if(year<8){
    #Year 1
    F1 <- randCross(Parents, nCrosses)
  }
  if(year<7){
    #Year 2
    DH <- makeDH(F1, nDH)
  }
  if(year<6){
    #Year 3
    
    HDRW <- DH
    
    # FieldSimR code to simulate phenotypes
    gv_tpe <- multi_asr_output(pop = HDRW, ntraits = 1, nenvs = nEnvTPE, nreps = nRepHDRW, cov.mat = covs_tpe)
    envs_met <- samples_met[[4]][1:nEnvHDRW]
    gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
    error_met <- field_trial_error(nenvs = nEnvHDRW, nblocks = nRepHDRW, varR = varE_HDRW,
                                   spatial.model = "AR1", ncols = nColHDRW, nrows = nRowHDRW)
    pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
    mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
    mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
    
    HDRW@pheno[,1] <- mean_pheno_met
    
    cat(" HDRW h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
  }
  if(year<5){
    #Year 4
  
    PYT <- selectWithinFam(HDRW, famMax)
    PYT <- selectInd(PYT, nInd = nPYT)
    
    # FieldSimR code to simulate phenotypes
    gv_tpe <- multi_asr_output(pop = PYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepPYT, cov.mat = covs_tpe)
    envs_met <- samples_met[[3]][1:nEnvPYT]
    gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
    error_met <- field_trial_error(nenvs = nEnvPYT, nblocks = nRepPYT, varR = varE_PYT,
                                   spatial.model = "AR1", ncols = nColPYT, nrows = nRowPYT)
    pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
    mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
    mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
    
    PYT@pheno[,1] <- mean_pheno_met
    
    cat(" PYT h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
  }
  if(year<4){
    #Year 5
    
    AYT = selectInd(PYT, nInd = nAYT)
    
    # FieldSimR code to simulate phenotypes
    gv_tpe <- multi_asr_output(pop = AYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepAYT, cov.mat = covs_tpe)
    envs_met <- samples_met[[2]][1:nEnvAYT]
    gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
    error_met <- field_trial_error(nenvs = nEnvAYT, nblocks = nRepAYT, varR = varE_AYT,
                                   spatial.model = "AR1", ncols = nColAYT, nrows = nRowAYT)
    pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
    mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
    mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
    
    AYT@pheno[,1] <- mean_pheno_met
    
    cat(" AYT h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
  }
  if(year<3){
    #Year 6
    
    EYT = selectInd(pop = AYT, nInd = nEYT)
    
    # FieldSimR code to simulate phenotypes
    gv_tpe <- multi_asr_output(pop = EYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepEYT, cov.mat = covs_tpe)
    envs_met <- samples_met[[1]][1:nEnvEYT]
    gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
    error_met <- field_trial_error(nenvs = nEnvEYT, nblocks = nRepEYT, varR = varE_EYT,
                                   spatial.model = "AR1", ncols = nColEYT, nrows = nRowEYT)
    pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
    mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))
    mean_gv_met <- with(gv_met, tapply(gv.Trait1, id, mean))
    
    EYT@pheno[,1] <- mean_pheno_met
    
    cat(" EYT h2 =", cor(mean_pheno_met, mean_gv_met)^2, "\n", sep = " ")
  }
  if(year<2){
    #Year 7
    Release = selectInd(pop = EYT, nInd = 1)
    cat(" Release", "\n")
  }
}

#------------------------
rm(envs_met, gv_tpe, gv_met, mean_gv_met, error_met, pheno_met, mean_pheno_met)