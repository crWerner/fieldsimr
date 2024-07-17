#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
cat(" Advance year \n")

#-- Year 7

#Release variety
Release <- selectInd(pop = EYT, nInd = 1)


#-- Year 6

# Perform selections
EYT <- selectInd(pop = AYT, nInd = nEYT)

# FieldSimR code to simulate phenotypes
gv_tpe <- multi_asr_output(pop = EYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepEYT, cov.mat = covs_tpe)
envs_met <- samples_met[[year]][1:nEnvEYT]
gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
error_met <- field_trial_error(nenvs = nEnvEYT, nblocks = nRepEYT, varR = varE_EYT,
                               spatial.model = "AR1", ncols = nColEYT, nrows = nRowEYT)
pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))

EYT@pheno[,1] <- mean_pheno_met

# Report parameters
output <- save_output(output, EYT, "EYT", year, nenvs = nEnvEYT, gv_tpe, pheno_met)


#-- Year 5

AYT <- selectInd(PYT, nInd = nAYT)

# FieldSimR code to simulate phenotypes
gv_tpe <- multi_asr_output(pop = AYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepAYT, cov.mat = covs_tpe)
envs_met <- samples_met[[year]][1:nEnvAYT]
gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
error_met <- field_trial_error(nenvs = nEnvAYT, nblocks = nRepAYT, varR = varE_AYT,
                               spatial.model = "AR1", ncols = nColAYT, nrows = nRowAYT)
pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))

AYT@pheno[,1] <- mean_pheno_met

# Report parameters
output <- save_output(output, AYT, "AYT", year, nenvs = nEnvAYT, gv_tpe, pheno_met)


#-- Year 4

PYT <- selectWithinFam(HDRW, famMax)
PYT <- selectInd(PYT, nInd = nPYT)

# FieldSimR code to simulate phenotypes
gv_tpe <- multi_asr_output(pop = PYT, ntraits = 1, nenvs = nEnvTPE, nreps = nRepPYT, cov.mat = covs_tpe)
envs_met <- samples_met[[year]][1:nEnvPYT]
gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
error_met <- field_trial_error(nenvs = nEnvPYT, nblocks = nRepPYT, varR = varE_PYT,
                               spatial.model = "AR1", ncols = nColPYT, nrows = nRowPYT)
pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))

PYT@pheno[,1] <- mean_pheno_met

# Report parameters
output <- save_output(output, PYT, "PYT", year, nenvs = nEnvPYT, gv_tpe, pheno_met)


#-- Year 3

HDRW <- DH

# FieldSimR code to simulate phenotypes
gv_tpe <- multi_asr_output(pop = HDRW, ntraits = 1, nenvs = nEnvTPE, nreps = nRepHDRW, cov.mat = covs_tpe)
envs_met <- samples_met[[year]][1:nEnvHDRW]
gv_met <- droplevels(gv_tpe[gv_tpe$env %in% envs_met,])
error_met <- field_trial_error(nenvs = nEnvHDRW, nblocks = nRepHDRW, varR = varE_HDRW,
                               spatial.model = "AR1", ncols = nColHDRW, nrows = nRowHDRW)
pheno_met <- make_phenotypes(gv.df = gv_met, error.df = error_met, randomise = TRUE)
mean_pheno_met <- with(pheno_met, tapply(y.Trait1, id, mean))

HDRW@pheno[,1] <- mean_pheno_met

# Report parameters
output <- save_output(output, HDRW, "HDRW", year, nenvs = nEnvHDRW, gv_tpe, pheno_met)


#-- Year 2

DH <- makeDH(F1, nDH)


#-- Year 1

F1 <- randCross(Parents, nCrosses)

#------------------------
rm(envs_met, gv_tpe, gv_met, error_met, pheno_met, mean_pheno_met)