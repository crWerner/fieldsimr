
#######################################
#
# Supplementary Script S2
#
# Working example: Application to the two
# hypothetical environments with fitting
# in ASReml-R
#
# Script Author: D.J. Tolhurst
#
#######################################


require(devtools)
require(FieldSimR)
require(ggplot2)
require(asreml)


# Firstly, we will add some noise to the true genetic effects
# by simulating spatial variation with FieldSimR
# assume an RCBD with 2 blocks
nreps <- nblocks <- 2
nplots <- nblocks*nids
ncols <- 2
nrows <- nplots/ncols

set.seed(123)
e_df <- field_trial_error(nenvs = 2,
                          nblocks = 2,
                          block.dir = "col",
                          ncols = ncols,
                          nrows = nrows,
                          varR = 3) # so plot-level h2 = 1/(1+3) = 0.25 and 3/(3+3) = 0.5
str(e_df)
FieldSimR::plot_effects(e_df[e_df$env == 1,], effect = "e.Trait1", blocks = T)
FieldSimR::plot_effects(e_df[e_df$env == 2,], effect = "e.Trait1", blocks = T)

# Construct true genetic effects
ge_df <- data.frame(env = factor(rep(1:nenvs, each = nids*nreps)),
                    id = factor(1:nids),
                    rep = factor(rep(1:nreps, each = nids)),
                    gv.Trait1 = rep(c(4,8), each = nids*nreps) + # add env means
                                c(ge_mat[rep(1:nids, 2),]))
head(ge_df)
with(ge_df, tapply(gv.Trait1, env, mean)) # 4, 8 means

# construct phenotypes
pheno_df <- make_phenotypes(ge_df, e_df)
head(pheno_df)

# check plot-level h2
tt <- with(pheno_df, tapply(y.Trait1, list(block, id, env), function(x) unique(x)))
gg <- with(ge_df, tapply(gv.Trait1, list(rep, id, env), function(x) unique(x)))
plot(tt[,,1], gg[,,1])
cor(c(tt[,,1]), c(gg[,,1]))^2 # 0.2164983 (target was 0.25)
plot(tt[,,2], gg[,,2])
cor(c(tt[,,2]), c(gg[,,2]))^2 # 0.541423  (target was 0.5)

# check line-mean r2 within envs
tt <- with(pheno_df, tapply(y.Trait1, list(id, env), function(x) mean(x)))
gg <- with(ge_df[ge_df$rep == 1,], tapply(gv.Trait1, list(id, env), function(x) unique(x)))
plot(tt[,1], gg[,1])
cor(c(tt[,1]), c(gg[,1]))^2 # 0.5440478 (target was 1/(1+3/2) = 0.4)
plot(tt[,2], gg[,2])
cor(c(tt[,2]), c(gg[,2]))^2 # 0.6610483  (target was 3/(3+3/2) = 0.6666667)
# overall
tt <- with(pheno_df, tapply(y.Trait1, id, function(x) mean(x)))
gg <- with(ge_df[ge_df$rep == 1,], tapply(gv.Trait1, id, function(x) mean(x)))
plot(tt, gg)
cor(tt, gg)^2 # 0.6588496

# Now we move to asreml for a reaction norm model,
# i.e., a random regression model with correlated intercepts and slopes
mbf_df <- data.frame(env = factor(1:nenvs),
                     int = 1,
                     cov = cov_mat)
asr <- asreml(y.Trait1 ~ env + env:block, # fixed block effects for the RCBD
              random = ~ str(~mbf(S):id, ~us(2):id(id)),
              residual = ~dsum(~id(col):ar1(row)|env),
              mbf = list(S = list(key = c("env", "env"), cov = "mbf_df")),
              data = pheno_df)
asr <- update(asr)
asr <- update(asr)
asr <- update(asr)
# -91.35956
summary(asr)$varcom

# variances
rr_mat <- matrix(NA, ncol = 2, nrow = 2)
rr_mat[upper.tri(rr_mat, diag = T)] <- asr$vparameters[grep("mbf", names(asr$vparameters))]
rr_mat[2,1] <- rr_mat[1,2]
rr_mat

# main effect and interaction variances
cov_mat <- as.matrix(mbf_df[,-1])
G_mat <- cov_mat %*% rr_mat %*% t(cov_mat)
G_mat # genetic variance and covariances among environments
mean(G_mat) # main effect variance
rr_mat[1,1] # equivalent to intercept variance
# 1.345293
# GxE variance within environments
mean(diag(G_mat - mean(G_mat)))
cbind(cov_mat[,2]) %*% rr_mat[2,2] %*% t(cov_mat[,2])
# 1.286667

# intercepts and slopes
(int <- asr$coefficients$random[grep("int", rownames(asr$coefficients$random))])
plot(int, g_mat) # with truth
cor(int, g_mat) # 0.8276225
(slope <- asr$coefficients$random[grep("cov", rownames(asr$coefficients$random))])
plot(slope, slope_mat) # with truth
cor(slope, slope_mat) # 0.8088664

# total effects
tot_eff <- cbind(int, slope) %*% t(cov_mat)
diag(cor(tot_eff, with(ge_df[ge_df$rep == 1,], tapply(gv.Trait1, list(id, env), mean)))) # with truth
# 0.7775920 0.8680046

# obtain non-crossover and crossover effects
(br <- rr_mat[2,2]/rr_mat[1,1]) # regression slope between intercepts and slopes
# non-crossover effects
ncross_eff <- int %*% rbind(1+ cov_mat[,2]*br)
plot_table(ncross_eff)
plot_effects(ncross_eff)
# crossover effects
cross_eff <- tot_eff - ncross_eff
slope %*% t(cov_mat[,2]) - int %*% t(cov_mat[,2]*br) # equivalence
plot_table(cross_eff)
plot_effects(cross_eff)

# obtain PEV matrix for the non-crossover and crossover effects
# consider the project given by
pev <- predict(asr,
               classify = "env:id",
               only = "mbf(S):id",
               vcov = TRUE, maxit = 1)
# check
range(pev$pvals$predicted.value - c(int, slope))
# ~0
(pev <- as.matrix(pev$vcov)) # nids*2 x nids*2

# pev for non-crossover effects
# consider the rotation given by
rot_ncross <- rbind(1 + cov_mat[,2]*br)
pev_int <- pev[1:nids, 1:nids]
K_ncross <- kronecker(t(rot_ncross), diag(nids))
pev_ncross <- K_ncross %*% pev_int %*% t(K_ncross)

# and we can see the mean pev across non-crossover effects is
pev_ncross_mean <- kronecker(t(c(1,1)), diag(nids)) %*% pev_ncross %*% kronecker(c(1,1), diag(nids))/nenvs^2
range(pev_ncross_mean - pev_int)
# which equals the pev for the intercepts, by design

# pev for crossover effects
# consider the rotation given by
rot_cross <- rbind(-cov_mat[, 2] * br, cov_mat[, 2])
K_cross <- kronecker(t(rot_cross), diag(nids))
pev_cross <- K_cross %*% pev %*% t(K_cross)

# compare to pev for interaction effects
K_gxe <- kronecker(cov_mat[,2], diag(nids))
pev_gxe <- K_gxe %*% pev[21:40,21:40] %*% t(K_gxe)
plot(pev_cross, pev_gxe); abline(a=0, b=1)
# those below the 1:1 line are for env2, where the crossover effects are
# smaller than interaction effects
# and vice versa for env1


#######################################
# end of script

save.image()

