
# multi_asr_input
ntraits <- sample(1:10, 1)
nenvs <- sample(1:100, 1)
corA <- rand_cor_mat(ntraits * nenvs, pos.def = T)
ntraits * nenvs
rankMatrix(corA, tol = 1e-12)
# corA <- cov2cor(eigen(corA)$vectors[,1:5] %*% diag(eigen(corA)$values[1:5]) %*% t(eigen(corA)$vectors[,1:5]))
# eigen(corA)$values
mean <- sample(seq(0,10,0.00001), ntraits * nenvs)
var <- sample(seq(0,10,0.00001), ntraits * nenvs)
(input_asr <- multi_asr_input(ntraits = ntraits,
                              nenvs = nenvs,
                              mean = mean,
                              var = var,
                              corA = corA,
                              nterms = ntraits * nenvs - 5))
plot(mean, input_asr$covs %*% input_asr$mean)

library(AlphaSimR)
FOUNDERPOP <- quickHaplo(nInd = 100,
                         nChr = 10,
                         segSites = 20)

SP <- SimParam$new(FOUNDERPOP)

SP$addTraitA(nQtlPerChr = 20,
             mean = input_asr$mean,
             var = input_asr$var,
             corA = input_asr$corA)

pop <- newPop(FOUNDERPOP)

nreps <- sample(1:4,1)

gv_ls <- multi_asr_output(pop = pop,
                          ntraits = ntraits,
                          nenvs = nenvs,
                          nreps = nreps,
                          covs = input_asr$covs,
                          return.effects = TRUE)

# make phenotypes
str(gv_df_unstr)
tt2 <- gv_df_unstr
tt2$env <- trimws(tt2$env)
tt2$rep <- trimws(tt2$rep)
tt2$id <- trimws(tt2$id)
str(error_df_bivar)
tt3 <- error_df_bivar
tt3$env <- trimws(tt3$env)
tt3$block <- trimws(tt3$block)
tt3$col <- trimws(tt3$col)
tt3$row <- trimws(tt3$row)
str(tt2)
str(tt3)
tt2$id
tt2 <- tt2[sample(rownames(tt2)),]
tt2 <- tt2[order(tt2$env, tt2$rep), ]
head(tt2)

tt2 <- data.frame(tt2, x = 1)
tt3 <- data.frame(tt3, x = 2)
# tt2$x <- NULL
tt4 <- make_phenotypes(gv.df = tt2,
                error.df = tt3,
                randomise = T,
                return.effects = T)
head(tt4$pheno.df)
head(tt4$gv.df)
head(tt4$error.df)
plot(tt4$gv.df$gv.trait1[order(tt4$gv.df$env, tt4$gv.df$rep, tt4$gv.df$id)], tt2$gv.Trait1[order(tt2$env, tt2$rep, as.numeric(tt2$id))])
plot(tt4$gv.df$gv.trait2[order(tt4$gv.df$env, tt4$gv.df$rep, tt4$gv.df$id)], tt2$gv.Trait2[order(tt2$env, tt2$rep, as.numeric(tt2$id))])
plot(tt4$pheno.df$phe.Trait1, tt4$gv.df$gv.trait1 + tt3$e.Trait1)
plot(tt4$pheno.df$phe.Trait2, tt4$gv.df$gv.trait2 + tt3$e.Trait2)

plot(tt4$phe.Trait2, tt2$gv.Trait2 + tt3$e.Trait2)
tt4$error.df
str(tt4)
tt4$id
tt2$id
plot(tt4$phe.Trait1, tt2$gv.Trait1 + tt3$e.Trait1)
plot(tt4$phe.Trait2, tt2$gv.Trait2 + tt3$e.Trait2)

tt5 <- make_phenotypes(gv.df = gv_df_unstr,
                       error.df = error_df_bivar,
                       randomise = F,
                       return.effects = F)
head(tt5)
head(tt5$gv.df)
head(tt5$error.df)

plot(tt5$pheno.df$phe.Trait1, tt5$gv.df$gv.trait1 + error_df_bivar$e.Trait1)
plot(tt5$pheno.df$phe.Trait2, gv_df_unstr$gv.Trait2 + error_df_bivar$e.Trait2)

plot(tt5$pheno.df$phe.Trait1, gv_df_unstr$gv.Trait1 + error_df_bivar$e.Trait1)
plot(tt5$pheno.df$phe.Trait2, gv_df_unstr$gv.Trait2 + error_df_bivar$e.Trait2)

plot(tt5$gv.df$gv.trait1, gv_df_unstr$gv.Trait1)
plot(tt5$gv.df$gv.trait2, gv_df_unstr$gv.Trait2)

plot(tt5$gv.df$gv.trait1[order(tt5$gv.df$env, tt5$gv.df$rep, tt5$gv.df$id)], gv_df_unstr$gv.Trait1[order(gv_df_unstr$env, gv_df_unstr$rep, as.numeric(gv_df_unstr$id))])
plot(tt5$gv.df$gv.trait2[order(tt5$gv.df$env, tt5$gv.df$rep, tt5$gv.df$id)], gv_df_unstr$gv.Trait2[order(gv_df_unstr$env, gv_df_unstr$rep, as.numeric(gv_df_unstr$id))])

str(tt5)
plot(tt2$gv.Trait1, gv_df_unstr$gv.Trait1)
plot(tt2$gv.Trait2, gv_df_unstr$gv.Trait2)
plot(error_df_bivar$e.Trait1, error_df_bivar$e.Trait1)
plot(error_df_bivar$e.Trait2, error_df_bivar$e.Trait2)

plot(trimws(tt4$id), tt5$id)
plot(tt4$pheno.df$phe.Trait1, tt5$phe.Trait1)
plot(tt4$pheno.df$phe.Trait2, tt5$phe.Trait2)


# rand_cor_mat
tt <- rand_cor_mat(n = 100,
                   min.cor = 0,
                   pos.def = T,
                   small.positive = 1e-8)
mean(tt[upper.tri(tt)])
solve(tt)
svd(tt)$d

# compsym_asr_input
ntraits <- sample(1:100,1)
nenvs <- sample(1:100,1)
(prop.mainA <- sample(seq(0,1,0.001),1))
(prop.mainDD <- sample(seq(0,1,0.001),1))
(relAA <- sample(seq(0,10,0.001),1))
(prop.mainAA <- sample(seq(0,1,0.001),1))

input_asr <- compsym_asr_input(ntraits = ntraits,
                               nenvs = nenvs,
                               prop.mainA = prop.mainA,
                               meanDD = 2,
                               varDD = 1,
                               prop.mainDD = prop.mainDD,
                               relAA = relAA,
                               prop.mainAA = prop.mainAA
                              )
input_asr$mean
input_asr$var
input_asr$corA[1:10,1:10]
input_asr$meanDD
input_asr$varDD
input_asr$corDD[1:10,1:10]
input_asr$relAA
input_asr$corAA[1:10,1:10]

library(AlphaSimR)
FOUNDERPOP <- quickHaplo(nInd = 100,
                         nChr = 10,
                         segSites = 20)

SP <- SimParam$new(FOUNDERPOP)

SP$addTraitA(nQtlPerChr = 20,
             mean = input_asr$mean,
             var = input_asr$var,
             corA = input_asr$corA)

pop <- newPop(FOUNDERPOP)

nreps <- sample(1:4,1)

# compsym_asr_output

gv_ls <- compsym_asr_output(pop = pop,
                            ntraits = ntraits,
                            nenvs = nenvs,
                            nreps = nreps,
                            return.effects = TRUE)
str(gv_ls$gv.df)
gv_ls$gv.df$env
head(gv_ls$gv.df)
head(gv_ls$Trait1)
head(gv_ls$Trait2)
