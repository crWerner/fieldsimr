

###################################################################
#
# Supplementary Script S11
#
# Linear mixed models fitted in Werner, Gemenet, and Tolhurst (2024)
#
# Script author: D.J. Tolhurst
#
###################################################################


library(FieldSimR)
library(asreml)
asreml.options(Cfixed = TRUE) # required to calculate AIC
library(ASExtras4) # available from https://mmade.org/tpsbits/
library(TPSbits) # available from https://mmade.org/tpsbits/
library(SpATS)
library(mvngGrAd)
library(ggplot2)

# data frames ----------
# The results in this script can be reproduced using the supplied object 'pheno_df'
load("pheno_df.RData")

# Grain yield phenotypes in Environment 1
df <- droplevels(pheno_df[pheno_df$env == 1,])
head(df)
dim(df)
#  200   6
# factorise
df$env <- factor(as.numeric(trimws(df$env)))
df$block <- factor(as.numeric(trimws(df$block)))
df$col <- factor(as.numeric(trimws(df$col)))
df$row <- factor(as.numeric(trimws(df$row)))
df$id <- factor(as.numeric(trimws(df$id)))
df <- df[order(df$col, df$row),]

# obtain true genetic values
true_gvs <- droplevels(gv_df_unstr[gv_df_unstr$env == 1,])
head(true_gvs)
dim(true_gvs)
# 200   5
# factorise
true_gvs$env <- factor(as.numeric(trimws(true_gvs$env)))
true_gvs$rep <- factor(as.numeric(trimws(true_gvs$rep)))
true_gvs$id <- factor(as.numeric(trimws(true_gvs$id)))
true_gvs <- true_gvs[order(true_gvs$env, true_gvs$rep, true_gvs$id),]


# Model fitting  ----------

# Note: users are encouraged to perform data cleaning and experimental design checks prior to model fitting,
# in addition to diagnostic checks during model fitting.


# Baseline model: ID error model, equivalent to a complete block analysis
id.asr0 <- asreml(y.Trait1 ~ 1,
                  random = ~ id + block,
                  residual = ~ idv(units),
                  data = df)
id.asr0 <- update(id.asr0)
# 35.21
# Residual deviance ---
-2 * id.asr0$loglik
# -70.43
# AIC ---
print(icREML(list(id.asr0)))
# -66.72
# Note: this is based on the full log-likelihood approach of Verbyla (2019)
# The icREML function can be obtained from https://doi.org/10.1111/anzs.12254
# Prediction accuracy ---
cor(id.asr0$coefficients$random[grep("id", rownames(id.asr0$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.652
# Estimated variance parameters
summary(id.asr0)$varcom
# bias = 0.061 - 0.087 = -0.026
# Note: the bias is calculated as the absolute difference between the estimated genetic variance
# and the true genetic variance of 0.087
# Some additional diagnostics
library(ASExtras4)
# ASExtras4 is available from https://mmade.org/tpsbits/
plot(id.asr0)


# Model 1: AR1 model
ar1.asr1 <- asreml(y.Trait1 ~ 1,
                   random = ~ id + block,
                   residual = ~ ar1(col):ar1(row),
                   data = df)
ar1.asr1 <- update(ar1.asr1)
# 49.53
# Residual deviance ---
-2 * ar1.asr1$loglik
# -99.06
# REMLRT
pchisq(-2*(id.asr0$loglik - ar1.asr1$loglik),2,lower.tail = F)
lrt(id.asr0, ar1.asr1, boundary = F) # standard REMLRT for testing correlation parameters
# p < 0.0001
# AIC ---
print(icREML(list(ar1.asr1)))
# -94.41
# Prediction accuracy ---
cor(ar1.asr1$coefficients$random[grep("id", rownames(ar1.asr1$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.722
# Estimated variance parameters
summary(ar1.asr1)$varcom
# bias = 0.073 - 0.087 = -0.014
# Some plots
vario.df1 <- data.frame(df[,c("col", "row")],
                        e.spat = c(ar1.asr1$residuals))
# Sample variogram
sample_variogram(vario.df1,
                 effect = "e.spat")
# Theoretical variogram
theoretical_variogram(ncols = 10, nrows = 20, varR = 0.20, col.cor = 0.51, row.cor = 0.23, prop.spatial = 1)
# Trellis plot
ggplot(vario.df1, aes(x = row, y = e.spat, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
# Strong evidence of row effects, but lets fit an additional ID error term first
# QQ plot
qq_plot(vario.df1, effect = "e.spat", labels = TRUE)
# Some additional diagnostics
plot(ar1.asr1)


# Model 2: AR1 + ID model
ar1.asr2 <- asreml(y.Trait1 ~ 1,
                   random = ~ id + block + units,
                   residual = ~ ar1(col):ar1(row),
                   data = df)
ar1.asr2 <- update(ar1.asr2)
# 56.87
# Residual LogLik ---
-2 * ar1.asr2$loglik
# -113.74
# REMLRT
pchisq(-2*(ar1.asr1$loglik - ar1.asr2$loglik),1,lower.tail = F)/2 # REMLRT based on the non-zero variance approach of Self and Liang (1987), Stram and Lee (1994)
lrt(ar1.asr1, ar1.asr2)
# p < 0.0001
# AIC ---
print(icREML(list(ar1.asr2)))
# -103.69
# Prediction accuracy ---
cor(ar1.asr2$coefficients$random[grep("id", rownames(ar1.asr2$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.742
# Estimated variance parameters
summary(ar1.asr2)$varcom
# bias = 0.066 - 0.087 = -0.021
# Some plots
vario.df2 <- data.frame(df[,c("col", "row")],
                        e.spat = c(ar1.asr2$residuals),
                        e.rand = ar1.asr2$coefficients$random[grep("units", rownames(ar1.asr2$coefficients$random)),])
vario.df2$e.total <- vario.df2$e.spat + vario.df2$e.rand
# Sample variogram
sample_variogram(vario.df2,
                 effect = "e.spat")
sample_variogram(vario.df2,
                 effect = "e.total")
# Strong evidence of row effects, we will fit row terms to account for this extraneous variation
# Theoretical variogram
theoretical_variogram(ncols = 10, nrows = 20, varR = 0.49, prop.spatial = 0.82, col.cor = 0.95, row.cor = 0.87)
# Trellis plot
ggplot(vario.df2, aes(x = row, y = e.spat, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
ggplot(vario.df2, aes(x = row, y = e.total, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
# QQ plot
qq_plot(vario.df2, effect = "e.spat", labels = TRUE)
qq_plot(vario.df2, effect = "e.total", labels = TRUE)
# Some additional diagnostics
plot(ar1.asr2)


# Model 3: AR1 + ID model, with fixed and random row
df$Frow <- 1
df$Frow[df$row %in% seq(1,20,2)] <- 2
df$Frow <- factor(as.numeric(trimws(df$Frow)))
ar1.asr3 <- asreml(y.Trait1 ~ 1 + Frow,
                   random = ~ id + block + row + units,
                   residual = ~ ar1(col):ar1(row),
                   data = df)
ar1.asr3 <- update(ar1.asr3)
# 66.71
# Fixed effects tests ---
wald(ar1.asr3)
wald(ar1.asr3, denDF = "algebraic") # Wald F-test with denominator degrees of freedom following Kenward and Roger (1997)
# p < 0.001 significant
# Residual deviance ---
-2 * ar1.asr3$loglik
# -133.43
# AIC ---
print(icREML(list(ar1.asr3)))
# -126.60
# Prediction accuracy ---
cor(ar1.asr3$coefficients$random[grep("id", rownames(ar1.asr3$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.757
# Estimated variance parameters
summary(ar1.asr3)$varcom
# bias = 0.075 - 0.087 = -0.012
# Some plots
ar1.asr3$coefficients$fixed
# Note: the fixed row effects are non-estimable, since they are confounded with the overall mean
# We will impose a constraint where they must sum-to-zero (for demonstration purposes)
ar1.asr3$coefficients$fixed[1:2] + ar1.asr3$coefficients$fixed[3] - mean(ar1.asr3$coefficients$fixed[1:2] + ar1.asr3$coefficients$fixed[3])
# -0.1715097  0.1715097
vario.df3 <- data.frame(df[,c("col", "row")],
                        e.spat = c(ar1.asr3$residuals),
                        e.rand = ar1.asr3$coefficients$random[grep("units", rownames(ar1.asr3$coefficients$random)),],
                        e.ext.row.fix = c(0.1715104, -0.1715104),
                        e.ext.row.rand = ar1.asr3$coefficients$random[grep("row", rownames(ar1.asr3$coefficients$random)),])
vario.df3$e.total1 <- vario.df3$e.spat + vario.df3$e.rand
vario.df3$e.total2 <- vario.df3$e.spat + vario.df3$e.rand + vario.df3$e.ext.row.fix + vario.df3$e.ext.row.rand
# Sample variogram
ggplot(vario.df3, aes(x = row, y = e.ext.row.fix + e.ext.row.rand, group = col)) + geom_line(linewidth = 0.5, colour = "steelblue") + geom_point(size = 2)
sample_variogram(vario.df3,
                 effect = "e.spat")
# There is not a linear trend in the column direction, the non-stationary
# variation observed here is because the ar1 process has not plateaued yet.
# We demonstrate this in the manuscript using the faces of the sample variogram,
# which show that the trend in the column direction is well within the 95% coverage
# intervals expected from the stochastic nature of the autoregressive process.
sample_variogram(vario.df3,
                 effect = "e.total1")
sample_variogram(vario.df3,
                 effect = "e.total2")
# Theoretical variogram
theoretical_variogram(ncols = 10, nrows = 20, varR = 0.18, prop.spatial = 0.57, col.cor = 0.75, row.cor = 0.89)
# Trellis plot
ggplot(vario.df3, aes(x = row, y = e.spat, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
ggplot(vario.df3, aes(x = row, y = e.total1, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
# QQ plot
qq_plot(vario.df3, effect = "e.spat", labels = TRUE)
qq_plot(vario.df3, effect = "e.total1", labels = TRUE)
qq_plot(vario.df3, effect = "e.total2", labels = TRUE)
# Some additional diagnostics
plot(ar1.asr3)


# Model 4: TPS + ID model
# First we fit in the SpATS package using the PSANOVA functionality,
# with 6 knots in the column direction and 12 knots in the row direction
library(SpATS)
df$col.num <- as.numeric(trimws(df$col))
df$row.num <- as.numeric(trimws(df$row))
tps4.spats <- SpATS(response = "y.Trait1",
                    fixed = ~ 1,
                    random = ~ block,
                    spatial = ~ PSANOVA(col.num, row.num, nseg = c(5,11), degree = 3, pord = 2),
                    genotype = "id",  genotype.as.random = TRUE,
                    data = df,
                    control =  list(tolerance = 1e-06))
tps4.spats$var.comp
plot(tps4.spats)
plot(variogram.SpATS(tps4.spats))

# Next we obtain identical results in ASReml-R with the help of
# the TPSbits helper functions tpsmmb, tpsfitted, and tpsurface
# TPSbits is available from https://mmade.org/tpsbits/
TPS.df <- tpsmmb("col.num", "row.num", df, stub="1", nsegments=c(5,11)) # to obtain the necessary basis functions
BcZ.df <- TPS.df$BcZ.df
BrZ.df <- TPS.df$BrZ.df
BcrZ.df <- TPS.df$BcrZ.df
tps4.asr <- asreml(y.Trait1 ~ 1 + TP.CR.2 + TP.CR.3 + TP.CR.4,
                   random= ~ id + block +
                             TP.C.1:mbf(TP.row) + TP.C.2:mbf(TP.row) +
                             TP.R.1:mbf(TP.col) + TP.R.2:mbf(TP.col) + mbf(TP.CxR),
                   mbf = list(TP.col = list(key = c("TP.col","TP.col"), cov = "BcZ.df"),
                              TP.row = list(key = c("TP.row","TP.row"), cov = "BrZ.df"),
                              TP.CxR = list(key = c("TP.CxR","TP.CxR"), cov = "BcrZ.df")),
                   data=TPS.df$data)
tps4.asr <- update(tps4.asr)
# 55.39
# Residual deviance ---
-2 * tps4.asr$loglik
# -110.789
# AIC ---
print(icREML(list(tps4.asr)))
# -93.8
# Prediction accuracy ---
cor(tps4.asr$coefficients$random[grep("id", rownames(tps4.asr$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.687
# Estimated variance parameters
summary(tps4.asr)$varcom
# bias = 0.037 - 0.087 = -0.050
# Some plots
vario.df4 <- data.frame(df[,c("col", "row")],
                        e.rand = c(tps4.asr$residuals))
# Sample variogram
sample_variogram(vario.df4,
                 effect = "e.rand")
# TPS surface plot
tpsurface(tpsfitted(tps4.asr, TPS.df))
# Trellis plot
ggplot(vario.df4, aes(x = row, y = e.rand, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
# Strong evidence of row effects
# QQ plot
qq_plot(vario.df4, effect = "e.rand", labels = TRUE)
# Some additional diagnostics
plot(tps4.asr)


# Model 5: TPS + ID model, with random column amd row terms
# This model is equivalent to the SPaTS approach of Rodriguez-Alvarez et al. (2018)
# First we fit in the SpATS package, using the same knot points as above
tps5.spats <- SpATS(response = "y.Trait1",
                    fixed = ~ 1,
                    random = ~ block + row + col,
                    spatial = ~ PSANOVA(col.num, row.num, nseg = c(5,11), degree = 3, pord = 2),
                    genotype = "id",  genotype.as.random = TRUE,
                    data = df,
                    control =  list(tolerance = 1e-06))
tps5.spats$var.comp
plot(tps5.spats)
plot(variogram.SpATS(tps5.spats))

# Next we obtain identical results in ASReml-R, using the
# same basis functions as above
tps5.asr <- asreml(y.Trait1 ~ 1 + TP.CR.2 + TP.CR.3 + TP.CR.4,
                   random= ~ id + block + col + row +
                             TP.C.1:mbf(TP.row) + TP.C.2:mbf(TP.row) +
                             TP.R.1:mbf(TP.col) + TP.R.2:mbf(TP.col) + mbf(TP.CxR),
                   mbf = list(TP.col = list(key = c("TP.col","TP.col"), cov = "BcZ.df"),
                              TP.row = list(key = c("TP.row","TP.row"), cov = "BrZ.df"),
                              TP.CxR = list(key = c("TP.CxR","TP.CxR"), cov = "BcrZ.df")),
                   data=TPS.df$data)
tps5.asr <- update(tps5.asr)
# 67.8372
# Residual deviance ---
-2 * tps5.asr$loglik
# -135.67
# REMLRT
pchisq(-2*(tps4.asr$loglik - tps5.asr$loglik),2,lower.tail = F)/2 + pchisq(-2*(tps4.asr$loglik - tps5.asr$loglik),1,lower.tail = F)/2
lrt(tps4.asr, tps5.asr)
# <0.00001
# AIC ---
print(icREML(list(tps5.asr)))
# -116.66
# Prediction accuracy
cor(tps5.asr$coefficients$random[grep("id", rownames(tps5.asr$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.720
# Estimated variance parameters
summary(tps5.asr)$varcom
# bias = 0.037 - 0.087 = -0.021
# Some plots
vario.df5 <- data.frame(df[,c("col", "row")],
                        e.rand = c(tps5.asr$residuals),
                        e.ext.row = c(tps5.asr$coefficients$random[grep("row_", rownames(tps5.asr$coefficients$random)),]),
                        e.ext.col = c(tps5.asr$coefficients$random[grep("col_", rownames(tps5.asr$coefficients$random)),]))
vario.df5$e.total <- vario.df5$e.rand + vario.df5$e.ext.row + vario.df5$e.ext.col
# Sample variogram
sample_variogram(vario.df5,
                 effect = "e.rand")
sample_variogram(vario.df5,
                 effect = "e.total")
# TPS surface plot
tpsurface(tpsfitted(tps5.asr, TPS.df))
# Trellis plot
ggplot(vario.df5, aes(x = row, y = e.rand, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
ggplot(vario.df5, aes(x = row, y = e.total, group = col)) + geom_point() + facet_wrap(~col) + geom_line(linewidth = 0.25, colour = "steelblue")
# QQ plot
qq_plot(vario.df5, effect = "e.rand", labels = TRUE)
qq_plot(vario.df5, effect = "e.total", labels = TRUE)
# Some additional diagnostics
plot(tps5.asr)


# Model 6: Nearest neighbour adjustment with 0 additional layers
library(mvngGrAd)
# This grid is taken from the mvngGrAd vignette
cross <- list(down = 1,
              up = 1,
              left = 1:4,
              right = 1:4)
# Plot the grid
sketchGrid(i = 10, rowLimit = 20,
           j = 5, colLimit = 10,
           shapeCross = cross,
           excludeCenter = TRUE,
           layers = NULL)
# Adjust the phenotypes
nn_adjustment0 <- movingGrid(rows = df$row.num, # requires numeric numbers
                             columns = df$col.num,
                             obsPhe = df$y.Trait1,
                             shapeCross = cross,
                             layers = NULL,
                             excludeCenter = TRUE)
df$y.Trait1.adjusted0 <- nn_adjustment0@adjustedPhe
# Fit the model in ASReml-R
nn.asr6 <- asreml(y.Trait1.adjusted0 ~ 1,
                  random = ~ id + block,
                  residual = ~ idv(col):id(row),
                  data = df)
nn.asr6 <- update(nn.asr6)
# 51.88
# Residual deviance ---
-2 * nn.asr6$loglik
# -103.76
# AIC ---
print(icREML(list(nn.asr6)))
# -102.0
# Prediction accuracy ---
cor(nn.asr6$coefficients$random[grep("id", rownames(nn.asr6$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.709
# Estimated variance parameters
summary(nn.asr6)$varcom
# bias = 0.084 - 0.087 = -0.003
# Some additional diagnostics
plot(nn.asr6)


# Model 7: Nearest neighbour adjustment with 1 additional layer
# This grid is an adaptation of the mvngGrAd vignette
# Plot the grid
sketchGrid(i = 10, rowLimit = 20,
           j = 5, colLimit = 10,
           shapeCross = cross,
           excludeCenter = TRUE,
           layers = 1)
# Adjust the phenotypes
nn_adjustment1 <- movingGrid(rows = df$row.num,
                             columns = df$col.num,
                             obsPhe = df$y.Trait1,
                             shapeCross = cross,
                             layers = 1,
                             excludeCenter = TRUE)
df$y.Trait1.adjusted1 <- nn_adjustment1@adjustedPhe
# Fit the model in ASReml-R
nn.asr7 <- asreml(y.Trait1.adjusted1 ~ 1,
                  random = ~ id + block,
                  residual = ~ idv(units),
                  data = df)
nn.asr7 <- update(nn.asr7)
# 47.11
# Residual deviance ---
-2 * nn.asr7$loglik
# -94.22
# AIC ---
print(icREML(list(nn.asr7)))
# -94.74
# Prediction accuracy ---
cor(nn.asr7$coefficients$random[grep("id", rownames(nn.asr7$coefficients$random)),], with(true_gvs, tapply(gv.Trait1, id, mean)))
# 0.702
# Estimated variance parameters
summary(nn.asr7)$varcom
# bias = 0.064 -0.087 = -0.023
# Some additional diagnostics
plot(nn.asr7)


save.image()

# end of script









