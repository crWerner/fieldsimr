### TRAIT VALUES: means and variances for grain yield (fresh weight) in tons and plant height in cm.

gy_fwt_mean <- 5.2 # tons/ha
gy_fwt_var <- 0.075

pht_mean <- 230 # cm
pht_var <- 13

rnorm_gy <- rnorm(10000, gy_fwt_mean, sqrt(gy_fwt_var))
hist(rnorm_gy, breaks = 100)
range(rnorm_gy)

rnorm_pht <- rnorm(10000, pht_mean, sqrt(pht_var))
hist(rnorm_pht, breaks = 100)
range(rnorm_pht)


### Unstructured model

n_env <- 3

gy_fwt_env_means <- rnorm(n_env, gy_fwt_mean, 0.15)
gy_fwt_env_means
pht_env_means <- rnorm(n_env, pht_mean, 20)
pht_env_means

gy_fwt_env_var <- rnorm(n_env, gy_fwt_var, 0.05)
gy_fwt_env_var
pht_env_var <- rnorm(n_env, pht_var, 5)
pht_env_var
