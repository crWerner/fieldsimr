# FieldSimR 1.1.0

* Added a `NEWS.md` file to track changes to the package.

* Added argument `pos_def` to function `rand_cor_mat` to make bending of a non-positive-definite correlation matrix to a positive-definite matrix optional.

* Added example data frame `df_error_bivar` with plot errors for two traits across three locations generated using `field_trial_error`.

* Added example data frame `df_gv_unstr ` with simulated genetic values for two traits across three environments generated using `unstr_asr_input` and `unstr_asr_output`.

* Added functionality to simulate extraneous variation to `field_trial_error`.

* Added internal functions `spline_interp` and `fill_matrix` to interpolate and extrapolate missing values if some NAs remain after bivariate interpolation.

* Added `make_phenotypes` to create phenotypes through combination of genetic values and plot errors.

* Added `qq_plot` to compare the theoretical quantiles of a normal distribution with the sample quantiles of the distribution of a user-defined effect.

* Added `sample_variogram` to create a variogram of a user-defined effect.

* Added `theoretical_variogram` to create a theoretical variogram.

* Added vignette `compound_symmetry_GxE_demo` to demonstrate the simulation of genetic values using a compound symmetry GxE model.

* Added vignette `spatial_error_demo` to demonstrate the simulation of plot errors and phenotypes for a multi-environment plant breeding trial.

* Added vignette `unstructured_GxE_demo` to demonstrate the simulation of genetic values using an unstructured GxE model.
 
* Removed argument `env` from function `plot_effects`.

* Replaced package `fields` for graphics in `plot_effects` by `ggplot2`.

* Set the `complexity` argument of `field_trial_error` by default to the maximum number of columns and rows in each environment.

* Updated Description in `DESCRIPTION`.


# FieldSimR 1.2.1

* Argument `ext_ord` replaced arguments `ext_col_cor` and `ext_row_cor` in function `field_trial_error`.

* Factorised argument `env`, `rep` and `id` in functions `field_trial_error`, `make_phenotypes`, `unstr_asr_output`, `compsym_asr_output`.

* Randomisation fixed in function `make_phenotypes`.

* Changed default parameters for `col_cor` and `row_cor`, `prop_spatial`, and `complexity` in function `field_trial_error`.

* Argument `plot_labels` added to function `plot_effects`.


