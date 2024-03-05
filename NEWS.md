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

* Added function `qq_plot` to create quantile-quantile (Q-Q) plots.

* Added function `sample_variogram` to create sample variograms.

* Added function `theoretical_variogram` to create theoretical variograms.

# FieldSimR 1.3.1

* Replaced all instances of `_` with `.` in all function arguments, e.g., `pos.def` replaced `pos_def`.

* Replaced `n_` with `n` in all function arguments, e.g., `ntraits` replaced `n_traits`.

* Added `multi_asr_input` and `multi_asr_output` wrapper functions for simulating genetic values based on a multiplicative model for GxE interaction.

* Added `small.positive` argument to function `rand_cor_mat`, which is passed to the `bend` function.

* Changed names of example data frames to `error_df_bivar` and `gv_df_unstr` from `df_error_bivar` and `df_gv_unstr`.

* Replaced `rel.main.eff.A` with `prop.main`, `rel.main.eff.DD` with `prop.mainD`, and `rel.main.eff.AA` with `prop.mainAA` in the `compsym_asr_input` function, since these arguments define the proportion of main effect variance, not the relative magnitude. 

* `prop.main` was implemented instead of `prop.mainA`, since this argument is aligned with `var` - i.e., it represents the proportion of additive or total main effect variance depending on whether `useVarA = TRUE` or `FALSE` in `AlphaSimR`.

* Added `return.effects` argument to the `make_phenotypes` function for returning the plot errors and genetic values. The latter will be returned in randomised order when `randomise = TRUE`.

* Added `plot_matrix` function for graphically displaying a symmetric matrix, e.g., correlation or covariance matrix.






