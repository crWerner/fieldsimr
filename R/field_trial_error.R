#' Simulate plot-level errors for a plant breeding trial
#'
#' Creates a data frame with simulated plot-level errors for one or more traits in plant breeding
#' trials across multiple environments. The simulated error consists of a spatial error term and
#' a random error term. The spatial error term can be simulated based on 1) bivariate
#' interpolation using the \link[interp]{interp} function of the package 'interp', or 2) a
#' separable first-order autoregressive process (AR1:AR1). The random error term is simulated
#' using an independent process. The spatial and random error terms are combined according to a
#' user-defined ratio. \cr
#' For multiple traits, correlated error terms can be generated assuming 1) correlated spatial
#' errors between traits, 2) correlated random error between traits, or 3) a combination of both. \cr
#' A separable covariance structure is assumed between traits and environments.
#'
#' @param n_envs Number of environments to be simulated (same as for \code{compsym_asr_input}
#'   or \code{unstr_asr_output}, where applicable).
#' @param n_traits Number of traits to be simulated.
#' @param n_cols A vector defining the total number of columns in each environment. If only one
#'   value is provided and \code{n_traits > 1}, all environments will be assigned the same number
#'   of columns.
#' @param n_rows A vector defining the total number of rows in each environment. If only one
#'   value is provided and \code{n_traits > 1}, all environments will be assigned the same number
#'   of rows.
#' @param plot_length A vector defining the plot length in each environment. If only one value is
#'   provided and \code{n_traits > 1}, the plots in all environments will be assigned the same
#'   plot length.
#' @param plot_width A vector defining the plot width in each environment. If only one value is
#'   provided and \code{n_traits > 1}, the plots in all environments will be assigned the same
#'   plot width.
#' @param n_reps A vector defining the number of complete replicates in each environment. If only
#'   one value is provided and \code{n_traits > 1}, all environments will be assigned the same
#'   number of replicates.
#' @param rep_dir A character string specifying the direction of replicate blocks. One of either
#'   "column" (side-by-side, the default) or "row" (above-and-below). \code{rep_dir} is ignored
#'   when \code{n_reps = 1}.
#' @param var_R A vector of error variances for each trait-by-environment combination (ordered
#'   as environments within traits). If the length of \code{var_R} is equal to \code{n_traits},
#'   all traits will be assigned the same error variance in each environment.
#' @param S_cor_R A matrix of spatial error correlations between more than one trait. If not
#'   defined and \code{n_traits > 1}, a diagonal matrix is constructed.
#' @param R_cor_R A matrix of random error correlations between more than one trait. If not
#'   defined and \code{n_traits > 1}, a diagonal matrix is constructed.
#' @param spatial_model A character string specifying the model used to simulate the spatial error
#'   term. One of either "Bivariate" (bivariate interpolation, the default) or "AR1:AR1"
#'   (separable first-order autoregressive process (AR1:AR1)).
#' @param prop_spatial A vector defining the proportion of spatial error variance to total
#'   error variance (spatial + random). If only one value is provided, all environments will be
#'   assigned the same proportion. By default, the spatial error variance accounts for half of the
#'   total error variance (\code{prop_spatial = 0.5}).
#' @param complexity A scalar defining the complexity of the bivariate interpolation model.
#'   By default, \code{complexity = 12}. Note that low values may lead to convergence problems.
#'   See \link[interp]{interp} for further details.
#' @param col_cor A vector of column autocorrelations for each environment used in the AR1:AR1
#'   spatial error model. If only one value is provided, all environments will be assigned the
#'   same column autocorrelation.
#' @param row_cor A vector of row autocorrelations for each environment used in the AR1:AR1
#'   spatial error model. If only one value is provided, all environments will be assigned the
#'   same row autocorrelation.
#' @param return_effects When true, a list is returned with additional entries for each trait
#'   containing the spatial and random errors.
#'
#' @return A data-frame containing the environment, block, column and row names, as well as the
#'   simulated error for each trait. When \code{return_effects = TRUE}, a list is returned with
#'   additional entries for each trait containing the spatial and random errors.
#'
#' @examples
#' # Simulation of plot-level errors for two traits in three environments using a bivariate
#' # interpolation model for spatial variation.
#'
#' n_envs <- 3 # Number of simulated environments.
#' n_traits <- 2 # Number of simulated traits.
#'
#' # Field layout
#' n_cols <- 10 # Total number of columns in each environment.
#' n_rows <- c(20, 30, 20) # Total number of rows in each environment.
#' plot_length <- 5 # Plot length set to 5 meters in each environment.
#' plot_width <- 2 # Plot width set to 2 meters in each environment.
#' n_reps <- c(2, 3, 2) # Number of complete replicates (blocks) per environment.
#'
#' # Error variances for traits 1 and 2.
#' var_R <- c(0.4, 15)
#'
#' # Spatial error correlations between traits 1 and 2.
#' S_cor_R <- matrix(c(
#'   1.0, 0.2,
#'   0.2, 1.0
#' ),
#' ncol = 2
#' )
#'
#' error_df <- field_trial_error(
#'   n_envs = n_envs,
#'   n_traits = n_traits,
#'   n_cols = n_cols,
#'   n_rows = n_rows,
#'   plot_length = plot_length,
#'   plot_width = plot_width,
#'   n_reps = n_reps,
#'   rep_dir = "row",
#'   var_R = var_R,
#'   S_cor_R = S_cor_R,
#'   spatial_model = "bivariate",
#'   prop_spatial = 0.6,
#'   complexity = 14,
#'   return_effects = TRUE
#' )
#' @export
field_trial_error <- function(n_envs,
                              n_traits,
                              n_reps,
                              n_cols,
                              n_rows,
                              plot_length,
                              plot_width,
                              rep_dir = "column",
                              var_R,
                              S_cor_R = NULL,
                              R_cor_R = NULL,
                              spatial_model = "bivariate",
                              prop_spatial = 0.5,
                              complexity = 12,
                              col_cor = NULL,
                              row_cor = NULL,
                              return_effects = FALSE) {
  if (n_envs < 1 | n_envs %% 1 != 0) stop("'n_envs' must be an integer > 0")
  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  if (min(n_cols) < 1 | any(n_cols %% 1 != 0)) {
    stop("'n_cols' must contain integers > 0")
  }
  if (length(n_cols) == 1) n_cols <- rep(n_cols, n_envs)
  if (length(n_cols) != n_envs) {
    stop("Length of vector 'n_cols' does not match total number of environments")
  }
  if (min(n_rows) < 1 | any(n_rows %% 1 != 0)) {
    stop("'n_rows' must contain integers > 0")
  }
  if (length(n_rows) == 1) n_rows <- rep(n_rows, n_envs)
  if (length(n_rows) != n_envs) {
    stop("Length of vector 'n_rows' does not match total number of environments")
  }

  if (min(plot_length) <= 0) stop("'plot_length' must contain values > 0")
  if (length(plot_length) == 1) plot_length <- rep(plot_length, n_envs)
  if (length(plot_length) != n_envs) {
    stop("Length of vector 'plot_length' does not match total number of environments")
  }
  if (min(plot_width) <= 0) stop("'plot_width' must be > 0")
  if (length(plot_width) == 1) plot_width <- rep(plot_width, n_envs)
  if (length(plot_width) != n_envs) {
    stop("Length of vector 'plot_width' does not match total number of environments")
  }

  if (min(n_reps) < 1 | any(n_reps %% 1 != 0)) {
    stop("'n_reps' must contain integers > 0")
  }

  rep_dir <- tolower(rep_dir)
  if (rep_dir == "column") {
    if (any((n_cols / n_reps) %% 1 != 0)) {
      stop("Number of columns not divisible by number of reps in at least one
         environment. Review your trial design!")
    }
  } else if (rep_dir == "row") {
    if (any((n_rows / n_reps) %% 1 != 0)) {
      stop("Number of rows not divisible by number of reps in at least one
           environment. Review your trial design!")
    }
  } else {
    stop("'rep_dir' must be 'row' or 'column'")
  }

  if (length(n_reps) == 1) n_reps <- rep(n_reps, n_envs)
  if (length(n_reps) != n_envs) {
    stop("Length of vector 'n_reps' does not match total number of environments")
  }

  if (length(var_R) == 1) var_R <- rep(var_R, n_traits)
  if (length(var_R) == n_traits) var_R <- rep(var_R, each = n_envs)
  if (any(var_R <= 0)) {
    stop("'var_R' must contain values greater than 0")
  }
  if (length(var_R) != (n_envs * n_traits)) {
    stop("Number of values in argument 'var_R' must either match number of traits
         or number of trait x environment combinations")
  }

  if (is.null(S_cor_R)) S_cor_R <- diag(n_traits)
  if (is.null(R_cor_R)) R_cor_R <- diag(n_traits)

  spatial_model <- tolower(spatial_model)
  if (spatial_model != "bivariate" & spatial_model != "ar1:ar1") {
    stop("'spatial_model' must be 'bivariate' or 'AR1:AR1'")
  }

  if (any(prop_spatial < 0) | any(prop_spatial > 1)) {
    stop("'prop_spatial' must contain values between 0 and 1")
  }
  if (length(prop_spatial) == 1) prop_spatial <- rep(prop_spatial, n_envs)
  if (length(prop_spatial) != n_envs) {
    stop("Length of vector 'prop_spatial' does not match total number of environments")
  }

  envs <- rep(1:n_envs, times = n_cols * n_rows)
  reps <- c(unlist(mapply(function(x, y, z) rep(1:z, each = c(x * y / z)), x = n_cols, y = n_rows, z = n_reps)))
  cols <- c(unlist(mapply(function(x, y) rep(1:x, each = y), x = n_cols, y = n_rows)))
  rows <- c(unlist(mapply(function(x, y) rep(1:x, times = y), y = n_cols, x = n_rows)))

  if (rep_dir == "row") {
    cols <- c(unlist(mapply(function(x, y) rep(1:x, times = y), x = n_cols, y = n_rows)))
    rows <- c(unlist(mapply(function(x, y) rep(1:x, each = y), y = n_cols, x = n_rows)))
  }

  plot_df <- data.frame(
    env = envs,
    block = reps,
    col = cols,
    row = rows
  )

  plot_df <- plot_df[order(plot_df$env, plot_df$col, plot_df$row), ]
  rownames(plot_df) <- NULL

  if (spatial_model == "ar1:ar1") {
    if (any(col_cor < 0) | any(col_cor > 1)) {
      stop("'col_cor' must contain values between 0 and 1'")
    }
    if (length(col_cor) == 1) col_cor <- rep(col_cor, n_envs)
    if (length(col_cor) != n_envs) {
      stop("Length of vector 'col_cor' does not match total number of environments")
    }
    if (any(row_cor < 0) | any(row_cor > 1)) stop("row_cor must be between 0 and 1")
    if (length(row_cor) == 1) row_cor <- rep(row_cor, n_envs)
    if (length(row_cor) != n_envs) {
      stop("Length of vector 'row_cor' does not match total number of environments")
    }

    power_lst <- lapply(n_cols, function(x) abs(outer(1:x, 1:x, "-")))
    col_ar1 <- mapply(function(x, y) x^y, x = col_cor, y = power_lst, SIMPLIFY = FALSE)

    power_lst <- lapply(n_rows, function(x) abs(outer(1:x, 1:x, "-")))
    row_ar1 <- mapply(function(x, y) x^y, x = row_cor, y = power_lst, SIMPLIFY = FALSE)

    cor_mat_lst <- mapply(function(x, y) kronecker(x, y), x = col_ar1, y = row_ar1, SIMPLIFY = FALSE)
    cor_mat_lst <- mapply(function(x) kronecker(S_cor_R, x), x = cor_mat_lst, SIMPLIFY = FALSE)

    l_lst <- lapply(cor_mat_lst, function(x) chol(x))
    plot_error_lst1 <- mapply(function(x, y) matrix(c(stats::rnorm(x) %*% y), ncol = n_traits),
                              x = n_cols * n_rows * n_traits, y = l_lst, SIMPLIFY = FALSE
    )
  }

  if (spatial_model == "bivariate") {
    if (complexity <= 0) stop("'complexity' must be an integer > 0")

    cols_lst <- with(plot_df, tapply(col, env, function(x) c(unique(x), max(x) + 1)))
    col_centres_lst <- mapply(function(x, y) round(y * (x - 0.5),8),
                              x = cols_lst, y = plot_length, SIMPLIFY = FALSE
    )

    rows_lst <- with(plot_df, tapply(row, env, function(x) c(unique(x), max(x) + 1)))
    row_centres_lst <- mapply(function(x, y) round(y * (x - 0.5),8),
                              x = rows_lst, y = plot_width, SIMPLIFY = FALSE
    )

    col_gap <- plot_length / 4
    row_gap <- plot_width / 4

    plot_error_lst1 <- NA
    while(sum(is.na(unlist(plot_error_lst1))) > 0){
    xInterp_list <- mapply(function(x, y, z) c(0 - z, x * y + z, 0 - z, x * y + z, sample(stats::runif(n = complexity, min = 0, max = (x * y)))),
                           x = n_cols, y = plot_length, z = col_gap, SIMPLIFY = FALSE
    )
    yInterp_list <- mapply(function(x, y, z) c(0 - z, 0 - z, x * y + z, x * y + z, sample(stats::runif(n = complexity, min = 0, max = (x * y)))),
                           x = n_rows, y = plot_width, z = row_gap, SIMPLIFY = FALSE
    )
    zInterp_list <- lapply(n_cols, function(x) scale(matrix(stats::rnorm((4 + complexity) * n_traits), ncol = n_traits)) %*% chol(S_cor_R))

    for (i in 1:n_traits) {
      tmp <- mapply(function(v, w, x, y, z, xo, yo) {
        c(t(interp::interp(x = x, y = y, z = z[, i], xo = c(xo), yo = c(yo), linear = F, extrap = T, duplicate = "mean")$z)[1:v, 1:w])
      },
      v = n_rows, w = n_cols, x = xInterp_list, y = yInterp_list, z = zInterp_list, xo = col_centres_lst, yo = row_centres_lst, SIMPLIFY = FALSE
      )
      # }
      if (i == 1) {
        plot_error_lst1 <- tmp
      }
      if (i > 1) {
        plot_error_lst1 <- Map("cbind", plot_error_lst1, tmp)
      }
     }
    }
  }

  plot_error_lst2 <- mapply(function(x) matrix(c(stats::rnorm(x)), ncol = n_traits) %*% chol(R_cor_R),
                            x = n_cols * n_rows * n_traits, SIMPLIFY = FALSE
  )

  var_R <- as.data.frame(t(matrix(var_R, ncol = n_traits, byrow = TRUE)))

  var_R <- lapply(X = var_R, FUN = c)
  e_spat <- mapply(function(w, x) (scale(w) * sqrt(x)), w = plot_error_lst1, x = prop_spatial, SIMPLIFY = F)
  e_rand <- mapply(function(x, y) (scale(y) * sqrt(1 - x)), x = prop_spatial, y = plot_error_lst2, SIMPLIFY = F)
  if (n_traits > 1) {
    e_scale <- mapply(function(x, y) sqrt(diag(1 / diag(as.matrix(stats::var(x + y))))), x = e_spat, y = e_rand, SIMPLIFY = F)
    e_spat <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_spat, y = e_scale, z = var_R, SIMPLIFY = F)
    e_rand <- mapply(function(x, y, z) x %*% y %*% diag(sqrt(z)), x = e_rand, y = e_scale, z = var_R, SIMPLIFY = F)
  }
  if (n_traits == 1) {
    e_scale <- mapply(function(x, y) sqrt(1 / stats::var(x + y)), x = e_spat, y = e_rand, SIMPLIFY = F)
    e_spat <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_spat, y = e_scale, z = var_R, SIMPLIFY = F)
    e_rand <- mapply(function(x, y, z) x %*% y %*% sqrt(z), x = e_rand, y = e_scale, z = var_R, SIMPLIFY = F)
  }
  plot_error_lst <- mapply(function(x, y) x + y, x = e_spat, y = e_rand, SIMPLIFY = F)
  plot_error <- do.call(what = "rbind", plot_error_lst)
  colnames(plot_error) <- paste0("e.Trait.", 1:n_traits)
  plot_df <- cbind(plot_df, plot_error)

  if (return_effects) {
    e_spat <- do.call("rbind", e_spat)
    e_rand <- do.call("rbind", e_rand)
    e_all <- lapply(seq_len(ncol(e_spat)), function(i) cbind(e_spat[, i], e_rand[, i]))
    resids <- lapply(e_all, function(x) {
      data.frame(plot_df[, 1:4],
                 e_spatial = x[, 1],
                 e_random = x[, 2]
      )
    })

    list_names <- c("plot_df", paste0("Trait.", 1:n_traits))
    plot_df <- list(plot_df)
    plot_df <- c(plot_df, resids)
    names(plot_df) <- list_names
  }
  plot_df
  return(plot_df)
}
