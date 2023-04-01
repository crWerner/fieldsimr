#' Simulate genetic values based on a compound symmetry model for GxE interaction - `AlphaSimR' input
#' parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR'} to simulate
#' genetic values for multiple traits across multiple environments based on a compound symmetry model
#' for genotype-by-environment (GxE) interaction. \cr
#' By default, `AlphaSimR' does not support complex models for GxE interaction. However, its
#' functionality to simulate correlated genetic values can be utilised for this purpose by
#' providing the required variance structures. \code{compsym_asr_input} is a wrapper function to
#' construct the variance structures required to simulate GxE interaction in `AlphaSimR' based on
#' a multi-trait compound symmetry model. This function assumes a separable structure between traits and
#' environments. After simulating the genetic values, the wrapper function
#' \link[FieldSimR]{compsym_asr_output} can be used to obtain data frames with the values.
#'
#' \strong{Note:} `AlphaSimR' can simulate different biological effects (see:
#' \code{\link[AlphaSimR]{SimParam}}).
#' \itemize{
#'   \item For additive traits use \code{addTraitA()}.
#'   \item For additive + dominance traits use \code{addTraitAD()}.
#'   \item For additive + epistatic traits use \code{addTraitAE()}.
#'   \item For additive + dominance + epistatic traits use \code{addTraitADE()}.
#'   }
#' If non-additive effects are to be simulated, check the \code{useVarA} argument of these
#' functions.
#'
#' @param n_envs Number of environments to be simulated. A minimum of two environments is required. 
#'   By default, \code{n_envs = 3}.
#' @param n_traits Number of traits to be simulated. By default, \code{n_traits = 2}.
#' @param mean A vector of mean genetic values for each environment-within-trait combination.
#'   If only one value is provided, all environment-within-trait combinations will be assigned the same mean.
#'   By default, \code{mean = 0}.
#' @param var A vector of genetic variances for each trait. Simulated traits are restricted by the
#'   compound symmetry model to having the same variance for each environment (i.e., main
#'   effect variance + GxE interaction variance) and the same covariance between each pair of
#'   environments (main effect variance). By default, \code{var = 1}. \cr
#'   \strong{Note:} When \code{useVarA = TRUE} is specified in `AlphaSimR' (default) the values in
#'   \code{var} represent the \code{additive} genetic variances, otherwise they will represent the
#'   \code{total} (additive + non-additive) genetic variances.
#' @param rel_main_eff_A  A vector defining the magnitude of the additive main effect variance
#'   relative to the additive main effect + GxE interaction variance for each trait. If only one
#'   value is provided, all traits will be assigned the same value. By default, \code{rel_main_eff_A = 0.5}. \cr
#'   \strong{Note:} \code{0 < rel_main_eff_A < 1}.
#' @param cor_A A matrix of additive genetic correlations between traits. If not
#'   defined, a diagonal matrix is constructed.
#' @param mean_DD A vector of mean dominance degrees for each environment-within-trait combination
#'   (similar to \code{mean}). If only one value is provided, all environment-within-trait combinations 
#'   will be assigned the same mean. By default, \code{mean_DD = NULL} and dominance is not simulated.
#' @param var_DD A vector of dominance degree variances for each trait. Simulated traits have the
#'   same dominance degree variance for each environment and the same dominance degree covariance
#'   between each pair of environments (similar to \code{var}). By default, \code{var_DD = NULL}.
#' @param rel_main_eff_DD A vector defining the magnitude of the dominance degree main effect
#'   variance relative to the main effect + GxE interaction variance for each trait (similar to
#'   \code{rel_main_eff_A}). By default, \code{rel_main_eff_DD = NULL}. \cr
#'   \strong{Note:} \code{0 < rel_main_eff_DD < 1}.
#' @param cor_DD A matrix of dominance degree correlations between traits (similar
#'   to \code{cor_A}). If not defined and dominance is simulated, a diagonal matrix is constructed. 
#'   By default, \code{cor_DD = NULL}.
#' @param rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each trait, that is in a diploid organism with
#'   allele frequency 0.5. Simulated traits have the same epistatic variance for each environment
#'   and the same epistatic covariance between each pair of environments (similar to \code{var}).
#'   By default, \code{rel_AA = NULL} and epistasis is not simulated.
#' @param rel_main_eff_AA A vector defining the magnitude of the epistatic main effect variance
#'   relative to the main effect + GxE interaction variance for each trait (similar to
#'   \code{rel_main_eff_A}). By default, \code{rel_main_eff_AA = NULL}. \cr
#'   \strong{Note:} \code{0 < rel_main_eff_AA < 1}.
#' @param cor_AA A matrix of epistatic correlations between traits (similar to
#'   \code{cor_A}). If not defined and epistasis is simulated, a diagonal matrix is constructed. By
#'   default, \code{cor_AA = NULL}.
#'
#' @return A list containing input parameters for `AlphaSimR', which is used to simulate
#'   correlated genetic effects based on a compound symmetry model.
#'
#' @examples
#' # Simulate genetic values in 'AlphaSimR' for two additive + dominance traits across
#' # three environments based on a compound symmetry model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(1, 3, 2, 80, 70, 100) # Trait 1 x 3 environments, trait 2 x 3 environments.
#' mean_DD <- c(0.1, 0.4) # Trait 1 and 2, same values set in all three environments.
#'
#' # Additive genetic variances (useVarA = TRUE) and dominance degree variances.
#' var <- c(0.2, 10) # Different values set for traits 1 and 2.
#' var_DD <- c(0.1, 0.2) # Different values set for traits 1 and 2.
#'
#' # Relative magnitude of the additive and dominance degree main effect variances.
#' rel_main_eff_A <- c(0.4, 0.6) # Different values set for traits 1 and 2.
#' rel_main_eff_DD <- 0.8 # Same value set for traits 1 and 2.
#'
#' # Additive and dominance degree correlations between the two simulated traits.
#' cor_A <- matrix(c(1.0, 0.3, 0.3, 1.0), ncol = 2) # Additive correlation matrix.
#' cor_DD <- diag(2) # Dominance correlation matrix - assume independence.
#'
#' input_asr <- compsym_asr_input(
#'   n_envs = 3,
#'   n_traits = 2,
#'   mean = mean,
#'   var = var,
#'   rel_main_eff_A = rel_main_eff_A,
#'   cor_A = cor_A,
#'   mean_DD = mean_DD,
#'   var_DD = var_DD,
#'   rel_main_eff_DD = rel_main_eff_DD,
#'   cor_DD = cor_DD
#' )
#' @export
compsym_asr_input <- function(n_envs = 3,
                              n_traits = 2,
                              mean = 0,
                              var = 1,
                              rel_main_eff_A = 0.5,
                              cor_A = NULL,
                              mean_DD = NULL,
                              var_DD = NULL,
                              rel_main_eff_DD = NULL,
                              cor_DD = NULL,
                              rel_AA = NULL,
                              rel_main_eff_AA = NULL,
                              cor_AA = NULL) {
  if (n_envs < 2) stop("'n_envs' must be > 1")
  if (n_envs %% 1 != 0) stop("'n_envs' must be an integer")
  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  if (is.null(mean_DD) & is.null(var_DD) & is.null(rel_AA)) {
    labels <- "A"
  } else if (is.null(mean_DD) & is.null(var_DD)) {
    labels <- c("A", "AA")
  } else if (is.null(rel_AA)) {
    labels <- c("A", "DD")
  } else {
    labels <- c("A", "DD", "AA")
  }

  for (i in labels) {
    if (i == "A") {
      if (length(mean) == n_traits) {
        mean_vals <- rep(mean, each = n_envs)
      } else if (length(mean) == (n_traits * n_envs)) {
        mean_vals <- mean
      } else {
        stop("Number of values in argument 'mean' must either match number of
               traits or number of trait x environment combinations")
      }

      if (length(var) != n_traits) {
        stop("Number of values in argument 'var' must match number of traits")
      }

      if (length(rel_main_eff_A) == n_traits) {
        rel_main_eff <- rel_main_eff_A
      } else if (length(rel_main_eff_A) == 1 & n_traits > 1) {
        rel_main_eff <- rep(rel_main_eff_A, n_traits)
      } else {
        stop("Number of values in 'rel_main_eff_A' has must either be 1 or must
                 match number of traits")
      }

      if (any(rel_main_eff_A <= 0) | any(rel_main_eff_A >= 1)) {
        stop("'rel_main_eff_A' must contain values between 0 and 1")
      }

      if (is.null(cor_A)) cor_A <- diag(n_traits)

      if (nrow(cor_A) != n_traits | ncol(cor_A) != n_traits) {
        stop("Dimension of 'cor_A' does not match number of traits")
      }

      if (any(unique(diag(cor_A)) != 1) | any(cor_A > 1) | any(cor_A < (-1))) {
        stop("'cor_A' is not a correlation matrix")
      }

      if (!isSymmetric(cor_A)) stop("'cor_A' is not symmetric")

      main_mean <- colMeans(matrix(mean_vals, ncol = n_traits))
      vars <- var
      T_cor <- as.matrix(cor_A)
    }

    if (i == "DD") {
      if (length(mean_DD) == n_traits) {
        mean_vals <- rep(mean_DD, each = n_envs)
      } else if (length(mean_DD) == (n_traits * n_envs)) {
        mean_vals <- mean_DD
      } else {
        stop("Number of values in argument 'mean_DD' must either match number
            of traits or number of trait x environment combinations")
      }

      if (length(var_DD) != n_traits) {
        stop("Number of values in argument 'var_DD' must match number of traits")
      }

      if (is.null(rel_main_eff_DD)) stop("'rel_main_eff_DD' is not defined")

      if (length(rel_main_eff_DD) == n_traits) {
        rel_main_eff <- rel_main_eff_DD
      } else if (length(rel_main_eff_DD) == 1 & n_traits > 1) {
        rel_main_eff <- rep(rel_main_eff_DD, n_traits)
      } else {
        stop("Number of values in 'rel_main_eff_DD' has must either be 1 or must
               match number of traits")
      }

      if (any(rel_main_eff_DD <= 0) | any(rel_main_eff_DD >= 1)) {
        stop("'rel_main_eff_DD' must contain values between 0 and 1")
      }

      if (is.null(cor_DD)) cor_DD <- diag(n_traits)

      if (nrow(cor_DD) != n_traits | ncol(cor_DD) != n_traits) {
        stop("Dimension of 'cor_DD' does not match number of traits")
      }

      if (any(unique(diag(cor_A)) != 1) | any(cor_DD > 1) | any(cor_DD < (-1))) {
        stop("'cor_DD' is not a correlation matrix")
      }

      if (!isSymmetric(cor_DD)) stop("'cor_DD' is not symmetric")

      main_mean <- colMeans(matrix(mean_vals, ncol = n_traits))
      vars <- var_DD
      T_cor <- as.matrix(cor_DD)
    }

    if (i == "AA") {
      if (length(rel_AA) != n_traits) {
        stop("Number of values in argument 'rel_AA' must match number of traits")
      }

      if (is.null(rel_main_eff_AA)) stop("'rel_main_eff_AA' is not defined")

      if (length(rel_main_eff_AA) == n_traits) {
        rel_main_eff <- rel_main_eff_AA
      } else if (length(rel_main_eff_AA) == 1 & n_traits > 1) {
        rel_main_eff <- rep(rel_main_eff_AA, n_traits)
      } else {
        stop("Number of values in 'rel_main_eff_AA' has must either be 1 or must
               match number of traits")
      }

      if (any(rel_main_eff_AA <= 0) | any(rel_main_eff_AA >= 1)) {
        stop("'rel_main_eff_AA' must contain values between 0 and 1")
      }

      if (is.null(cor_AA)) cor_AA <- diag(n_traits)

      if (nrow(cor_AA) != n_traits | ncol(cor_AA) != n_traits) {
        stop("Dimension of 'cor_AA' does not match number of traits")
      }

      if (any(unique(diag(cor_A)) != 1) | any(cor_AA > 1) | any(cor_AA < (-1))) {
        stop("'cor_AA' is not a correlation matrix")
      }

      if (!isSymmetric(cor_AA)) stop("'cor_AA' is not symmetric")


      mean_vals <- rep(0, n_traits * n_envs)
      main_mean <- rep(0, n_traits)
      vars <- rel_AA
      T_cor <- as.matrix(cor_AA)
    }

    gxe_devs <- as.list(as.data.frame(matrix(mean_vals, ncol = n_traits) - rep(main_mean, each = n_envs)))
    main_var <- vars * rel_main_eff
    gxe_var <- vars - main_var

    mean_pseudo <- c(mapply(function(x, y) c(x, y), x = main_mean, y = gxe_devs))
    var_pseudo <- c(mapply(function(x, y) c(x, rep(y, n_envs)), x = main_var, y = gxe_var))

    cor_pseudo <- kronecker(T_cor, (diag(n_envs + 1)))

    if (i == "A") {
      input_asr <- list(
        mean = mean_pseudo,
        var = var_pseudo,
        cor_A = cor_pseudo
      )
    }

    if (i == "DD") {
      input_asr <- c(
        input_asr,
        list(
          mean_DD = mean_pseudo,
          var_DD = var_pseudo,
          cor_DD = cor_pseudo
        )
      )
    }

    if (i == "AA") {
      input_asr <- c(
        input_asr,
        list(
          mean_AA = mean_pseudo,
          rel_AA = var_pseudo,
          cor_AA = cor_pseudo
        )
      )
    }
  }

  return(input_asr)
}

#' Simulate genetic values based on a compound symmetry model for GxE interaction -
#' Simulation using `AlphaSimR'
#'
#' Creates a data frame of simulated genetic values for multiple traits across multiple environments
#' based on a compound symmetry model for genotype-by-environment (GxE) interaction. This function
#' requires an `AlphaSimR' population object generated using the function
#' \link[FieldSimR]{compsym_asr_input}.
#'
#' @param pop An `AlphaSimR' population object (\code{\link[AlphaSimR]{Pop-class}} or
#'   \code{\link[AlphaSimR]{HybridPop-class}}) generated using \link[FieldSimR]{compsym_asr_input}.
#' @param n_envs Number of simulated environments (same number used in
#'   \link[FieldSimR]{compsym_asr_input}).
#' @param n_traits Number of simulated traits (same number used in \link[FieldSimR]{compsym_asr_input}).
#' @param n_reps A vector defining the number of complete replicates in each environment. If only
#'   one value is provided, all environments will be assigned the same number.
#' @param effects When TRUE, a list is returned with additional entries containing the total
#'   (additive + dominance + epistatic) main effects and GxE interaction effects for each
#'   environment-within-trait combination. By default, effects = FALSE.
#'
#' @return A data frame with columns `env', `rep', genotype `id', and the
#'   simulated genetic values for each trait. When \code{effects = TRUE}, a list is returned with
#'   additional entries containing the total (additive + dominance + epistatic) main effects and
#'   GxE interaction effects for each environment-within-trait combination.
#'
#' @examples
#' # Simulate genetic values in 'AlphaSimR' for two additive + dominance traits across
#' # three environments based on a compound symmetry model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(1, 3, 2, 80, 70, 100) # Trait 1 x 3 environments, trait 2 x 3 environments.
#' mean_DD <- c(0.1, 0.4) # Trait 1 and 2, same values set in all three environments.
#'
#' # Additive genetic variances (useVarA = TRUE) and dominance degree variances.
#' var <- c(0.2, 10) # Different values set for traits 1 and 2.
#' var_DD <- c(0.1, 0.2) # Different values set for traits 1 and 2.
#'
#' # Relative magnitude of additive and dominance degree main effect variances.
#' rel_main_eff_A <- c(0.4, 0.6) # Different values set for traits 1 and 2.
#' rel_main_eff_DD <- 0.8 # Same value set for traits 1 and 2.
#'
#' # Additive and dominance degree correlations between the two simulated traits.
#' cor_A <- matrix(c(1.0, 0.3, 0.3, 1.0), ncol = 2) # Additive correlation matrix.
#' cor_DD <- diag(2) # Dominance correlation matrix - assume independence.
#'
#' input_asr <- compsym_asr_input(
#'   n_envs = 3,
#'   n_traits = 2,
#'   mean = mean,
#'   var = var,
#'   rel_main_eff_A = rel_main_eff_A,
#'   cor_A = cor_A,
#'   mean_DD = mean_DD,
#'   var_DD = var_DD,
#'   rel_main_eff_DD = rel_main_eff_DD,
#'   cor_DD = cor_DD
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on a compound symmetry model
#' # for GxE interaction.
#'
#' library("AlphaSimR")
#' FOUNDERPOP <- quickHaplo(
#'   nInd = 100,
#'   nChr = 6,
#'   segSites = 100
#' )
#'
#' SP <- SimParam$new(FOUNDERPOP)
#'
#' SP$addTraitAD(
#'   nQtlPerChr = 100,
#'   mean = input_asr$mean,
#'   var = input_asr$var,
#'   meanDD = input_asr$mean_DD,
#'   varDD = input_asr$var_DD,
#'   corA = input_asr$cor_A,
#'   corDD = input_asr$cor_DD,
#'   useVarA = TRUE
#' )
#'
#' # By default, the value provided in 'var' represents the additive variance.
#' # If useVarA=FALSE, 'var' represents the total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for the two traits across the
#' # three environments.
#'
#' n_reps <- c(2, 3, 2) # Vector containing the number of complete replicates in each environment.
#'
#' gv_df <- compsym_asr_output(
#'   pop = pop,
#'   n_envs = 3,
#'   n_traits = 2,
#'   n_reps = n_reps,
#'   effects = TRUE
#' )
#' @export
compsym_asr_output <- function(pop,
                               n_envs,
                               n_traits,
                               n_reps,
                               effects = FALSE) {
  if (n_envs < 2) stop("'n_envs' must be > 1")
  if (n_envs %% 1 != 0) stop("'n_envs' must be an integer")

  if ((sum(n_reps < 1) > 0) | (sum(n_reps %% 1 != 0) > 0)) {
    stop("'n_reps' must contain positive integer values")
  }
  if (length(n_reps) == 1) n_reps <- rep(n_reps, n_envs)

  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  envs <- rep(1:n_envs, times = length(pop@id) * n_reps)
  reps <- unlist(lapply(n_reps, function(x) rep(1:x, each = length(pop@id))))

  g_main <- as.list(as.data.frame(pop@gv[, seq(1, (n_traits + n_traits * n_envs), (n_envs + 1))]))
  gxe_env <- pop@gv[, -seq(1, (n_traits + n_traits * n_envs), n_envs + 1)]
  index <- as.list(as.data.frame(t(matrix(1:(n_traits * n_envs), ncol = n_traits))))
  gxe_env <- lapply(index, function(x) gxe_env[, x])

  gv <- lapply(gxe_env, function(x) x + do.call(cbind, g_main))
  gv <- do.call(rbind, mapply(function(x, y) x[rep(1:nrow(x), y), ], x = gv, y = as.list(n_reps), SIMPLIFY = F))
  colnames(gv) <- paste0("Trait.", 1:n_traits)

  compsym_asr <- data.frame(
    env = envs,
    rep = reps,
    id = pop@id,
    gv = gv
  )

  if (effects) {
    g_main <- mapply(cbind, list(data.frame(id = pop@id)), g_main = g_main, SIMPLIFY = F)

    gxe_env <- pop@gv[, -seq(1, (n_traits + n_traits * n_envs), n_envs + 1)]
    index_eff <- as.list(as.data.frame(t(matrix(1:(n_traits * n_envs), nrow = n_traits, byrow = TRUE))))
    gxe_env <- lapply(index_eff, function(x) gxe_env[, x])

    eff_comps <- mapply(cbind, g_main, gxe_env = gxe_env, SIMPLIFY = F)

    listNames <- c("Trial.df", paste0("Trait.", 1:n_traits))
    compsym_asr <- list(compsym_asr)
    compsym_asr <- c(compsym_asr, eff_comps)
    names(compsym_asr) <- listNames
  }

  return(compsym_asr)
}
