#' Simulate genetic values based on an unstructured model for GxE interaction - 'AlphaSimR' input
#' parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} to simulate
#' genetic values for multiple traits across multiple environments based on an unstructured model
#' for genotype-by-environment (GxE) interaction. \cr
#' By default, 'AlphaSimR' does not support complex models for GxE interaction. However, its
#' functionality to simulate correlated genetic values can be utilised for this purpose by
#' providing the required variance structures. \code{unstr_asr_input} is a wrapper function to
#' construct the variance structures required to simulate GxE interaction in 'AlphaSimR' based on
#' a multi-trait unstructured model. This function can handle separable and non-separable structures between traits and
#' environments (see below). After simulating the genetic values, the wrapper function
#' \link[FieldSimR]{unstr_asr_output} can be used to obtain a data frame with the relevant values.
#'
#' \code{unstr_asr_input} can handle separable and non-separable structures between traits and
#' environments.
#' \itemize{
#'   \item For separable structures, provide (1) \code{T_var} & \code{E_var}, and (2)
#'   \code{T_cor_A} & \code{E_cor_A}.
#'   \item For non-separable structures, provide (1) \code{var}, and (2) \code{cor_A}. \cr
#'   }
#'
#' \strong{Note:} 'AlphaSimR' can simulate different biological effects (see:
#' \href{https://gaynorr.github.io/AlphaSimR/reference/SimParam.html}{SimParam}).
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
#' @param n_traits Number of traits to be simulated.
#' @param mean A vector of mean genetic values for each environment-within-trait combination.
#'   If only one value is specified, all environment-within-trait combinations will be assigned
#'   the same mean.
#' @param var A vector of genetic variances for each environment-within-trait combination.
#'   If only one value is specified, all environment-within-trait combinations will be assigned
#'   the same mean. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{T_var} and \code{E_var} can be specified.
#' @param T_var A vector of genetic variances for each trait. Must be provided in combination with
#'   \code{E_var}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param E_var A vector of genetic variances for each environment. Must be provided in
#'   combination with \code{T_var}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param cor_A A matrix of additive genetic correlations between all environment-within-trait
#'   combinations. By default, a diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{T_cor_A} and \code{E_cor_A} can be specified.
#' @param T_cor_A A matrix of additive genetic correlations between traits. Must be provided in
#'   combination with \code{E_cor_A}. \cr
#'   \strong{Alternatively}, \code{cor_A} can be specified.
#' @param E_cor_A A matrix of additive genetic correlations between environments.
#'   Must be provided in combination with \code{T_cor_A}. \cr
#'   \strong{Alternatively}, \code{cor_A} can be specified.
#' @param mean_DD A vector of mean dominance degrees for each environment-within-trait combination
#'   (similar to \code{mean}). If only one value is specified, all environment-within-trait
#'   combinations will be assigned the same mean. By default, \code{mean_DD = NULL} and dominance
#'   is not simulated.
#' @param var_DD A vector of dominance degree variances for each environment-within-trait
#'   combination (similar to \code{var}). If only one value is specified, all environment-within-trait
#'   combinations will be assigned the same variance. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{T_var_DD} and \code{E_var_DD} can be specified.
#' @param T_var_DD A vector of dominance degree variances for each trait (similar to \code{T_var}).
#'   Must be provided in combination with \code{E_var_DD}. \cr
#'   \strong{Alternatively}, \code{var_DD} can be specified.
#' @param E_var_DD A vector of dominance degree genetic variances for each environment (similar to
#'   \code{E_var}). Must be provided in combination with \code{T_var_DD}. \cr
#'   \strong{Alternatively}, \code{var_DD} can be specified.
#' @param cor_DD A matrix of dominance degree correlations between all environment-within-trait
#'   combinations (similar to \code{cor_A}). If not specified and dominance is simulated, a
#'   diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{T_cor_DD} and \code{E_cor_DD} can be specified.
#' @param T_cor_DD A matrix of dominance degree correlations between traits (similar to
#'   \code{T_cor_A}). Must be provided in combination with \code{E_cor_DD}. \cr
#'   \strong{Alternatively}, \code{cor_DD} can be specified.
#' @param E_cor_DD A matrix of dominance degree correlations between environments (similar to
#'   \code{E_cor_A}). Must be provided in combination with \code{T_cor_DD}. \cr
#'   \strong{Alternatively}, \code{cor_DD} can be specified.
#' @param rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to additive genetic variance for each environment-within-trait combination,
#'   that is in a diploid organism with allele frequency 0.5. If only one value is specified,
#'   all environment-within-trait combinations will be assigned the same value. By default,
#'   \code{rel_AA = NULL} and epistasis is not simulated. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{T_rel_AA} and \code{E_rel_AA} can be specified.
#' @param T_rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each trait. Must be provided in combination
#'   with \code{E_rel_AA}. \cr
#'   \strong{Alternatively}, \code{rel_AA} can be specified.
#' @param E_rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each environment. Must be provided in
#'   combination with \code{T_rel_AA}. \cr
#'   \strong{Alternatively}, \code{rel_AA} can be specified.
#' @param cor_AA A matrix of epistatic correlations between all environment-within-trait
#'   combinations (similar to \code{cor_A}). If not specified and epistasis is simulated, a
#'   diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{T_cor_AA} and \code{E_cor_AA} can be specified.
#' @param T_cor_AA A matrix of epistatic correlations between traits (similar to \code{T_cor_A}).
#'   Must be provided in combination with \code{E_cor_AA}. \cr
#'   \strong{Alternatively}, \code{cor_AA} can be specified.
#' @param E_cor_AA A matrix of epistatic correlations between environments (similar to
#'   \code{E_cor_A}). Must be provided in combination with \code{T_cor_AA}. \cr
#'   \strong{Alternatively}, \code{cor_AA} can be specified.
#'
#' @return A list containing input parameters for 'AlphaSimR', which is used to simulate
#'   correlated genetic effects based on an unstructured model.
#'
#' @examples
#' # Simulate genetic values in 'AlphaSimR' for two additive + dominance traits across
#' # two environments based on an unstructured model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' mean_DD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same values in the two environments.
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.085, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' var_DD <- rep(0.2, 4) # The same value set for traits 1 and 2.
#'
#' # Additive genetic correlations between the two simulated traits.
#' T_cor_A <- matrix(
#'   c(
#'     1.0, 0.6,
#'     0.6, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the two simulated environments.
#' E_cor_A <- matrix(
#'   c(
#'     1.0, 0.2,
#'     0.2, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' cor_DD <- diag(4) # Assuming independence between traits
#'
#' input_asr <- unstr_asr_input(
#'   n_envs = 2,
#'   n_traits = 2,
#'   mean = mean,
#'   var = var,
#'   T_cor_A = T_cor_A,
#'   E_cor_A = E_cor_A,
#'   mean_DD = mean_DD,
#'   var_DD = var_DD,
#'   cor_DD = cor_DD
#' )
#' @export
unstr_asr_input <- function(n_envs = 3,
                            n_traits = 2,
                            mean = 0,
                            var = 1,
                            T_var = NULL,
                            E_var = NULL,
                            cor_A = NULL,
                            T_cor_A = NULL,
                            E_cor_A = NULL,
                            mean_DD = NULL,
                            var_DD = NULL,
                            T_var_DD = NULL,
                            E_var_DD = NULL,
                            cor_DD = NULL,
                            T_cor_DD = NULL,
                            E_cor_DD = NULL,
                            rel_AA = NULL,
                            T_rel_AA = NULL,
                            E_rel_AA = NULL,
                            cor_AA = NULL,
                            T_cor_AA = NULL,
                            E_cor_AA = NULL) {
  if (n_envs < 2) stop("'n_envs' must be > 1")
  if (n_envs %% 1 != 0) stop("'n_envs' must be an integer")
  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")


  if (is.null(mean_DD) & is.null(var_DD) & is.null(E_var_DD) & is.null(rel_AA) & is.null(E_rel_AA)) {
    labels <- "A"
  } else if (is.null(mean_DD) & is.null(var_DD) & is.null(E_var_DD)) {
    labels <- c("A", "AA")
  } else if (is.null(rel_AA) & is.null(E_rel_AA)) {
    labels <- c("A", "DD")
  } else {
    labels <- c("A", "DD", "AA")
  }


  for (i in labels) {
    if (i == "A") {
      if (length(mean) == 1) {
        mean_pseudo <- rep(mean, each = n_traits * n_envs)
      } else if (length(mean) == (n_traits * n_envs)) {
        mean_pseudo <- mean
      } else {
        stop("Number of values in 'mean' must be 1 or match
              number of environment-within-trait combinations")
      }


      if (is.null(var) & !(!is.null(E_var) & !is.null(T_var))) {
        stop("Either 'var' or 'E_var' and 'T_var' must be specified")
      } else if (!is.null(var)) {
        if (!is.null(E_var) | !is.null(T_var)) {
          stop("'E_var' and 'T_var' must be NULL if 'var' is specified")
        }
        if (length(var) == 1) {
          var_pseudo <- rep(var, each = n_traits * n_envs)
        } else if (length(var) == (n_traits * n_envs)) {
          var_pseudo <- var
        } else {
          stop("Number of values in 'var' must be 1 or match number of
                 environment-within-trait combinations")
        }
      } else if (!is.null(E_var)) {
        if (length(E_var) != n_envs) {
          stop("Number of values in 'E_var' does not match number of
                 environments (n_envs)")
        }
        if (is.null(T_var)) stop("'T_var' must be specified")
        if (length(T_var) != n_traits) {
          stop("Number of values in 'T_var' does not match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_var, each = n_envs) * rep(E_var, n_traits)
      } else {
        stop("Either 'var' or 'E_var' must be specified")
      }

      if (is.null(cor_A) & !(!is.null(E_cor_A) & !is.null(T_cor_A))) {
        stop("Either 'cor_A' or 'E_cor_A' and 'T_cor_A' must be specified")
      } else if (!is.null(cor_A)) {
        if (!is.null(E_cor_A) | !is.null(T_cor_A)) {
          stop("'E_cor_A' and 'T_cor_A' must be NULL if 'cor_A' is specified")
        }
        if (nrow(cor_A) != length(mean_pseudo) | ncol(cor_A) != length(mean_pseudo)) {
          stop("Dimension of 'cor_A' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(cor_A)) != 1) | any(cor_A > 1) | any(cor_A < (-1))) {
          stop("'cor_A' is not a correlation matrix")
        }
        if (!isSymmetric(cor_A)) stop("cor_A is not symmetric")
        cor_pseudo <- cor_A
      } else if (!is.null(E_cor_A)) {
        if (is.null(T_cor_A) & n_traits == 1) {
          T_cor_A <- matrix(1)
        } else if (is.null(T_cor_A) & n_traits > 1) {
          stop("'T_cor_A' must be specified in combination with
                   'E_cor_A' if 'n_traits' > 1")
        } else if (length(T_cor_A) == 1) {
          T_cor_A <- matrix(1)
        } else {
          if (any(unique(diag(T_cor_A)) != 1) | any(T_cor_A > 1) | any(T_cor_A < (-1))) {
            stop("'T_cor_A' is not a correlation matrix")
          }
          if (!isSymmetric(T_cor_A)) stop("'T_cor_A' is not symmetric")
        }
        if ((nrow(E_cor_A) * nrow(T_cor_A)) != length(mean_pseudo) | (ncol(E_cor_A) * ncol(T_cor_A)) != length(mean_pseudo)) {
          stop("Dimension of the 'E_cor_A' x 'T_cor_A' Kronecker product does not
                 match number of environment-within-trait combinations")
        }
        if (any(unique(diag(E_cor_A)) != 1) | any(E_cor_A > 1) | any(E_cor_A < (-1))) {
          stop("'E_cor_A' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_A)) stop("'E_cor_A' is not symmetric")
        cor_pseudo <- kronecker(T_cor_A, E_cor_A)
      } else {
        stop("Either 'cor_A' or 'E_cor_A' must be specified")
      }


      input_asr <- list(
        mean = mean_pseudo,
        var = var_pseudo,
        cor_A = cor_pseudo
      )
    }

    if (i == "DD") {
      if (length(mean_DD) == 1) {
        mean_pseudo <- rep(mean_DD, each = n_traits * n_envs)
      } else if (length(mean_DD) == (n_traits * n_envs)) {
        mean_pseudo <- mean_DD
      } else {
        stop("Number of values in 'mean_DD' must be 1 or match number of
              environment-within-trait combinations")
      }

      if (is.null(var_DD) & !(!is.null(E_var_DD) & !is.null(T_var_DD))) {
        stop("Either 'var_DD' or 'E_var_DD' and 'T_var_DD' must be specified")
      } else if (!is.null(var_DD)) {
        if (!is.null(E_var_DD) | !is.null(T_var_DD)) {
          stop("'E_var_DD' and 'T_var_DD' must be NULL if 'var_DD' is specified")
        }
        if (length(var_DD) == 1) {
          var_pseudo <- rep(var_DD, each = n_traits * n_envs)
        } else if (length(var_DD) == (n_traits * n_envs)) {
          var_pseudo <- var_DD
        } else {
          stop("Number of values in 'var_DD' must be 1 or match number of
                environment-within-trait combinations")
        }
      } else if (!is.null(E_var_DD)) {
        if (length(E_var) != n_envs) {
          stop("Number of values in 'E_var_DD' does not match number of
                 environments (n_envs)")
        }
        if (is.null(T_var_DD)) stop("'T_var_DD' must be specified")
        if (length(T_var_DD) != n_traits) {
          stop("Number of values in 'T_var_DD' does not must match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_var_DD, each = n_envs) * rep(E_var_DD, n_traits)
      } else {
        stop("Either 'var_DD' or 'E_var_DD' must be specified")
      }


      if (is.null(cor_DD) & !(!is.null(E_cor_DD) & !is.null(T_cor_DD))) {
        stop("Either 'cor_DD' or 'E_cor_DD' and 'T_cor_DD' must be specified")
      } else if (!is.null(cor_DD)) {
        if (!is.null(E_cor_DD) | !is.null(T_cor_DD)) {
          stop("'E_cor_DD' and 'T_cor_DD' must be NULL if 'cor_DD' is specified")
        }
        if (nrow(cor_DD) != length(mean_pseudo) | ncol(cor_DD) != length(mean_pseudo)) {
          stop("Dimension of 'cor_DD' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(cor_DD)) != 1) | any(cor_DD > 1) | any(cor_DD < (-1))) {
          stop("'cor_DD' is not a correlation matrix")
        }
        if (!isSymmetric(cor_DD)) stop("cor_DD is not symmetric")
        cor_pseudo <- cor_DD
      } else if (!is.null(E_cor_DD)) {
        if (is.null(T_cor_DD) & n_traits == 1) {
          T_cor_DD <- matrix(1)
        } else if (is.null(T_cor_DD) & n_traits > 1) {
          stop("'T_cor_DD' must be specified in combination with
                   'E_cor_DD'  if 'n_traits' > 1")
        } else if (length(T_cor_DD) == 1) {
          T_cor_DD <- matrix(1)
        } else {
          if (any(unique(diag(T_cor_DD)) != 1) | any(T_cor_DD > 1) | any(T_cor_DD < (-1))) {
            stop("'T_cor_DD' is not a correlation matrix")
          }
          if (!isSymmetric(T_cor_DD)) stop("'T_cor_DD' is not symmetric")
        }
        if ((nrow(E_cor_DD) * nrow(T_cor_DD)) != length(mean_pseudo) | (ncol(E_cor_DD) * ncol(T_cor_DD)) != length(mean_pseudo)) {
          stop("Dimension of the 'E_cor_DD' x 'T_cor_DD' Kronecker product does
                 not match number of environment-within-trait combinations")
        }
        if (any(unique(diag(E_cor_DD)) != 1) | any(E_cor_DD > 1) | any(E_cor_DD < (-1))) {
          stop("'E_cor_DD' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_DD)) stop("'E_cor_DD' is not symmetric")
        cor_pseudo <- kronecker(T_cor_DD, E_cor_DD)
      } else {
        stop("Either 'cor_DD' or 'E_cor_DD' must be specified")
      }


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
      if (is.null(rel_AA) & !(!is.null(E_rel_AA) & !is.null(T_rel_AA))) {
        stop("Either 'rel_AA' or 'E_rel_AA' and 'T_rel_AA' must be specified")
      } else if (!is.null(rel_AA)) {
        if (!is.null(E_rel_AA) | !is.null(T_rel_AA)) {
          stop("'E_rel_AA' and 'T_rel_AA' must be NULL if 'rel_AA' is specified")
        }
        if (length(rel_AA) == n_traits) {
          var_pseudo <- rep(rel_AA, each = n_envs)
        } else if (length(rel_AA) == (n_traits * n_envs)) {
          var_pseudo <- rel_AA
        } else {
          stop("Number of values in argument 'rel_AA' must match number of
               environment-within-trait combinations")
        }
      } else if (!is.null(E_rel_AA)) {
        if (length(E_rel_AA) != n_envs) {
          stop("Number of values in argument 'E_rel_AA' must match number of
                 environments (n_envs)")
        }
        if (is.null(T_rel_AA)) stop("Argument 'T_rel_AA' must be specified")
        if (length(T_rel_AA) != n_traits) {
          stop("Number of values in argument 'T_rel_AA' must match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_rel_AA, each = n_envs) * rep(E_rel_AA, n_traits)
      } else {
        stop("Either argument 'rel_AA' or 'E_rel_AA' must be specified")
      }

      if (is.null(cor_AA) & !(!is.null(E_cor_AA) & !is.null(T_cor_AA))) {
        stop("Either 'cor_AA' or 'E_cor_AA' and 'T_cor_AA' must be specified")
      } else if (!is.null(cor_AA)) {
        if (!is.null(E_cor_AA) | !is.null(T_cor_AA)) {
          stop("'E_cor_AA' and 'T_cor_AA' must be NULL if 'cor_AA' is specified")
        }
        if (nrow(cor_AA) != length(mean_pseudo) | ncol(cor_AA) != length(mean_pseudo)) {
          stop("Dimension of 'cor_AA' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(cor_AA)) != 1) | any(cor_AA > 1) | any(cor_AA < (-1))) {
          stop("'cor_AA' is not a correlation matrix")
        }
        if (!isSymmetric(cor_AA)) stop("cor_AA is not symmetric")
        cor_pseudo <- cor_AA
      } else if (!is.null(E_cor_AA)) {
        if (is.null(T_cor_AA) & n_traits == 1) {
          T_cor_AA <- matrix(1)
        } else if (is.null(T_cor_AA) & n_traits > 1) {
          stop("Argument 'T_cor_AA' must be specified in combination with
                   'E_cor_AA' if n_traits > 1")
        } else if (length(T_cor_AA) == 1) {
          T_cor_AA <- matrix(1)
        } else {
          if (any(unique(diag(T_cor_AA)) != 1) | any(T_cor_AA > 1) | any(T_cor_AA < (-1))) {
            stop("'T_cor_AA' is not a correlation matrix")
          }
          if (!isSymmetric(T_cor_AA)) stop("'T_cor_AA' is not symmetric")
        }
        if ((nrow(E_cor_AA) * nrow(T_cor_AA)) != length(mean_pseudo) | (ncol(E_cor_AA) * ncol(T_cor_AA)) != length(mean_pseudo)) {
          stop("Dimension of the 'E_cor_AA' x 'T_cor_AA' Kronecker product does not
                 match number of environment-within-trait combinations")
        }
        if (any(unique(diag(E_cor_AA)) != 1) | any(E_cor_AA > 1) | any(E_cor_AA < (-1))) {
          stop("'E_cor_AA' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_AA)) stop("'E_cor_AA' is not symmetric")
        cor_pseudo <- kronecker(T_cor_AA, E_cor_AA)
      } else {
        stop("Either argument 'cor_AA' or 'E_cor_AA' must be specified")
      }


      input_asr <- c(
        input_asr,
        list(
          meanAA = rep(0, n_envs * n_traits),
          rel_AA = var_pseudo,
          cor_AA = cor_pseudo
        )
      )
    }
  }

  return(input_asr)
}

#' Simulate genetic values based on an unstructured model for GxE interaction -
#' Simulation using 'AlphaSimR'
#'
#' Creates a data frame of simulated genetic values for multiple traits across multiple environments
#' based on an unstructured model for genotype-by-environment (GxE) interaction. This function
#' requires an \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} population object
#' generated using the \link[FieldSimR]{unstr_asr_input} function.
#'
#' @param pop An \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} population object
#'   (\href{https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html}{Pop-class} or
#'   \href{https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html}{HybridPop-class})
#'   generated using \link[FieldSimR]{unstr_asr_input}.
#' @param n_envs Number of simulated environments (same number used in \link[FieldSimR]{unstr_asr_input}).
#' @param n_traits Number of simulated traits (same number used in \link[FieldSimR]{unstr_asr_input}).
#' @param n_reps A vector defining the number of replicates in each environment. If only one value
#'   is specified, all environments will be assigned the same number.
#'
#' @return  A data frame with columns 'env', 'rep', genotype 'id', and the
#'   simulated genetic values for each trait.
#'
#' @examples
#' # Simulate genetic values in 'AlphaSimR' for two additive + dominance traits across
#' # two environments based on an unstructured model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' mean_DD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same values in the two environments.
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.085, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' var_DD <- rep(0.2, 4) # The same value set for traits 1 and 2.
#'
#' # Additive genetic correlations between the two simulated traits.
#' T_cor_A <- matrix(
#'   c(
#'     1.0, 0.6,
#'     0.6, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the two simulated environments.
#' E_cor_A <- matrix(
#'   c(
#'     1.0, 0.2,
#'     0.2, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' cor_DD <- diag(4) # Assuming independence between traits
#'
#' input_asr <- unstr_asr_input(
#'   n_envs = 2,
#'   n_traits = 2,
#'   mean = mean,
#'   var = var,
#'   T_cor_A = T_cor_A,
#'   E_cor_A = E_cor_A,
#'   mean_DD = mean_DD,
#'   var_DD = var_DD,
#'   cor_DD = cor_DD
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on an unstructured model for
#' # GxE interaction.
#'
#' library("AlphaSimR")
#' FOUNDERPOP <- quickHaplo(
#'   nInd = 10,
#'   nChr = 1,
#'   segSites = 20
#' )
#'
#' SP <- SimParam$new(FOUNDERPOP)
#'
#' \dontshow{
#' SP$nThreads <- 1L
#' }
#'
#' SP$addTraitAD(
#'   nQtlPerChr = 20,
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
#' # If useVarA = FALSE, 'var' represents the total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for the two traits across the
#' # two environments.
#'
#' n_reps <- c(2, 2) # Vector containing the number of complete replicates in each environment.
#'
#' gv_df <- unstr_asr_output(
#'   pop = pop,
#'   n_envs = 2,
#'   n_traits = 2,
#'   n_reps = n_reps
#' )
#' @export
unstr_asr_output <- function(pop,
                             n_envs,
                             n_traits,
                             n_reps) {
  if (n_envs < 2) stop("'n_envs' must be > 1")
  if (n_envs %% 1 != 0) stop("'n_envs' must be an integer")

  if ((sum(n_reps < 1) > 0) | (sum(n_reps %% 1 != 0) > 0)) {
    stop("'n_reps' must contain positive integer values")
  }
  if (length(n_reps) == 1) n_reps <- rep(n_reps, n_envs)

  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  envs <- factor(rep(1:n_envs, times = length(pop@id) * n_reps))
  reps <- factor(unlist(lapply(n_reps, function(x) rep(1:x, each = length(pop@id)))))
  ids <- factor(as.numeric(as.character(pop@id)))

  index <- as.list(as.data.frame(t(matrix(1:(n_traits * n_envs), ncol = n_traits))))
  gv <- lapply(index, function(x) cbind(pop@gv[, x]))
  gv <- do.call(rbind, mapply(function(x, y) cbind(x[rep(1:nrow(x), y), ]), x = gv, y = as.list(n_reps), SIMPLIFY = F))
  colnames(gv) <- paste0("Trait.", 1:n_traits)

  unstr_asr <- data.frame(
    env = envs,
    rep = reps,
    id = ids,
    gv
  )
  unstr_asr <- unstr_asr[order(unstr_asr$env, unstr_asr$rep, unstr_asr$id), ]

  return(unstr_asr)
}
