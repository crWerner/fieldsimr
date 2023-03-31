#' Genetic values based on an unstructured model for GxE interaction using 'AlphaSimR' - Input
#' parameters
#'
#' Creates a list of simulation parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} to simulate
#' genetic values for multiple traits in multiple environments based on an unstructured model for
#' genotype-by-environment (GxE) interaction. \cr
#' By default, 'AlphaSimR' does not support complex models for GxE interaction. However, its
#' functionality to simulate correlated genetic values can be utilised for this purpose by
#' providing the required variance structures. \code{unstr_asr_input} is a wrapper function to
#' construct the variance structures required to simulate GxE interaction in 'AlphaSimR' based on
#' an unstructured model. This function is also used in combination with the
#' wrapper function \link[FieldSimR]{compsym_asr_output}.
#'
#' \code{unstr_asr_input} can handle non-separable and separable structures between traits and
#' environments.
#' \itemize{
#'   \item For non-separable structures, provide (1) \code{var}, and (2) \code{cor_A}.
#'   \item For separable structures, provide (1) \code{T_var} & \code{E_var}, and (2)
#'   \code{T_cor_A} & \code{E_cor_A}. \cr
#'   }
#'
#' \strong{Note:} 'AlphaSimR' can simulate different biological effects (see:
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
#' @param n_traits Number of traits to be simulated.
#' @param mean A vector of mean genetic values for each trait-by-environment combination (ordered
#'   as environments within traits). Simulated traits can have a different mean for each
#'   environment. If the length of \code{mean} is equal to to \code{n_traits}, all traits will be
#'   assigned the same mean for each environment.
#' @param var A vector of genetic variances for each trait-by-environment combination (ordered as
#'   environments within traits). If the length of \code{var} is equal to \code{n_traits}, all
#'   traits will be assigned the same variance in each environment. \cr
#'   Alternatively, if a separable structure between traits and environments is desired,
#'   \code{T_var} and \code{E_var} can be provided. By default, \code{var = NULL}.
#' @param T_var A vector of genetic variances for each trait. Must be provided in combination with
#'   \code{E_var}. \cr
#'   Alternatively, \code{var} can be provided. By default, \code{T_var = NULL}.
#' @param E_var A vector of genetic variances for each environment. Must be provided in
#'   combination with \code{T_var}. \cr
#'   Alternatively, \code{var} can be provided. By default, \code{E_var = NULL}.
#' @param cor_A A matrix of additive genetic correlations between all trait-by-environment
#'   combinations. If not defined and \code{n_traits > 1}, a diagonal matrix is constructed. \cr
#'   Alternatively, \code{T_cor_A} and \code{E_cor_A} can be provided.
#' @param T_cor_A A matrix of additive genetic correlations between more than one trait. Must be
#'   provided in combination with \code{E_cor_A}. \cr
#'   Alternatively, \code{cor_A} can be provided. By default, \code{T_cor_A = NULL}.
#' @param E_cor_A A matrix of additive genetic correlations between more than one environment.
#'   Must be provided in combination with \code{T_cor_A}. \cr
#'   Alternatively, \code{cor_A} can be provided. By default, \code{E_cor_A = NULL}.
#' @param mean_DD A vector of mean dominance degrees for each trait-by-environment combination
#'   (ordered as environments within traits), similar to \code{mean}. If the length of
#'   \code{mean_DD} is equal to \code{n_traits}, all traits will be assigned the same mean_DD for
#'   each environment. By default, \code{mean_DD = NULL} and dominance is not simulated.
#' @param var_DD A vector of dominance degree variances for each trait-by-environment combination
#'   (ordered as environments within traits), similar to \code{var}. If the length of
#'   \code{var_DD} is equal to to \code{n_traits}, all traits will be assigned the same var_DD for
#'   each environment.\cr
#'   Alternatively, if a separable structure between traits and environments is desired,
#'   \code{T_var_DD} and \code{E_var_DD} can be provided. By default, \code{var_DD = NULL}.
#' @param T_var_DD A vector of dominance degree variances for each trait, similar to \code{T_var}.
#'   Must be provided in combination with \code{E_var_DD}. \cr
#'   Alternatively, \code{var_DD} can be provided. By default, \code{T_var_DD = NULL}.
#' @param E_var_DD A vector of dominance degree genetic variances for each environment, similar to
#'   \code{E_var}. Must be provided in combination with \code{T_var_DD}. \cr
#'   Alternatively, \code{var_DD} can be provided. By default, \code{E_var_DD = NULL}.
#' @param cor_DD A matrix of dominance degree correlations between all trait-by-environment
#'   combinations, similar to \code{cor_A}. If not defined and \code{n_traits > 1}, a diagonal
#'   matrix is constructed. \cr
#'   Alternatively, \code{T_cor_DD} and \code{E_cor_DD} can be provided. By default,
#'   \code{cor_DD = NULL}.
#' @param T_cor_DD A matrix of dominance degree correlations between more than one trait, similar
#'   to \code{T_cor_A}. Must be provided in combination with \code{E_cor_DD}. \cr
#'   Alternatively, \code{cor_DD} can be provided. By default, \code{T_cor_DD = NULL}.
#' @param E_cor_DD A matrix of dominance degree correlations between more than one environment,
#'   similar to \code{E_cor_A}. Must be provided in combination with \code{T_cor_DD}. \cr
#'   Alternatively, \code{cor_DD} can be provided. By default, \code{E_cor_DD = NULL}.
#' @param rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to additive genetic variance for each trait-by-environment combination (ordered as
#'   environments within traits), that is in a diploid organism with allele frequency 0.5. If the
#'   length of \code{rel_AA} is equal to to \code{n_traits}, all traits will be assigned the same
#'   rel_AA for each environment.\cr
#'   Alternatively, if a separable structure between traits and environments is desired,
#'   \code{T_rel_AA} and \code{E_rel_AA} can be provided. By default, \code{rel_AA = NULL} and
#'   epistasis is not simulated.
#' @param T_rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each trait, that is in a diploid organism with
#'   allele frequency 0.5. Must be provided in combination with \code{E_rel_AA}. \cr
#'   Alternatively, \code{rel_AA} can be provided. By default, \code{T_rel_AA = NULL}.
#' @param E_rel_AA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each environment, that is in a diploid organism
#'   with allele frequency 0.5. Must be provided in combination with \code{T_rel_AA}. \cr
#'   Alternatively, \code{rel_AA} can be provided. By default, \code{E_rel_AA = NULL}.
#' @param cor_AA A matrix of epistatic correlations between all trait-by-environment
#'   combinations, similar to \code{cor_A}. If not defined and \code{n_traits > 1},
#'   a diagonal matrix is constructed. \cr
#'   Alternatively, \code{T_cor_AA} and \code{E_cor_AA} can be provided. By default,
#'   \code{cor_AA = NULL}.
#' @param T_cor_AA A matrix of epistatic correlations between more than one trait, similar to
#'   \code{T_cor_A}. Must be provided in combination with \code{E_cor_AA}. \cr
#'   Alternatively, \code{cor_AA} can be provided. By default, \code{T_cor_AA = NULL}.
#' @param E_cor_AA A matrix of epistatic correlations between more than one environment, similar
#'   to \code{E_cor_A}. Must be provided in combination with \code{T_cor_AA}. \cr
#'   Alternatively, \code{cor_AA} can be provided. By default, \code{E_cor_AA = NULL}.
#'
#' @return A list containing input parameters for 'AlphaSimR', which is used to simulate
#'   correlated genetic effects based on an unstructured model.
#'
#' @examples
#' # Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested in
#' # three environments based on an unstructured model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees for trait 1 in all 3 environments and trait 2
#' # in all 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # Trait 1 x 3 environments, trait 2 x 3 environments.
#' mean_DD <- c(0.1, 0.4) # Trait 1 and 2, same values set in all 3 environments for each trait.
#'
#' # Additive genetic variances (useVarA = TRUE) and dominance degree variances for traits 1 and 2,
#' # assuming a separable structure between traits and environments.
#' T_var <- c(0.2, 10) # Genetic variances defined for the two traits.
#' E_var <- c(0.5, 1, 1.5) # Genetic variances defined for the three environments.
#'
#' # Dominance degree variances for trait 1 in 3 environments and for trait 2 in 3 environments,
#' # assuming a non-separable structure between traits and environments.
#' var_DD <- c(0.1, 0.15, 0.2, 0.2, 0.3, 0.2)
#'
#' # Additive genetic correlations between the two simulated traits.
#' T_cor_A <- matrix(
#'   c(
#'     1.0, 0.3,
#'     0.3, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the three simulated environments.
#' E_cor_A <- stats::cov2cor(matrix(
#'   c(
#'     0.5, 0.4, 0.6,
#'     0.4, 1.0, 0.5,
#'     0.6, 0.5, 1.5
#'   ),
#'   ncol = 3
#' ))
#'
#' # Dominance degree correlation between all six trait-by-environment combinations.
#' cor_DD <- diag(6) # Assuming independence between traits
#'
#' input_asr <- unstr_asr_input(
#'   n_envs = 3,
#'   n_traits = 2,
#'   mean = mean,
#'   T_var = T_var,
#'   E_var = E_var,
#'   T_cor_A = T_cor_A,
#'   E_cor_A = E_cor_A,
#'   mean_DD = mean_DD,
#'   var_DD = var_DD,
#'   cor_DD = cor_DD
#' )
#' @export
unstr_asr_input <- function(n_envs,
                            n_traits,
                            mean,
                            var = NULL,
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
      if (length(mean) == n_traits) {
        mean_pseudo <- rep(mean, each = n_envs)
      } else if (length(mean) == (n_traits * n_envs)) {
        mean_pseudo <- mean
      } else {
        stop("Number of values in argument 'mean' must either match number of
               traits or number of trait x environment combinations")
      }


      if (!is.null(var) & !is.null(E_var)) {
        stop("Either argument 'var' or 'E_var' must be provided")
      } else if (!is.null(var)) {
        if (!is.null(T_var)) {
          stop("Argument 'T_var' must be NULL if 'var' is provided")
        }
        if (length(var) == n_traits) {
          var_pseudo <- rep(var, each = n_envs)
        } else if (length(var) == (n_traits * n_envs)) {
          var_pseudo <- var
        } else {
          stop("Number of values in argument 'var' must match number of
                 trait x environment combinations")
        }
      } else if (!is.null(E_var)) {
        if (length(E_var) != n_envs) {
          stop("Number of values in argument 'E_var' must match number of
                 environments (n_envs)")
        }
        if (is.null(T_var)) stop("Argument 'T_var' is not provided")
        if (length(T_var) != n_traits) {
          stop("Number of values in argument 'T_var' must match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_var, each = n_envs) * rep(E_var, n_traits)
      } else {
        stop("Either argument 'var' or 'E_var' must be provided")
      }


      if (!is.null(cor_A) & !is.null(E_cor_A)) {
        stop("Either argument 'cor_A' or 'E_cor_A' must be provided")
      } else if (!is.null(cor_A)) {
        if (!is.null(T_cor_A)) {
          stop("Argument 'T_cor_A' must be NULL if 'cor_A' is provided")
        }
        if (nrow(cor_A) != length(mean_pseudo) | ncol(cor_A) != length(mean_pseudo)) {
          stop("Dimension of 'cor_A' does not match number of trait x
               environment combinations")
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
          stop("Argument 'T_cor_A' must be provided if in combination with
                   'E_cor_A' if n_traits > 1")
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
                 match number of trait x environment combinations")
        }
        if (any(unique(diag(E_cor_A)) != 1) | any(E_cor_A > 1) | any(E_cor_A < (-1))) {
          stop("'E_cor_A' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_A)) stop("'E_cor_A' is not symmetric")
        cor_pseudo <- kronecker(T_cor_A, E_cor_A)
      } else {
        stop("Either argument 'cor_A' or 'E_cor_A' must be provided")
      }


      input_asr <- list(
        mean = mean_pseudo,
        var = var_pseudo,
        cor_A = cor_pseudo
      )
    }

    if (i == "DD") {
      if (length(mean_DD) == n_traits) {
        mean_pseudo <- rep(mean_DD, each = n_envs)
      } else if (length(mean_DD) == (n_traits * n_envs)) {
        mean_pseudo <- mean_DD
      } else {
        stop("Number of values in argument 'mean_DD' must either match number
               of traits or number of trait x environment combinations")
      }


      if (!is.null(var_DD) & !is.null(E_var_DD)) {
        stop("Either argument 'var_DD' or 'E_var_DD' must be provided")
      } else if (!is.null(var_DD)) {
        if (!is.null(T_var_DD)) {
          stop("Argument 'T_var_DD' must be NULL if 'var_DD' is provided")
        }
        if (length(var_DD) == n_traits) {
          var_pseudo <- rep(var_DD, each = n_envs)
        } else if (length(var_DD) == (n_traits * n_envs)) {
          var_pseudo <- var_DD
        } else {
          stop("Number of values in argument 'var_DD' must match number of
                 trait x environment combinations")
        }
      } else if (!is.null(E_var_DD)) {
        if (length(E_var) != n_envs) {
          stop("Number of values in argument 'E_var_DD' must match number of
                 environments (n_envs)")
        }
        if (is.null(T_var_DD)) stop("Argument 'T_var_DD' is not provided")
        if (length(T_var_DD) != n_traits) {
          stop("Number of values in argument 'T_var_DD' must match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_var_DD, each = n_envs) * rep(E_var_DD, n_traits)
      } else {
        stop("Either argument 'var_DD' or 'E_var_DD' must be provided")
      }


      if (!is.null(cor_DD) & !is.null(E_cor_DD)) {
        stop("Either argument 'cor_DD' or 'E_cor_DD' must be provided")
      } else if (!is.null(cor_DD)) {
        if (!is.null(T_cor_DD)) {
          stop("Argument 'T_cor_DD' must be NULL if 'cor_DD' is provided")
        }
        if (nrow(cor_DD) != length(mean_pseudo) | ncol(cor_DD) != length(mean_pseudo)) {
          stop("Dimension of 'cor_DD' does not match number of trait x
               environment combinations")
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
          stop("Argument 'T_cor_DD' must be provided if in combination with
                   'E_cor_DD' if n_traits > 1")
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
                 not match number of trait x environment combinations")
        }
        if (any(unique(diag(E_cor_DD)) != 1) | any(E_cor_DD > 1) | any(E_cor_DD < (-1))) {
          stop("'E_cor_DD' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_DD)) stop("'E_cor_DD' is not symmetric")
        cor_pseudo <- kronecker(T_cor_DD, E_cor_DD)
      } else {
        stop("Either argument 'cor_DD' or 'E_cor_DD' must be provided")
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
      if (!is.null(rel_AA) & !is.null(E_rel_AA)) {
        stop("Either argument 'rel_AA' or 'E_rel_AA' must be provided")
      } else if (!is.null(rel_AA)) {
        if (!is.null(T_rel_AA)) {
          stop("Argument 'T_rel_AA' must be NULL if 'rel_AA' is provided")
        }
        if (length(rel_AA) == n_traits) {
          var_pseudo <- rep(rel_AA, each = n_envs)
        } else if (length(rel_AA) == (n_traits * n_envs)) {
          var_pseudo <- rel_AA
        } else {
          stop("Number of values in argument 'rel_AA' must match number of
               trait x environment combinations")
        }
      } else if (!is.null(E_rel_AA)) {
        if (length(E_rel_AA) != n_envs) {
          stop("Number of values in argument 'E_rel_AA' must match number of
                 environments (n_envs)")
        }
        if (is.null(T_rel_AA)) stop("Argument 'T_rel_AA' is not provided")
        if (length(T_rel_AA) != n_traits) {
          stop("Number of values in argument 'T_rel_AA' must match number of
                 traits (n_traits)")
        }
        var_pseudo <- rep(T_rel_AA, each = n_envs) * rep(E_rel_AA, n_traits)
      } else {
        stop("Either argument 'rel_AA' or 'E_rel_AA' must be provided")
      }


      if (!is.null(cor_AA) & !is.null(E_cor_AA)) {
        stop("Either argument 'cor_AA' or 'E_cor_AA' must be provided")
      } else if (!is.null(cor_AA)) {
        if (!is.null(T_cor_AA)) {
          stop("Argument 'T_cor_AA' must be NULL if 'cor_AA' is provided")
        }
        if (nrow(cor_AA) != length(mean_pseudo) | ncol(cor_AA) != length(mean_pseudo)) {
          stop("Dimension of 'cor_AA' does not match number of trait x
               environment combinations")
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
          stop("Argument 'T_cor_AA' must be provided if in combination with
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
                 match number of trait x environment combinations")
        }
        if (any(unique(diag(E_cor_AA)) != 1) | any(E_cor_AA > 1) | any(E_cor_AA < (-1))) {
          stop("'E_cor_AA' is not a correlation matrix")
        }
        if (!isSymmetric(E_cor_AA)) stop("'E_cor_AA' is not symmetric")
        cor_pseudo <- kronecker(T_cor_AA, E_cor_AA)
      } else {
        stop("Either argument 'cor_AA' or 'E_cor_AA' must be provided")
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

#' Genetic values based on an unstructured model for GxE interaction using 'AlphaSimR' - Simulated
#' genetic values
#'
#' Creates a data frame of correlated genetic values for multiple traits in multiple environments
#' based on an unstructured model for genotype-by-environment (GxE) interaction. This function
#' requires an 'AlphaSimR' population object generated using the \link[FieldSimR]{unstr_asr_input}
#' function.
#'
#' @param pop An 'AlphaSimR' population object (\code{\link[AlphaSimR]{Pop-class}} or
#'   \code{\link[AlphaSimR]{HybridPop-class}}) generated using \link[FieldSimR]{unstr_asr_input}.
#' @param n_envs Number of simulated environments (same as in \link[FieldSimR]{unstr_asr_input}).
#' @param n_reps A vector defining the number of complete replicates in each environment. If only
#'   one value is provided and \code{n_traits > 1}, all environments will be assigned the same
#'   number of replicates.
#' @param n_traits Number of simulated traits (same as in \link[FieldSimR]{unstr_asr_input}).
#'
#' @return A data frame containing the environment id, replicate number, genotype id, and the
#'   simulated genetic values for each trait.
#'
#' @examples
#' # Simulation of genetic values in 'AlphaSimR' for two additive + dominance traits tested in
#' # three environments based on an unstructured model for GxE interaction.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees for trait 1 in all 3 environments and trait 2
#' # in all 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # Trait 1 x 3 environments, trait 2 x 3 environments.
#' mean_DD <- c(0.1, 0.4) # Trait 1 and 2, same values set in all 3 environments for each trait.
#'
#' # Additive genetic variances (useVarA = TRUE) and dominance degree variances for traits 1 and 2,
#' # assuming a separable structure between traits and environments.
#' T_var <- c(0.2, 10) # Genetic variances defined for the two traits.
#' E_var <- c(0.5, 1, 1.5) # Genetic variances defined for the three environments.
#'
#' # Dominance degree variances for trait 1 in 3 environments and for trait 2 in 3 environments,
#' # assuming a non-separable structure between traits and environments.
#' var_DD <- c(0.1, 0.15, 0.2, 0.2, 0.3, 0.2)
#'
#' # Additive genetic correlations between the two simulated traits.
#' T_cor_A <- matrix(
#'   c(
#'     1.0, 0.3,
#'     0.3, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the three simulated environments.
#' E_cor_A <- stats::cov2cor(matrix(
#'   c(
#'     0.5, 0.4, 0.6,
#'     0.4, 1.0, 0.5,
#'     0.6, 0.5, 1.5
#'   ),
#'   ncol = 3
#' ))
#'
#' # Dominance degree correlation between all six trait-by-environment combinations.
#' cor_DD <- diag(6) # Assuming independence between traits
#'
#' input_asr <- unstr_asr_input(
#'   n_envs = 3,
#'   n_traits = 2,
#'   mean = mean,
#'   T_var = T_var,
#'   E_var = E_var,
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
#' # 3. Create a data frame containing the simulated genetic values for the two traits in the
#' # three environments.
#'
#' n_reps <- c(2, 3, 2) # Vector containing the number of complete replicates in each environment.
#'
#' gv_df <- unstr_asr_output(
#'   pop = pop,
#'   n_envs = 3,
#'   n_reps = n_reps,
#'   n_traits = 2
#' )
#' @export
unstr_asr_output <- function(pop,
                             n_envs,
                             n_reps,
                             n_traits) {
  if (n_envs < 2) stop("'n_envs' must be > 1")
  if (n_envs %% 1 != 0) stop("'n_envs' must be an integer")

  if ((sum(n_reps < 1) > 0) | (sum(n_reps %% 1 != 0) > 0)) {
    stop("'n_reps' must contain positive integer values")
  }
  if (length(n_reps) == 1) n_reps <- rep(n_reps, n_envs)

  if (n_traits < 1 | n_traits %% 1 != 0) stop("'n_traits' must be an integer > 0")

  envs <- rep(1:n_envs, times = length(pop@id) * n_reps)
  reps <- unlist(lapply(n_reps, function(x) rep(1:x, each = length(pop@id))))

  index <- as.list(as.data.frame(t(matrix(1:(n_traits * n_envs), ncol = n_traits))))
  gv <- lapply(index, function(x) pop@gv[, x])
  gv <- do.call(rbind, mapply(function(x, y) x[rep(1:nrow(x), y), ], x = gv, y = as.list(n_reps), SIMPLIFY = F))
  colnames(gv) <- paste0("Trait.", 1:n_traits)

  unstr_asr <- data.frame(
    env = envs,
    rep = reps,
    id = pop@id,
    gv = gv
  )

  return(unstr_asr)
}
