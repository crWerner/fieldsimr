#' Simulate genetic values based on an unstructured model for GxE interaction - 'AlphaSimR' input
#' parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} to simulate
#' genetic values for multiple traits across multiple environments based on an unstructured model
#' for genotype-by-environment (GxE) interaction. \cr
#' By default, 'AlphaSimR' does not support complex models for GxE interaction. However, its
#' functionality to simulate correlated genetic values can be utilised for this purpose by
#' providing the required variance structures. \code{output_asr_input} is a wrapper function to
#' construct the variance structures required to simulate GxE interaction in 'AlphaSimR' based on
#' a multi-trait unstructured model. This function can handle separable and non-separable structures between traits and
#' environments (see below). After simulating the genetic values, the wrapper function
#' \link[FieldSimR]{output_asr_output} can be used to obtain a data frame with the relevant values.
#'
#' \code{output_asr_input} can handle separable and non-separable structures between traits and
#' environments.
#' \itemize{
#'   \item For separable structures, provide (1) \code{Tvar} & \code{Evar}, and (2)
#'   \code{TcorA} & \code{EcorA}.
#'   \item For non-separable structures, provide (1) \code{var}, and (2) \code{corA}. \cr
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
#' @param nenvs Number of environments to be simulated. A minimum of two environments is required.
#' @param ntraits Number of traits to be simulated.
#' @param mean A vector of mean genetic values for each environment-within-trait combination.
#'   If only one value is specified, all environment-within-trait combinations will be assigned
#'   the same mean.
#' @param var A vector of genetic variances for each environment-within-trait combination.
#'   If only one value is specified, all environment-within-trait combinations will be assigned
#'   the same mean. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{Tvar} and \code{Evar} can be specified.
#' @param Tvar A vector of genetic variances for each trait. Must be provided in combination with
#'   \code{Evar}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param Evar A vector of genetic variances for each environment. Must be provided in
#'   combination with \code{Tvar}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param corA A matrix of additive genetic correlations between all environment-within-trait
#'   combinations. By default, a diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{TcorA} and \code{EcorA} can be specified.
#' @param TcorA A matrix of additive genetic correlations between traits. Must be provided in
#'   combination with \code{EcorA}. \cr
#'   \strong{Alternatively}, \code{corA} can be specified.
#' @param EcorA A matrix of additive genetic correlations between environments.
#'   Must be provided in combination with \code{TcorA}. \cr
#'   \strong{Alternatively}, \code{corA} can be specified.
#' @param meanDD A vector of mean dominance degrees for each environment-within-trait combination
#'   (similar to \code{mean}). If only one value is specified, all environment-within-trait
#'   combinations will be assigned the same mean. By default, \code{meanDD = NULL} and dominance
#'   is not simulated.
#' @param varDD A vector of dominance degree variances for each environment-within-trait
#'   combination (similar to \code{var}). If only one value is specified, all environment-within-trait
#'   combinations will be assigned the same variance. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{TvarDD} and \code{EvarDD} can be specified.
#' @param TvarDD A vector of dominance degree variances for each trait (similar to \code{Tvar}).
#'   Must be provided in combination with \code{EvarDD}. \cr
#'   \strong{Alternatively}, \code{varDD} can be specified.
#' @param EvarDD A vector of dominance degree genetic variances for each environment (similar to
#'   \code{Evar}). Must be provided in combination with \code{TvarDD}. \cr
#'   \strong{Alternatively}, \code{varDD} can be specified.
#' @param corDD A matrix of dominance degree correlations between all environment-within-trait
#'   combinations (similar to \code{corA}). If not specified and dominance is simulated, a
#'   diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{TcorDD} and \code{EcorDD} can be specified.
#' @param TcorDD A matrix of dominance degree correlations between traits (similar to
#'   \code{TcorA}). Must be provided in combination with \code{EcorDD}. \cr
#'   \strong{Alternatively}, \code{corDD} can be specified.
#' @param EcorDD A matrix of dominance degree correlations between environments (similar to
#'   \code{EcorA}). Must be provided in combination with \code{TcorDD}. \cr
#'   \strong{Alternatively}, \code{corDD} can be specified.
#' @param relAA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to additive genetic variance for each environment-within-trait combination,
#'   that is in a diploid organism with allele frequency 0.5. If only one value is specified,
#'   all environment-within-trait combinations will be assigned the same value. By default,
#'   \code{relAA = NULL} and epistasis is not simulated. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{TrelAA} and \code{ErelAA} can be specified.
#' @param TrelAA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each trait. Must be provided in combination
#'   with \code{ErelAA}. \cr
#'   \strong{Alternatively}, \code{relAA} can be specified.
#' @param ErelAA A vector defining the magnitude of additive-by-additive (epistatic) variance
#'   relative to the additive genetic variance for each environment. Must be provided in
#'   combination with \code{TrelAA}. \cr
#'   \strong{Alternatively}, \code{relAA} can be specified.
#' @param corAA A matrix of epistatic correlations between all environment-within-trait
#'   combinations (similar to \code{corA}). If not specified and epistasis is simulated, a
#'   diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{TcorAA} and \code{EcorAA} can be specified.
#' @param TcorAA A matrix of epistatic correlations between traits (similar to \code{TcorA}).
#'   Must be provided in combination with \code{EcorAA}. \cr
#'   \strong{Alternatively}, \code{corAA} can be specified.
#' @param EcorAA A matrix of epistatic correlations between environments (similar to
#'   \code{EcorA}). Must be provided in combination with \code{TcorAA}. \cr
#'   \strong{Alternatively}, \code{corAA} can be specified.
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
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same values in the two environments.
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.085, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' varDD <- rep(0.2, 4) # The same value set for traits 1 and 2.
#'
#' # Additive genetic correlations between the two simulated traits.
#' TcorA <- matrix(
#'   c(
#'     1.0, 0.6,
#'     0.6, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the two simulated environments.
#' EcorA <- matrix(
#'   c(
#'     1.0, 0.2,
#'     0.2, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' corDD <- diag(4) # Assuming independence between traits
#'
#' input_asr <- output_asr_input(
#'   nenvs = 2,
#'   ntraits = 2,
#'   mean = mean,
#'   var = var,
#'   TcorA = TcorA,
#'   EcorA = EcorA,
#'   meanDD = meanDD,
#'   varDD = varDD,
#'   corDD = corDD
#' )
#' @export
output_asr_input <- function(ntraits = 1,
                             nenvs = 2,
                            mean = 0,
                            var = 1,
                            Tvar = NULL,
                            Evar = NULL,
                            corA = NULL,
                            TcorA = NULL,
                            EcorA = NULL,
                            meanDD = NULL,
                            varDD = NULL,
                            TvarDD = NULL,
                            EvarDD = NULL,
                            corDD = NULL,
                            TcorDD = NULL,
                            EcorDD = NULL,
                            relAA = NULL,
                            TrelAA = NULL,
                            ErelAA = NULL,
                            corAA = NULL,
                            TcorAA = NULL,
                            EcorAA = NULL) {
  if (nenvs < 2) stop("'nenvs' must be > 1")
  if (nenvs %% 1 != 0) stop("'nenvs' must be an integer")
  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be an integer > 0")


  if (is.null(meanDD) & is.null(varDD) & is.null(EvarDD) & is.null(relAA) & is.null(ErelAA)) {
    labels <- "A"
  } else if (is.null(meanDD) & is.null(varDD) & is.null(EvarDD)) {
    labels <- c("A", "AA")
  } else if (is.null(relAA) & is.null(ErelAA)) {
    labels <- c("A", "DD")
  } else {
    labels <- c("A", "DD", "AA")
  }


  for (i in labels) {
    if (i == "A") {
      if (length(mean) == 1) {
        mean_pseudo <- rep(mean, each = ntraits * nenvs)
      } else if (length(mean) == (ntraits * nenvs)) {
        mean_pseudo <- mean
      } else {
        stop("Number of values in 'mean' must be 1 or match
              number of environment-within-trait combinations")
      }


      if (is.null(var) & !(!is.null(Evar) & !is.null(Tvar))) {
        stop("Either 'var' or 'Evar' and 'Tvar' must be specified")
      } else if (!is.null(var)) {
        if (!is.null(Evar) | !is.null(Tvar)) {
          stop("'Evar' and 'Tvar' must be NULL if 'var' is specified")
        }
        if (length(var) == 1) {
          var_pseudo <- rep(var, each = ntraits * nenvs)
        } else if (length(var) == (ntraits * nenvs)) {
          var_pseudo <- var
        } else {
          stop("Number of values in 'var' must be 1 or match number of
                 environment-within-trait combinations")
        }
      } else if (!is.null(Evar)) {
        if (length(Evar) != nenvs) {
          stop("Number of values in 'Evar' does not match number of
                 environments (nenvs)")
        }
        if (is.null(Tvar)) stop("'Tvar' must be specified")
        if (length(Tvar) != ntraits) {
          stop("Number of values in 'Tvar' does not match number of
                 traits (ntraits)")
        }
        var_pseudo <- rep(Tvar, each = nenvs) * rep(Evar, ntraits)
      } else {
        stop("Either 'var' or 'Evar' must be specified")
      }

      if (is.null(corA) & !(!is.null(EcorA) & !is.null(TcorA))) {
        stop("Either 'corA' or 'EcorA' and 'TcorA' must be specified")
      } else if (!is.null(corA)) {
        if (!is.null(EcorA) | !is.null(TcorA)) {
          stop("'EcorA' and 'TcorA' must be NULL if 'corA' is specified")
        }
        if (nrow(corA) != length(mean_pseudo) | ncol(corA) != length(mean_pseudo)) {
          stop("Dimension of 'corA' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(corA)) != 1) | any(corA > 1) | any(corA < (-1))) {
          stop("'corA' is not a correlation matrix")
        }
        if (!isSymmetric(corA)) stop("corA is not symmetric")
        cor_pseudo <- corA
      } else if (!is.null(EcorA)) {
        if (is.null(TcorA) & ntraits == 1) {
          TcorA <- matrix(1)
        } else if (is.null(TcorA) & ntraits > 1) {
          stop("'TcorA' must be specified in combination with
                   'EcorA' if 'ntraits' > 1")
        } else if (length(TcorA) == 1) {
          TcorA <- matrix(1)
        } else {
          if (any(unique(diag(TcorA)) != 1) | any(TcorA > 1) | any(TcorA < (-1))) {
            stop("'TcorA' is not a correlation matrix")
          }
          if (!isSymmetric(TcorA)) stop("'TcorA' is not symmetric")
        }
        if ((nrow(EcorA) * nrow(TcorA)) != length(mean_pseudo) | (ncol(EcorA) * ncol(TcorA)) != length(mean_pseudo)) {
          stop("Dimension of the 'EcorA' x 'TcorA' Kronecker product does not
                 match number of environment-within-trait combinations")
        }
        if (any(unique(diag(EcorA)) != 1) | any(EcorA > 1) | any(EcorA < (-1))) {
          stop("'EcorA' is not a correlation matrix")
        }
        if (!isSymmetric(EcorA)) stop("'EcorA' is not symmetric")
        cor_pseudo <- kronecker(TcorA, EcorA)
      } else {
        stop("Either 'corA' or 'EcorA' must be specified")
      }


      input_asr <- list(
        mean = mean_pseudo,
        var = var_pseudo,
        corA = cor_pseudo
      )
    }

    if (i == "DD") {
      if (length(meanDD) == 1) {
        mean_pseudo <- rep(meanDD, each = ntraits * nenvs)
      } else if (length(meanDD) == (ntraits * nenvs)) {
        mean_pseudo <- meanDD
      } else {
        stop("Number of values in 'meanDD' must be 1 or match number of
              environment-within-trait combinations")
      }

      if (is.null(varDD) & !(!is.null(EvarDD) & !is.null(TvarDD))) {
        stop("Either 'varDD' or 'EvarDD' and 'TvarDD' must be specified")
      } else if (!is.null(varDD)) {
        if (!is.null(EvarDD) | !is.null(TvarDD)) {
          stop("'EvarDD' and 'TvarDD' must be NULL if 'varDD' is specified")
        }
        if (length(varDD) == 1) {
          var_pseudo <- rep(varDD, each = ntraits * nenvs)
        } else if (length(varDD) == (ntraits * nenvs)) {
          var_pseudo <- varDD
        } else {
          stop("Number of values in 'varDD' must be 1 or match number of
                environment-within-trait combinations")
        }
      } else if (!is.null(EvarDD)) {
        if (length(Evar) != nenvs) {
          stop("Number of values in 'EvarDD' does not match number of
                 environments (nenvs)")
        }
        if (is.null(TvarDD)) stop("'TvarDD' must be specified")
        if (length(TvarDD) != ntraits) {
          stop("Number of values in 'TvarDD' does not must match number of
                 traits (ntraits)")
        }
        var_pseudo <- rep(TvarDD, each = nenvs) * rep(EvarDD, ntraits)
      } else {
        stop("Either 'varDD' or 'EvarDD' must be specified")
      }


      if (is.null(corDD) & !(!is.null(EcorDD) & !is.null(TcorDD))) {
        stop("Either 'corDD' or 'EcorDD' and 'TcorDD' must be specified")
      } else if (!is.null(corDD)) {
        if (!is.null(EcorDD) | !is.null(TcorDD)) {
          stop("'EcorDD' and 'TcorDD' must be NULL if 'corDD' is specified")
        }
        if (nrow(corDD) != length(mean_pseudo) | ncol(corDD) != length(mean_pseudo)) {
          stop("Dimension of 'corDD' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(corDD)) != 1) | any(corDD > 1) | any(corDD < (-1))) {
          stop("'corDD' is not a correlation matrix")
        }
        if (!isSymmetric(corDD)) stop("corDD is not symmetric")
        cor_pseudo <- corDD
      } else if (!is.null(EcorDD)) {
        if (is.null(TcorDD) & ntraits == 1) {
          TcorDD <- matrix(1)
        } else if (is.null(TcorDD) & ntraits > 1) {
          stop("'TcorDD' must be specified in combination with
                   'EcorDD'  if 'ntraits' > 1")
        } else if (length(TcorDD) == 1) {
          TcorDD <- matrix(1)
        } else {
          if (any(unique(diag(TcorDD)) != 1) | any(TcorDD > 1) | any(TcorDD < (-1))) {
            stop("'TcorDD' is not a correlation matrix")
          }
          if (!isSymmetric(TcorDD)) stop("'TcorDD' is not symmetric")
        }
        if ((nrow(EcorDD) * nrow(TcorDD)) != length(mean_pseudo) | (ncol(EcorDD) * ncol(TcorDD)) != length(mean_pseudo)) {
          stop("Dimension of the 'EcorDD' x 'TcorDD' Kronecker product does
                 not match number of environment-within-trait combinations")
        }
        if (any(unique(diag(EcorDD)) != 1) | any(EcorDD > 1) | any(EcorDD < (-1))) {
          stop("'EcorDD' is not a correlation matrix")
        }
        if (!isSymmetric(EcorDD)) stop("'EcorDD' is not symmetric")
        cor_pseudo <- kronecker(TcorDD, EcorDD)
      } else {
        stop("Either 'corDD' or 'EcorDD' must be specified")
      }


      input_asr <- c(
        input_asr,
        list(
          meanDD = mean_pseudo,
          varDD = var_pseudo,
          corDD = cor_pseudo
        )
      )
    }

    if (i == "AA") {
      if (is.null(relAA) & !(!is.null(ErelAA) & !is.null(TrelAA))) {
        stop("Either 'relAA' or 'ErelAA' and 'TrelAA' must be specified")
      } else if (!is.null(relAA)) {
        if (!is.null(ErelAA) | !is.null(TrelAA)) {
          stop("'ErelAA' and 'TrelAA' must be NULL if 'relAA' is specified")
        }
        if (length(relAA) == ntraits) {
          var_pseudo <- rep(relAA, each = nenvs)
        } else if (length(relAA) == (ntraits * nenvs)) {
          var_pseudo <- relAA
        } else {
          stop("Number of values in argument 'relAA' must match number of
               environment-within-trait combinations")
        }
      } else if (!is.null(ErelAA)) {
        if (length(ErelAA) != nenvs) {
          stop("Number of values in argument 'ErelAA' must match number of
                 environments (nenvs)")
        }
        if (is.null(TrelAA)) stop("Argument 'TrelAA' must be specified")
        if (length(TrelAA) != ntraits) {
          stop("Number of values in argument 'TrelAA' must match number of
                 traits (ntraits)")
        }
        var_pseudo <- rep(TrelAA, each = nenvs) * rep(ErelAA, ntraits)
      } else {
        stop("Either argument 'relAA' or 'ErelAA' must be specified")
      }

      if (is.null(corAA) & !(!is.null(EcorAA) & !is.null(TcorAA))) {
        stop("Either 'corAA' or 'EcorAA' and 'TcorAA' must be specified")
      } else if (!is.null(corAA)) {
        if (!is.null(EcorAA) | !is.null(TcorAA)) {
          stop("'EcorAA' and 'TcorAA' must be NULL if 'corAA' is specified")
        }
        if (nrow(corAA) != length(mean_pseudo) | ncol(corAA) != length(mean_pseudo)) {
          stop("Dimension of 'corAA' does not match number of environment-within-trait
               combinations")
        }
        if (any(unique(diag(corAA)) != 1) | any(corAA > 1) | any(corAA < (-1))) {
          stop("'corAA' is not a correlation matrix")
        }
        if (!isSymmetric(corAA)) stop("corAA is not symmetric")
        cor_pseudo <- corAA
      } else if (!is.null(EcorAA)) {
        if (is.null(TcorAA) & ntraits == 1) {
          TcorAA <- matrix(1)
        } else if (is.null(TcorAA) & ntraits > 1) {
          stop("Argument 'TcorAA' must be specified in combination with
                   'EcorAA' if ntraits > 1")
        } else if (length(TcorAA) == 1) {
          TcorAA <- matrix(1)
        } else {
          if (any(unique(diag(TcorAA)) != 1) | any(TcorAA > 1) | any(TcorAA < (-1))) {
            stop("'TcorAA' is not a correlation matrix")
          }
          if (!isSymmetric(TcorAA)) stop("'TcorAA' is not symmetric")
        }
        if ((nrow(EcorAA) * nrow(TcorAA)) != length(mean_pseudo) | (ncol(EcorAA) * ncol(TcorAA)) != length(mean_pseudo)) {
          stop("Dimension of the 'EcorAA' x 'TcorAA' Kronecker product does not
                 match number of environment-within-trait combinations")
        }
        if (any(unique(diag(EcorAA)) != 1) | any(EcorAA > 1) | any(EcorAA < (-1))) {
          stop("'EcorAA' is not a correlation matrix")
        }
        if (!isSymmetric(EcorAA)) stop("'EcorAA' is not symmetric")
        cor_pseudo <- kronecker(TcorAA, EcorAA)
      } else {
        stop("Either argument 'corAA' or 'EcorAA' must be specified")
      }


      input_asr <- c(
        input_asr,
        list(
          meanAA = rep(0, nenvs * ntraits),
          relAA = var_pseudo,
          corAA = cor_pseudo
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
#' generated using the \link[FieldSimR]{output_asr_input} function.
#'
#' @param pop An \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'} population object
#'   (\href{https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html}{Pop-class} or
#'   \href{https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html}{HybridPop-class})
#'   generated using \link[FieldSimR]{output_asr_input}.
#' @param nenvs Number of simulated environments (same number used in \link[FieldSimR]{output_asr_input}).
#' @param ntraits Number of simulated traits (same number used in \link[FieldSimR]{output_asr_input}).
#' @param nreps A vector defining the number of replicates in each environment. If only one value
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
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same values in the two environments.
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.085, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, trait 2 x 2 environments.
#' varDD <- rep(0.2, 4) # The same value set for traits 1 and 2.
#'
#' # Additive genetic correlations between the two simulated traits.
#' TcorA <- matrix(
#'   c(
#'     1.0, 0.6,
#'     0.6, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Additive genetic correlations between the two simulated environments.
#' EcorA <- matrix(
#'   c(
#'     1.0, 0.2,
#'     0.2, 1.0
#'   ),
#'   ncol = 2
#' )
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' corDD <- diag(4) # Assuming independence between traits
#'
#' input_asr <- output_asr_input(
#'   nenvs = 2,
#'   ntraits = 2,
#'   mean = mean,
#'   var = var,
#'   TcorA = TcorA,
#'   EcorA = EcorA,
#'   meanDD = meanDD,
#'   varDD = varDD,
#'   corDD = corDD
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
#'   meanDD = input_asr$meanDD,
#'   varDD = input_asr$varDD,
#'   corA = input_asr$corA,
#'   corDD = input_asr$corDD,
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
#' nreps <- c(2, 2) # Vector containing the number of complete replicates in each environment.
#'
#' gv_df <- output_asr_output(
#'   pop = pop,
#'   nenvs = 2,
#'   ntraits = 2,
#'   nreps = nreps
#' )
#' @export
output_asr_output <- function(pop,
                             nenvs,
                             ntraits,
                             nreps) {
  if (nenvs < 2) stop("'nenvs' must be > 1")
  if (nenvs %% 1 != 0) stop("'nenvs' must be an integer")

  if ((sum(nreps < 1) > 0) | (sum(nreps %% 1 != 0) > 0)) {
    stop("'nreps' must contain positive integer values")
  }
  if (length(nreps) == 1) nreps <- rep(nreps, nenvs)

  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be an integer > 0")

  envs <- factor(rep(1:nenvs, times = length(pop@id) * nreps))
  reps <- factor(unlist(lapply(nreps, function(x) rep(1:x, each = length(pop@id)))))
  ids <- factor(as.numeric(as.character(pop@id)))

  index <- as.list(as.data.frame(t(matrix(1:(ntraits * nenvs), ncol = ntraits))))
  gv <- lapply(index, function(x) cbind(pop@gv[, x]))
  gv <- do.call(rbind, mapply(function(x, y) cbind(x[rep(1:nrow(x), y), ]), x = gv, y = as.list(nreps), SIMPLIFY = F))
  colnames(gv) <- paste0("gv.Trait", 1:ntraits)

  output_asr <- data.frame(
    env = envs,
    rep = reps,
    id = ids,
    gv
  )
  output_asr <- output_asr[order(output_asr$env, output_asr$rep, output_asr$id), ]

  return(output_asr)
}
