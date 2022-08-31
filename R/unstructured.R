#' Unstructured model for genotype-by-environment interaction using AlphaSimR -
#' input parameters
#'
#' Creates a list of simulation parameters for use with AlphaSimR to simulate
#' genetic values across multiple environments and traits based on an
#' unstructured model for genotype-by-environment (GxE) interaction.
#' By default AlphaSimR does not support complex models for GxE interaction.
#' However, its functionality to simulate correlated genetic values can be
#' utilised for this purpose by providing the required variance structure.
#' \code{unstr_asr_input} constructs the required variance structure
#' for an unstructured model.
#'
#' \strong{NOTE:} \code{unstr_asr_input} handles non-separable and separable
#' structures between traits and environments. For the non-separable structure,
#' provide (1) \code{var} and (2) \code{corA}. For the separable structure,
#' provide (1) \code{Tvar}, \code{Evar} and (2) \code{TcorA}, \code{EcorA}.
#' Similar groups of terms are required for dominance and epistatic traits.
#'
#' \strong{NOTE:} for additive traits use AlphaSimR's \code{addTraitA()};
#' for additive + dominance traits use \code{addTraitAD()}; for additive +
#' epistatic traits use \code{addTraitAE()}; and for additive + dominance +
#' epistatic traits use \code{addTraitADE()}. If non-additive effects are required,
#' check the \code{useVarA} argument of AlphaSimR's \code{addTrait} methods.
#'
#'
#' @param nEnvs Number of environments to be simulated. A minimum of two
#'   environments is required.
#' @param nTraits Number of traits to be simulated.
#' @param mean A vector of desired mean genetic values for each trait-by-environment
#'   combination (ordered as environments within traits). Simulated traits can
#'   have a different mean for each environment. If the length of \code{mean}
#'   corresponds to \code{nTraits}, however, then traits will be assigned the
#'   same mean for each environment.
#' @param var A vector of desired genetic variances for each trait-by-environment
#'   combination (ordered as environments within traits). If the length of
#'   \code{var} corresponds to \code{nTraits}, however, then traits will be
#'   assigned the same genetic variance for each environment. Alternatively,
#'   \code{Tvar} and \code{Evar} can be provided if a separable structure between
#'   traits and environments is required. By default, \code{var = NULL}.
#' @param Tvar A vector of desired genetic variances for each trait. Must be
#'   provided in combination with \code{Evar}.
#'   Alternatively, \code{var} can be provided. By default, \code{Tvar = NULL}.
#' @param Evar A vector of desired genetic variances for each environment. Must
#'   be provided in combination with \code{Tvar}. Alternatively, \code{var} can
#'   be provided. By default, \code{Evar = NULL}.
#' @param corA A matrix of additive genetic correlations between all
#'   trait-by-environment combinations. If not defined and \code{nTraits > 1}, a
#'   diagonal matrix is assigned. Alternatively, \code{TcorA} and \code{EcorA}
#'   can be provided.
#' @param TcorA A matrix of additive genetic correlations between more than one
#'   traits. Must be provided in combination with \code{EcorA}. Alternatively,
#'   \code{corA} can be provided. By default, \code{TcorA = NULL}.
#' @param EcorA A matrix of additive genetic correlations between more than one
#'   environment. Must be provided in combination with \code{TcorA}. Alternatively,
#'   \code{corA} can be provided. By default, \code{EcorA = NULL}.
#' @param meanDD A vector of mean dominance degrees for each trait-by-environment
#'   combination (ordered as environments within traits), similar to \code{mean}.
#'   By default, \code{meanDD = NULL} and dominance is not simulated.
#' @param varDD A vector of dominance degree variances for each trait-by-environment
#'   combination (ordered as environments within traits), similar to \code{var}.
#'   Alternatively, \code{TvarDD} and \code{EvarDD} can be provided if a separable
#'   structure between traits and environments is required. By default,
#'   \code{varDD = NULL}.
#' @param TvarDD A vector of dominance degree variances for each trait, similar
#'   to \code{Tvar}. Must be provided in combination with \code{EvarDD}.
#'   Alternatively, \code{varDD} can be provided. By default, \code{TvarDD = NULL}.
#' @param EvarDD A vector of dominance degree genetic variances for each
#'   environment, similar to \code{Evar}. Must be provided in combination with
#'   \code{TvarDD}. Alternatively, \code{varDD} can be provided. By default,
#'   \code{EvarDD = NULL}.
#' @param corDD A matrix of dominance degree correlations between all
#'   trait-by-environment combinations, similar to \code{corA}. If not defined
#'   and \code{nTraits > 1}, a diagonal matrix is assigned. Alternatively,
#'   \code{TcorDD} and \code{EcorDD} can be provided. By default,
#'   \code{corDD = NULL}.
#' @param TcorDD A matrix of dominance degree correlations between more than one
#'   traits, similar to \code{TcorA}. Must be provided in combination with
#'   \code{EcorDD}. Alternatively, \code{corDD} can be provided. By default,
#'   \code{TcorDD = NULL}.
#' @param EcorDD A matrix of dominance degree correlations between more than one
#'   environment, similar to \code{EcorA}. Must be provided in combination with
#'   \code{TcorDD}. Alternatively, \code{corDD} can be provided. By default,
#'   \code{EcorDD = NULL}.
#' @param relAA A vector with the magnitude of additive-by-additive (epistatic)
#'   variance relative to additive genetic variance for each trait-by-environment
#'   combination (ordered as environments within traits), that is in a diploid
#'   organism with allele frequency 0.5. Alternatively, \code{TvarAA} and
#'   \code{EvarAA} can be provided if a separable structure between traits and
#'   environments is required. By default, \code{relAA = NULL}.
#' @param TrelAA A vector with the magnitude of epistatic variance relative to
#'   additive genetic variance for each trait, that is in a diploid organism with
#'   allele frequency 0.5. Must be provided in combination with \code{ErelAA}.
#'   Alternatively, \code{relAA} can be provided. By default, \code{TrelAA = NULL}.
#' @param ErelAA A vector with the magnitude of epistatic variance relative to
#'   additive genetic variance for each environment, that is in a diploid organism
#'   with allele frequency 0.5. Must be provided in combination with \code{TrelAA}.
#'   Alternatively, \code{relAA} can be provided. By default, \code{ErelAA = NULL}.
#' @param corAA A matrix of epistatic correlations between all trait-by-environment
#'   combinations, similar to \code{corA}. If not defined and \code{nTraits > 1},
#'   a diagonal matrix is assigned. Alternatively, \code{TcorAA} and \code{EcorAA}
#'   can be provided. By default, \code{corAA = NULL}.
#' @param TcorAA A matrix of epistatic correlations between more than one traits,
#'   similar to \code{TcorA}. Must be provided in combination with \code{EcorAA}.
#'   Alternatively, \code{corAA} can be provided. By default, \code{TcorAA = NULL}.
#' @param EcorAA A matrix of epistatic correlations between more than one
#'   environment, similar to \code{EcorA}. Must be provided in combination with
#'   \code{TcorAA}. Alternatively, \code{corAA} can be provided. By default,
#'   \code{EcorAA = NULL}.
#'
#' @return A list containing input parameters for AlphaSimR, which is then used
#'   to simulate correlated genetic effects based on an unstructured model for
#'   GxE interaction.
#'
#' @examples
#' # Simulation of genetic values for two additive + dominance traits in three
#' # environments in AlphaSimR based on an unstructured GxE interaction model.
#'
#' # 1. Assign genetic architecture of traits
#' # Mean genetic values and mean dominance degrees for trait 1 in 3 environments
#' # and trait 2 in 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # trait 1 by 3 envs, trait 2 by 3 envs.
#' meanDD <- c(0.1, 0.4) # trait 1 and trait 2, same values assigned to all 3 envs for each trait
#'
#' # Additive genetic variances (set usevarA=TRUE) and dominance degree variances for the two traits,
#' # that is assuming a separable structure between traits and environments
#' Tvar <- c(0.2, 10)
#' Evar <- c(0.5, 1, 1.5)
#'
#' # Dominance degree variances for trait 1 in 3 environments and trait 2 in 3 environments,
#' # that is assuming a non-separable structure between traits and environments
#' varDD <- c(0.1, 0.15, 0.2, 0.2, 0.3, 0.2)
#'
#'
#' # Additive genetic correlations between traits.
#' TcorA <- matrix(c(1.0,  0.3,
#'                   0.3,  1.0), ncol = 2)
#'
#' # Additive genetic correlations between environments.
#' EcorA <- cov2cor(matrix(c(0.5, 0.4,  0.6,
#'                           0.4, 1,  0.5,
#'                           0.6, 0.5,  1.5), ncol = 3))
#'
#' # Dominance degree correlation between all six trait-by-environment cobinations
#' corDD <- diag(6) # assuming independence between traits
#'
#' input_asr <- unstr_asr_input(nEnvs = 3, nTraits = 2, mean = mean,
#'                              Tvar = Tvar, Evar = Evar,
#'                              TcorA = TcorA, EcorA = EcorA,
#'                              meanDD = meanDD, varDD = varDD,
#'                              corDD = corDD)
#'
#'
#' # 2. Use input_asr to simulate genetic values in AlphaSimR based on an unstructured model for
#' #    GxE interaction.
#'
#' library(AlphaSimR)
#' FOUNDERPOP <- quickHaplo(nInd = 100,
#'                          nChr = 6,
#'                          segSites = 100)
#'
#' SP <- SimParam$new(FOUNDERPOP)
#'
#' SP$addTraitAD(nQtlPerChr = 100,
#'               mean = input_asr$mean,
#'               var = input_asr$var,
#'               meanDD = input_asr$meanDD,
#'               varDD = input_asr$varDD,
#'               corA = input_asr$corA,
#'               corDD = input_asr$corDD,
#'               useVarA = TRUE) # Variance in var is used as additive variance.
#'                               # If FALSE, var = total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for each of the two traits
#' #    and three environments.
#'
#' nReps <- c(2, 3, 2)  # Vector with the number of complete replicates in each environment
#'
#' trial_df <- unstr_asr_output(pop = pop, nEnvs = 3, nReps = nReps, nTraits = 2)
#'
#' @export
unstr_asr_input <- function(nEnvs, nTraits, mean,
                            var = NULL, Tvar = NULL, Evar = NULL,
                            corA = NULL, EcorA = NULL, TcorA = NULL,
                            meanDD = NULL, varDD = NULL,
                            EvarDD = NULL, TvarDD = NULL,
                            corDD = NULL, EcorDD = NULL, TcorDD = NULL,
                            relAA = NULL, ErelAA = NULL, TrelAA = NULL,
                            corAA = NULL, EcorAA = NULL, TcorAA = NULL) {


  if (nEnvs < 2) stop("'nEnvs' must be > 1")
  if (nEnvs %% 1 != 0) stop("'nEnvs' must be an integer")
  if (nTraits < 1 | nTraits %% 1 != 0) stop("'nTraits' must be an integer > 0")


  if (is.null(meanDD) & is.null(varDD) & is.null(EvarDD) & is.null(relAA) & is.null(ErelAA)) {
    labels <- "A"
    } else if (is.null(meanDD) & is.null(varDD) & is.null(EvarDD)) {
      labels <- c("A", "AA")
      } else if (is.null(relAA) & is.null(ErelAA)) {
        labels <- c("A", "DD")
        } else {
          labels <- c("A", "DD", "AA")
        }


  for(i in labels){

    if(i == "A"){

      if (length(mean) == nTraits) {
        meanPseudo <- rep(mean, each = nEnvs)
        } else if (length(mean) == (nTraits * nEnvs)) {
          meanPseudo <- mean
          } else {
            stop("Number of values in argument 'mean' must either match number of
               traits or number of trait x environment combinations")
          }


      if (!is.null(var) & !is.null(Evar)) {
        stop("Either argument 'var' or 'Evar' must be provided")
        } else if (!is.null(var)) {
          if (!is.null(Tvar)) {
            stop("Argument 'Tvar' must be NULL if 'var' is provided")
            }
          if (length(var) == (nTraits * nEnvs)) {
            varPseudo <- var
            } else {
            stop("Number of values in argument 'var' must match number of
                 trait x environment combinations")
              }
          } else if (!is.null(Evar)) {
            if (length(Evar) != nEnvs) {
            stop("Number of values in argument 'Evar' must match number of
                 environments (nEnvs)")
              }
            if (is.null(Tvar)) stop("Argument 'Tvar' is not provided")
            if (length(Tvar) != nTraits) {
            stop("Number of values in argument 'Tvar' must match number of
                 traits (nTraits)")
              }
            varPseudo <- rep(Tvar, each = nEnvs) * rep(Evar, nTraits)
            } else {
              stop("Either argument 'var' or 'Evar' must be provided")
            }


      if (!is.null(corA) & !is.null(EcorA)) { ### CHECK @Chrstian
        stop("Either argument 'corA' or 'EcorA' must be provided")
        } else if (!is.null(corA)) {
          if (!is.null(TcorA)) {
            stop("Argument 'TcorA' must be NULL if 'corA' is provided")
            }
          if (nrow(corA) != length(meanPseudo) | ncol(corA) != length(meanPseudo)) {
          stop("Dimension of 'corA' does not match number of trait x
               environment combinations")
            }
          if (any(unique(diag(corA)) != 1) | any(corA > 1) | any(corA < (-1))) {
            stop("'corA' is not a correlation matrix")
            }
          if (!isSymmetric(corA)) stop("corA is not symmetric")
          corPseudo <- corA
          } else if (!is.null(EcorA)) {
            if (is.null(TcorA) & nTraits == 1) {
              TcorA <- matrix(1)
              } else if (is.null(TcorA) & nTraits > 1) {
              stop("Argument 'TcorA' must be provided if in combination with
                   'EcorA' if nTraits > 1")
                } else if(length(TcorA) == 1) {
                  TcorA <- matrix(1)
                  } else {
                    if(any(unique(diag(TcorA)) != 1) | any(TcorA > 1) | any(TcorA < (-1))) {
                      stop("'TcorA' is not a correlation matrix")
                      }
                    if(!isSymmetric(TcorA)) stop("'TcorA' is not symmetric")
                    }
            if ((nrow(EcorA) * nrow(TcorA)) != length(meanPseudo) | (ncol(EcorA) * ncol(TcorA)) != length(meanPseudo)) {
            stop("Dimension of the 'EcorA' x 'TcorA' Kronecker product does not
                 match number of trait x environment combinations")
              }
            if (any(unique(diag(EcorA)) != 1) | any(EcorA > 1) | any(EcorA < (-1))) {
              stop("'EcorA' is not a correlation matrix")
              }
            if (!isSymmetric(EcorA)) stop("'EcorA' is not symmetric")
            corPseudo <- kronecker(TcorA, EcorA)
            } else {
              stop("Either argument 'corA' or 'EcorA' must be provided")
              }


      inputAsr <- list(mean = meanPseudo,
                       var = varPseudo,
                       corA = corPseudo)

    }

    if(i == "DD"){

      if (length(meanDD) == nTraits) {
        meanPseudo <- rep(meanDD, each = nEnvs)
        } else if (length(meanDD) == (nTraits * nEnvs)) {
          meanPseudo <- meanDD
          } else {
          stop("Number of values in argument 'meanDD' must either match number
               of traits or number of trait x environment combinations")
          }


      if (!is.null(varDD) & !is.null(EvarDD)) {
        stop("Either argument 'varDD' or 'EvarDD' must be provided")
        } else if (!is.null(varDD)) {
          if (!is.null(TvarDD)) {
            stop("Argument 'TvarDD' must be NULL if 'varDD' is provided")
            }
          if (length(varDD) == (nTraits * nEnvs)) {
            varPseudo <- varDD
            } else {
            stop("Number of values in argument 'varDD' must match number of
                 trait x environment combinations")
              }
          } else if (!is.null(EvarDD)) {
            if (length(Evar) != nEnvs) {
            stop("Number of values in argument 'EvarDD' must match number of
                 environments (nEnvs)")
              }
            if (is.null(TvarDD)) stop("Argument 'TvarDD' is not provided")
            if (length(TvarDD) != nTraits) {
            stop("Number of values in argument 'TvarDD' must match number of
                 traits (nTraits)")
              }
            varPseudo <- rep(TvarDD, each = nEnvs) * rep(EvarDD, nTraits)
            } else {
              stop("Either argument 'varDD' or 'EvarDD' must be provided")
            }


      if (!is.null(corDD) & !is.null(EcorDD)) {
        stop("Either argument 'corDD' or 'EcorDD' must be provided")
        } else if (!is.null(corDD)) {
          if (!is.null(TcorDD)) {
            stop("Argument 'TcorDD' must be NULL if 'corDD' is provided")
            }
          if (nrow(corDD) != length(meanPseudo) | ncol(corDD) != length(meanPseudo)) {
          stop("Dimension of 'corDD' does not match number of trait x
               environment combinations")
            }
          if (any(unique(diag(corDD)) != 1) | any(corDD > 1) | any(corDD < (-1))) {
            stop("'corDD' is not a correlation matrix")
            }
          if (!isSymmetric(corDD)) stop("corDD is not symmetric")
          corPseudo <- corDD
          } else if (!is.null(EcorDD)) {
            if (is.null(TcorDD) & nTraits == 1) {
              TcorDD <- matrix(1)
              } else if (is.null(TcorDD) & nTraits > 1) {
              stop("Argument 'TcorDD' must be provided if in combination with
                   'EcorDD' if nTraits > 1")
                } else if(length(TcorDD) == 1) {
                  TcorDD <- matrix(1)
                  } else {
                    if(any(unique(diag(TcorDD)) != 1) | any(TcorDD > 1) | any(TcorDD < (-1))) {
                      stop("'TcorDD' is not a correlation matrix")
                      }
                    if(!isSymmetric(TcorDD)) stop("'TcorDD' is not symmetric")
                    }
            if ((nrow(EcorDD) * nrow(TcorDD)) != length(meanPseudo) | (ncol(EcorDD) * ncol(TcorDD)) != length(meanPseudo)) {
            stop("Dimension of the 'EcorDD' x 'TcorDD' Kronecker product does
                 not match number of trait x environment combinations")
              }
            if (any(unique(diag(EcorDD)) != 1) | any(EcorDD > 1) | any(EcorDD < (-1))) {
              stop("'EcorDD' is not a correlation matrix")
              }
            if (!isSymmetric(EcorDD)) stop("'EcorDD' is not symmetric")
            corPseudo <- kronecker(TcorDD, EcorDD)
            } else {
              stop("Either argument 'corDD' or 'EcorDD' must be provided")
            }


      inputAsr <- c(inputAsr,
                    list(meanDD = meanPseudo,
                         varDD = varPseudo,
                         corDD = corPseudo))

    }

    if(i == "AA"){

      if (!is.null(relAA) & !is.null(ErelAA)) {
        stop("Either argument 'relAA' or 'ErelAA' must be provided")
        } else if (!is.null(relAA)) {
          if (!is.null(TrelAA)) {
            stop("Argument 'TrelAA' must be NULL if 'relAA' is provided")
            }
          if (length(relAA) == (nTraits * nEnvs)) {
          varPseudo <- relAA
          } else {
          stop("Number of values in argument 'relAA' must match number of
               trait x environment combinations")
            }
          } else if (!is.null(ErelAA)) {
            if (length(ErelAA) != nEnvs) {
            stop("Number of values in argument 'ErelAA' must match number of
                 environments (nEnvs)")
              }
            if (is.null(TrelAA)) stop("Argument 'TrelAA' is not provided")
            if (length(TrelAA) != nTraits) {
            stop("Number of values in argument 'TrelAA' must match number of
                 traits (nTraits)")
              }
            varPseudo <- rep(TrelAA, each = nEnvs) * rep(ErelAA, nTraits)
            } else {
              stop("Either argument 'relAA' or 'ErelAA' must be provided")
            }


      if (!is.null(corAA) & !is.null(EcorAA)) {
        stop("Either argument 'corAA' or 'EcorAA' must be provided")
        } else if (!is.null(corAA)) {
          if (!is.null(TcorAA)) {
            stop("Argument 'TcorAA' must be NULL if 'corAA' is provided")
            }
          if (nrow(corAA) != length(meanPseudo) | ncol(corAA) != length(meanPseudo)) {
          stop("Dimension of 'corAA' does not match number of trait x
               environment combinations")
            }
          if (any(unique(diag(corAA)) != 1) | any(corAA > 1) | any(corAA < (-1))) {
            stop("'corAA' is not a correlation matrix")
            }
          if (!isSymmetric(corAA)) stop("corAA is not symmetric")
          corPseudo <- corAA
          } else if (!is.null(EcorAA)) {
            if (is.null(TcorAA) & nTraits == 1) {
              TcorAA <- matrix(1)
              } else if (is.null(TcorAA) & nTraits > 1) {
              stop("Argument 'TcorAA' must be provided if in combination with
                   'EcorAA' if nTraits > 1")
                } else if(length(TcorAA) == 1) {
                  TcorAA <- matrix(1)
                  } else {
                    if(any(unique(diag(TcorAA)) != 1) | any(TcorAA > 1) | any(TcorAA < (-1))) {
                      stop("'TcorAA' is not a correlation matrix")
                      }
                    if(!isSymmetric(TcorAA)) stop("'TcorAA' is not symmetric")
                    }
            if ((nrow(EcorAA) * nrow(TcorAA)) != length(meanPseudo) | (ncol(EcorAA) * ncol(TcorAA)) != length(meanPseudo)) {
            stop("Dimension of the 'EcorAA' x 'TcorAA' Kronecker product does not
                 match number of trait x environment combinations")
              }
            if (any(unique(diag(EcorAA)) != 1) | any(EcorAA > 1) | any(EcorAA < (-1))) {
              stop("'EcorAA' is not a correlation matrix")
              }
            if (!isSymmetric(EcorAA)) stop("'EcorAA' is not symmetric")
            corPseudo <- kronecker(TcorAA, EcorAA)
            } else {
              stop("Either argument 'corAA' or 'EcorAA' must be provided")
            }


      inputAsr <- c(inputAsr,
                    list(meanAA = rep(0, nEnvs * nTraits),
                         relAA = varPseudo,
                         corAA = corPseudo))
    }
  }

  return(inputAsr)

}



#' Unstructured genotype-by-environment interaction model in AlphaSimR -
#' simulated genetic values
#'
#' Creates a data frame of correlated genetic values across multiple traits
#' and multiple environments based on an unstructured model for
#' genotype-by-environment (GxE) interaction. This function requires an AlphaSimR
#' population object that was generated using input parameters from
#' \link[FieldSimR]{unstr_asr_input}.
#'
#' @param pop An AlphaSimR population object
#'   (\code{\link[AlphaSimR]{Pop-class}} or \code{\link[AlphaSimR]{HybridPop-class}})
#'    generated using \link[FieldSimR]{unstr_asr_input}.
#' @param nEnvs Number of simulated environments (same as in \link[FieldSimR]{unstr_asr_input}.
#' @param nReps A vector with the number of complete replicates in each
#'   environment. If only one value is provided, then all environments will be assigned the same value.
#' @param nTraits Number of simulated traits (same as in \link[FieldSimR]{unstr_asr_input}.
#'
#' @return A data-frame containing environment number, replicate number, genotype
#'   ID and simulated genetic values for each trait.
#'
#' @examples
#' # Simulation of genetic values for two additive + dominance traits in three
#' # environments in AlphaSimR based on an unstructured GxE interaction model.
#'
#' # 1. Assign genetic architecture of traits
#' # Mean genetic values and mean dominance degrees for trait 1 in 3 environments
#' # and trait 2 in 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # trait 1 by 3 envs, trait 2 by 3 envs.
#' meanDD <- c(0.1, 0.4) # trait 1 and trait 2, same values assigned to all 3 envs for each trait
#'
#' # Additive genetic variances (set usevarA=TRUE) and dominance degree variances for the two traits,
#' # that is assuming a separable structure between traits and environments
#' Tvar <- c(0.2, 10)
#' Evar <- c(0.5, 1, 1.5)
#'
#' # Dominance degree variances for trait 1 in 3 environments and trait 2 in 3 environments,
#' # that is assuming a non-separable structure between traits and environments
#' varDD <- c(0.1, 0.15, 0.2, 0.2, 0.3, 0.2)
#'
#'
#' # Additive genetic correlations between traits.
#' TcorA <- matrix(c(1.0,  0.3,
#'                   0.3,  1.0), ncol = 2)
#'
#' # Additive genetic correlations between environments.
#' EcorA <- cov2cor(matrix(c(0.5, 0.4,  0.6,
#'                           0.4, 1,  0.5,
#'                           0.6, 0.5,  1.5), ncol = 3))
#'
#' # Dominance degree correlation between all six trait-by-environment cobinations
#' corDD <- diag(6) # assuming independence between traits
#'
#' input_asr <- unstr_asr_input(nEnvs = 3, nTraits = 2, mean = mean,
#'                              Tvar = Tvar, Evar = Evar,
#'                              TcorA = TcorA, EcorA = EcorA,
#'                              meanDD = meanDD, varDD = varDD,
#'                              corDD = corDD)
#'
#'
#' # 2. Use input_asr to simulate genetic values in AlphaSimR based on an unstructured model for
#' #    GxE interaction.
#'
#' library(AlphaSimR)
#' FOUNDERPOP <- quickHaplo(nInd = 100,
#'                          nChr = 6,
#'                          segSites = 100)
#'
#' SP <- SimParam$new(FOUNDERPOP)
#'
#' SP$addTraitAD(nQtlPerChr = 100,
#'               mean = input_asr$mean,
#'               var = input_asr$var,
#'               meanDD = input_asr$meanDD,
#'               varDD = input_asr$varDD,
#'               corA = input_asr$corA,
#'               corDD = input_asr$corDD,
#'               useVarA = TRUE) # Variance in var is used as additive variance.
#'                               # If FALSE, var = total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for each of the two traits
#' #    and three environments.
#'
#' nReps <- c(2, 3, 2)  # Vector with the number of complete replicates in each environment
#'
#' trial_df <- unstr_asr_output(pop = pop, nEnvs = 3, nReps = nReps, nTraits = 2)
#'
#' @export
unstr_asr_output <- function(pop, nEnvs, nReps, nTraits) {

  if (nEnvs < 2) stop("'nEnvs' must be > 1")
  if (nEnvs %% 1 != 0) stop("'nEnvs' must be an integer")

  if ((sum(nReps < 1) > 0) | (sum(nReps %% 1 != 0) > 0)) {
    stop("'nReps' must contain positive integer values")
  }
  if (length(nReps) == 1) nReps <- rep(nReps, nEnvs)

  if (nTraits < 1 | nTraits %% 1 != 0) stop("'nTraits' must be an integer > 0")

  envs <- rep(1:nEnvs, times = length(pop@id) * nReps)
  reps <- unlist(lapply(nReps, function(x) rep(1:x, each = length(pop@id))))

  index <- as.list(as.data.frame(t(matrix(1:(nTraits * nEnvs), ncol = nTraits))))
  gv <- lapply(index, function(x) pop@gv[ , x])
  gv <- do.call(rbind, mapply(function(x, y) x[rep(1:nrow(x), y), ], x = gv, y = as.list(nReps), SIMPLIFY = F))
  colnames(gv) <- paste0("Trait.", 1:nTraits)

  unstrAsr <- data.frame(env = envs,
                         rep = reps,
                         id = pop@id,
                         gv = gv)

 return(unstrAsr)

}
