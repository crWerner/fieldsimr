#' Simulate genetic values based on an unstructured model for GxE interaction - `AlphaSimR` input
#' parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} to simulate
#' genetic values in multiple environments for one or more traits based on an unstructured
#' model for genotype-by-environment (GxE) interaction. \cr
#' This function utilises the ability of `AlphaSimR` to simulate correlated traits.
#' The wrapper function \code{unstr_asr_input()} is used to specify the input parameters required in `AlphaSimR`,
#' and can handle separable and non-separable structures between traits and
#' environments (see below).
#' After simulating the genetic values, the wrapper function \link[FieldSimR]{unstr_asr_output} can be used to
#' generate a data frame with output values.
#'
#' \code{unstr_asr_input} can handle separable and non-separable structures between traits and
#' environments.
#' \itemize{
#'   \item For separable structures, provide (1) \code{Tvar} & \code{Evar}, and (2)
#'   \code{TcorA} & \code{EcorA}.
#'   \item For non-separable structures, provide (1) \code{var}, and (2) \code{corA}. \cr
#'   }
#'
#' \strong{Note:} `AlphaSimR` can simulate different biological effects (see:
#' \href{https://gaynorr.github.io/AlphaSimR/reference/SimParam.html}{SimParam}).
#' \itemize{
#'   \item For additive traits use \code{addTraitA()}.
#'   \item For additive + dominance traits use \code{addTraitAD()}.
#'   \item For additive + epistatic traits use \code{addTraitAE()}.
#'   \item For additive + dominance + epistatic traits use \code{addTraitADE()}.
#'   }
#' Check the \code{useVarA} argument of these functions when simulating non-additive traits.
#'
#' @param ntraits Number of traits to be simulated.
#' @param nenvs Number of environments to be simulated (minimum of two).
#' @param mean A vector of mean genetic values for each environment-within-trait combination.
#'   If only one value is specified, all combinations will be assigned the same mean.
#' @param var A vector of genetic variances for each environment-within-trait combination.
#'   If only one value is specified, all combinations will be assigned
#'   the same variance. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{Tvar} and \code{Evar} can be specified.
#' @param Tvar A vector of genetic variances for each trait. Must be provided in combination with
#'   \code{Evar}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param Evar A vector of genetic variances for each environment. Must be provided in
#'   combination with \code{Tvar}. \cr
#'   \strong{Alternatively}, \code{var} can be specified.
#' @param corA A matrix of additive genetic correlations between environment-within-trait
#'   combinations. By default, a diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{TcorA} and \code{EcorA} can be specified.
#' @param TcorA A matrix of additive genetic correlations between traits. Must be provided in
#'   combination with \code{EcorA}. \cr
#'   \strong{Alternatively}, \code{corA} can be specified.
#' @param EcorA A matrix of additive genetic correlations between environments.
#'   Must be provided in combination with \code{TcorA}. \cr
#'   \strong{Alternatively}, \code{corA} can be specified.
#' @param meanDD A vector of mean dominance degrees for each environment-within-trait combination
#'   (similar to \code{mean}). If only one value is specified, all
#'   combinations will be assigned the same mean. By default, \code{meanDD = NULL} and dominance
#'   is not simulated.
#' @param varDD A vector of dominance degree variances for each environment-within-trait
#'   combination (similar to \code{var}). If only one value is specified, all
#'   combinations will be assigned the same variance. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{TvarDD} and \code{EvarDD} can be specified.
#' @param TvarDD A vector of dominance degree variances for each trait (similar to \code{Tvar}).
#'   Must be provided in combination with \code{EvarDD}. \cr
#'   \strong{Alternatively}, \code{varDD} can be specified.
#' @param EvarDD A vector of dominance degree variances for each environment (similar to
#'   \code{Evar}). Must be provided in combination with \code{TvarDD}. \cr
#'   \strong{Alternatively}, \code{varDD} can be specified.
#' @param corDD A matrix of dominance degree correlations between environment-within-trait
#'   combinations (similar to \code{corA}). If not specified and dominance is simulated, a
#'   diagonal matrix is constructed. \cr
#'   \strong{Alternatively}, \code{TcorDD} and \code{EcorDD} can be specified.
#' @param TcorDD A matrix of dominance degree correlations between traits (similar to
#'   \code{TcorA}). Must be provided in combination with \code{EcorDD}. \cr
#'   \strong{Alternatively}, \code{corDD} can be specified.
#' @param EcorDD A matrix of dominance degree correlations between environments (similar to
#'   \code{EcorA}). Must be provided in combination with \code{TcorDD}. \cr
#'   \strong{Alternatively}, \code{corDD} can be specified.
#' @param relAA A vector defining the relative magnitude of additive-by-additive (epistatic) variance
#'   to additive genetic variance for each environment-within-trait combination,
#'   that is in a diploid organism with allele frequency of 0.5. If only one value is specified,
#'   all environment-within-trait combinations will be assigned the same value. By default,
#'   \code{relAA = NULL} and epistasis is not simulated. \cr
#'   \strong{Alternatively}, if a separable structure between traits and environments is desired,
#'   \code{TrelAA} and \code{ErelAA} can be specified.
#' @param corAA A matrix of epistatic correlations between environment-within-trait
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
#' @return A list with input parameters for `AlphaSimR`, which are used to simulate
#'   correlated genetic values based on an unstructured model for GxE interaction.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
#' # in two environments based on an unstructured model.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' varDD <- 0.2 # Same value for all environment-within-trait combinations
#'
#' # Additive genetic correlations between the two simulated traits.
#' TcorA <- matrix(c(
#'   1.0, 0.6,
#'   0.6, 1.0
#' ), ncol = 2)
#'
#' # Additive genetic correlations between the two simulated environments.
#' EcorA <- matrix(c(
#'   1.0, 0.2,
#'   0.2, 1.0
#' ), ncol = 2)
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' corDD <- diag(4) # Assuming independence
#'
#' input_asr <- unstr_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
#'   mean = mean,
#'   var = var,
#'   TcorA = TcorA,
#'   EcorA = EcorA,
#'   meanDD = meanDD,
#'   varDD = varDD,
#'   corDD = corDD
#' )
#'
#' @export
unstr_asr_input <- function(ntraits = 1,
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
                            corAA = NULL,
                            TcorAA = NULL,
                            EcorAA = NULL) {
if (!(is.atomic(ntraits) && length(ntraits) == 1L)) 
    stop("'ntraits' must be a scalar")
  if (!ntraits > 0 || ntraits%%1 != 0) 
    stop("'ntraits' must be a positive integer")
  if (!(is.atomic(nenvs) && length(nenvs) == 1L)) 
    stop("'nenvs' must be a scalar")
  if (!nenvs > 1 || nenvs%%1 != 0) 
    stop("'nenvs' must be an integer greater than 1")
  
  labelDD <- labelAA <- NULL
  if (!is.null(meanDD) || !is.null(varDD) || !is.null(EvarDD) || !is.null(TvarDD) ||
      !is.null(corDD) || !is.null(TcorDD) || !is.null(EcorDD)) {
      labelDD <- "DD"
  }
  if (!is.null(relAA) ||
      !is.null(corAA) || !is.null(TcorAA) || !is.null(EcorAA)) {
      labelAA <- "AA"
  }
  labels <- c("A", labelDD, labelAA)
  
  var_fn <- function (var, Evar = NULL, Tvar = NULL, ntraits, nenvs) {
            var_name <- "var_name"
            Evar_name <- "Evar_name"
            Tvar_name <- "Tvar_name"
            if (is.null(var) && is.null(Evar) && is.null(Tvar)) {
               message(paste0("'", var_name, "' has been set to 1 for all environment-within-trait combinations"))
               var <- 1
            }
            if (!is.null(var)) {
              if(!is.null(Evar) || !is.null(Tvar)) message(paste0("'", var_name, "' is set, '", Evar_name, "' and '", Tvar_name, "' will be ignored"))
              var <- round(c(var), 8)
              if (!is.vector(var)) stop (paste0("'", var_name, "' must be a vector"))
              Evar <- Tvar <- NULL
              if (length(var) == 1) {
                var <- rep(var, each = ntraits * nenvs)
              } else if (length(var) != (ntraits * nenvs)) {
                stop(paste0("Number of values in '", var_name, "' must match number of environment-within-trait combinations"))
              }
            } else {
              if (is.null(Evar)) {
                message(paste0("'", Evar_name, "' has been set to 1 for all environments"))
                Evar <- 1
              }
              Evar <- round(c(Evar), 8)
              if (!is.vector(Evar)) stop (paste0("'", Evar_name, "' must be a vector"))
              if (length(Evar) == 1) {
                Evar <- rep(Evar, each = nenvs)
              } else if (length(Evar) != nenvs) stop(paste0("Number of values in '", Evar_name, "' must be 1 or match number of environments"))
              if (is.null(Tvar)) {
                message(paste0("'", Tvar_name, "' has been set to 1 for all traits"))
                Tvar <- 1
              }
              Tvar <- round(c(Tvar), 8)
              if (!is.vector(Tvar)) stop (paste0("'", Tvar_name, "' must be a vector"))
              if (length(Tvar) == 1) {
                Tvar <- rep(Tvar, each = ntraits)
              } else if (length(Tvar) != ntraits) stop(paste0("Number of values in '", Tvar_name, "' must be 1 or match number of traits"))
              var <- rep(Tvar, each = nenvs) * rep(Evar, ntraits)
            }
            return(var)
  }
  
  cor_fn <- function (cor, Ecor = NULL, Tcor = NULL, ntraits, nenvs) {
            cor_name <- "cor_name"
            Ecor_name <- "Ecor_name"
            Tcor_name <- "Tcor_name"
            if (is.null(cor) && is.null(Ecor) && is.null(Tcor)) {
              message(paste0("'", cor_name, "' has been set to a diagonal matrix"))
              cor <- diag(1, nrow = ntraits * nenvs)
            }
            if (!is.null(cor)) {
              cor <- round(cbind(cor), 8)
              if (!is.matrix(cor)) stop (paste0("'", cor_name, "' must be a matrix"))
              if (!is.null(Ecor) || !is.null(Tcor)) message (paste0("'", cor_name, "' is set, '", Ecor_name, "' and '", Tcor_name, "' will be ignored"))
              Ecor <- Tcor <- NULL
              if (any(unique(diag(cor)) != 1) || any(abs(cor) > 1) || !isSymmetric(cor)) {
                stop(paste0("'", cor_name, "' must be a symmetric correlation matrix"))
              }
              if (ncol(cor) == 1) {
                 cor <- diag(1, nrows = ntraits * nenvs)
              } else if (ncol(cor) != (ntraits * nenvs)) {
                stop(paste0("Dimensions of '", cor_name, "' must match number of environment-within-trait combinations"))
              }
              } else {
              if (is.null(Ecor)) {
                message(paste0("'", Ecor_name, "' has been set to a diagonal matrix"))
                Ecor <- diag(1, nrows = nenvs)
              }
                Ecor <- round(cbind(Ecor), 8)
              if (!is.matrix(Ecor)) stop (paste0("'", Ecor_name, "' must be a matrix"))
              if (ncol(Ecor) == 1) {
                Ecor <- diag(1, nrows = nenvs)
              } else stop(paste0("Dimensions of '", Ecor_name, "' must match number of environments"))
              if (is.null(Tcor)) {
                message(paste0("'", Tcor_name, "' has been set to a diagonal matrix"))
                Tcor <- diag(1, nrows = ntraits)
              }
                Tcor <- round(cbind(Tcor), 8)
                if (!is.matrix(Tcor)) stop (paste0("'", Tcor_name, "' must be a matrix"))
                if (ncol(Tcor) == 1) {
                  Tcor <- diag(1, nrows = ntraits)
                } else stop(paste0("Dimensions of '", Tcor_name, "' must match number of traits"))
                cor_mat <- kronecker(TcorA, EcorA)
              }
            return(cor)
  }
  
  if ("A" %in% labels) {
      mean <- round(c(mean), 8)
      if (!is.vector(mean)) stop (paste0("'mean' must be a vector"))
      if (length(mean) == 1) {
         if (mean == 0) {
             message(paste0("'mean' has been set to 0 for all environment-within-trait combinations"))
             var <- 1
        }
        mean <- rep(mean, each = ntraits * nenvs)
      }
      else if (length(mean) != (ntraits * nenvs)) {
        stop("Number of values in 'mean' must match number of environment-within-trait combinations")
      }
      
      var <- var_fn(var = var, Tvar = Tvar, Evar = Evar, ntraits = ntraits, nenvs = nenvs)
      corA <- cor_fn(cor = corA, Tcor = TcorA, Ecor = EcorA, ntraits = ntraits, nenvs = nenvs)
      
      input_asr <- list(mean = mean, var = var, 
                        corA = corA)
    }
    if ("DD" %in% labels) {
      if (is.null(meanDD)) {
        message("'meanDD' has been set to 1 for all environment-within-trait combinations")
        meanDD <- 1
      }
      meanDD <- round(c(meanDD), 8)
      if (!is.vector(meanDD)) stop ("'meanDD' must be a vector")
      if (length(meanDD) == 1) {
        meanDD <- rep(meanDD, each = ntraits * nenvs)
      }
      else if (length(meanDD) != (ntraits * nenvs)) {
        stop("Number of values in 'meanDD' must match number of environment-within-trait combinations")
      }
      
      varDD <- var_fn(var = varDD, Tvar = TvarDD, Evar = EvarDD, ntraits = ntraits, nenvs = nenvs)
      corDD <- cor_fn(cor = corDD, Tcor = TcorDD, Ecor = EcorDD, ntraits = ntraits, nenvs = nenvs)
      
      input_asr <- c(input_asr, list(meanDD = meanDD, 
                                     varDD = varDD, corDD = corDD))
    }
    if ("AA" %in% labels) {
        if (is.null(relAA)) {
            message("'relAA' has been set to 1 for all environment-within-trait combinations")
            relAA <- 1
        }
        relAA <- round(c(relAA), 8)
        if (!is.vector(relAA)) stop ("'relAA' must be a vector")
        if (length(relAA) == 1) {
            relAA <- rep(relAA, each = ntraits * nenvs)
        }
        else if (length(relAA) != (ntraits * nenvs)) {
          stop("Number of values in 'relAA' must match number of environment-within-trait combinations")
        }
      if(any(relAA < 0)) stop ("'relAA' must contain positive values")
      
      corAA <- cor_fn(cor = corAA, Tcor = TcorAA, Ecor = EcorAA, ntraits = ntraits, nenvs = nenvs)
      input_asr <- c(input_asr, list(relAA = relAA, corAA = corAA))
  }
  return(input_asr)
}

#' Simulate genetic values based on an unstructured model for GxE interaction -
#' Simulation with `AlphaSimR`
#'
#' Creates a data frame of simulated genetic values in multiple environments for one or more traits
#' based on an unstructured model for genotype-by-environment (GxE) interaction. The wrapper function
#' \code{unstr_asr_output} requires an \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`}
#' population object generated with \link[FieldSimR]{unstr_asr_input}.
#'
#' @param pop An \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} population object
#'   (\href{https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html}{Pop-class} or
#'   \href{https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html}{HybridPop-class})
#'   generated with \link[FieldSimR]{unstr_asr_input}.
#' @param ntraits Number of simulated traits specified in \link[FieldSimR]{unstr_asr_input}.
#' @param nenvs Number of simulated environments specified in \link[FieldSimR]{unstr_asr_input}.
#' @param nreps A vector defining the number of replicates in each environment. If only one value
#'   is specified, all environments will be assigned the same number.
#'
#' @return A data frame with columns 'env', genotype 'id', and 'rep', followed by the
#'   simulated genetic values for each trait.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
#' # in two environments based on an unstructured model.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' varDD <- 0.2 # Same value for all environment-within-trait combinations
#'
#' # Additive genetic correlations between the two simulated traits.
#' TcorA <- matrix(c(
#'   1.0, 0.6,
#'   0.6, 1.0
#' ), ncol = 2)
#'
#' # Additive genetic correlations between the two simulated environments.
#' EcorA <- matrix(c(
#'   1.0, 0.2,
#'   0.2, 1.0
#' ), ncol = 2)
#'
#' # Dominance degree correlations between the four environment-within-trait combinations.
#' corDD <- diag(4) # Assuming independence
#'
#' input_asr <- unstr_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
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
#' # 2. Use input_asr to simulate genetic values with 'AlphaSimR' based on an
#' # unstructured model.
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
#'   corA = input_asr$corA,
#'   meanDD = input_asr$meanDD,
#'   varDD = input_asr$varDD,
#'   corDD = input_asr$corDD,
#'   useVarA = TRUE
#' )
#'
#' # By default, the variances in 'var' represent additive genetic variances.
#' # When useVarA = FALSE, the values represent total genetic variances.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame with simulated genetic values for the two traits in
#' # the two environments, with two replicates of each genotype.
#'
#' gv_df <- unstr_asr_output(
#'   pop = pop,
#'   ntraits = 2,
#'   nenvs = 2,
#'   nreps = 2
#' )
#'
#' @export
unstr_asr_output <- function(pop,
                             ntraits = 1,
                             nenvs,
                             nreps = 1) {
  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
  if (nenvs < 1 | nenvs %% 1 != 0) stop("'nenvs' must be a positive integer")

  if ((sum(nreps < 1) > 0) | (sum(nreps %% 1 != 0) > 0)) {
    stop("'nreps' must contain positive integers")
  }
  if (length(nreps) == 1) nreps <- rep(nreps, nenvs)

  envs <- factor(rep(1:nenvs, times = length(pop@id) * nreps))
  reps <- factor(unlist(lapply(nreps, function(x) rep(1:x, each = length(pop@id)))))
  if (all(!grepl("\\D", pop@id))) {
    ids <- factor(as.numeric(as.character(pop@id)))
  } else {
    ids <- factor(as.character(pop@id))
  }


  index <- as.list(as.data.frame(t(matrix(1:(ntraits * nenvs), ncol = ntraits))))
  gv <- lapply(index, function(x) cbind(pop@gv[, x]))
  gv <- do.call(rbind, mapply(function(x, y) cbind(x[rep(1:nrow(x), y), ]), x = gv, y = as.list(nreps), SIMPLIFY = F))
  colnames(gv) <- paste0("gv.Trait", 1:ntraits)

  output_asr <- data.frame(
    env = envs,
    id = ids,
    rep = reps,
    gv
  )
  output_asr <- output_asr[order(output_asr$env, output_asr$rep, output_asr$id), ]

  return(output_asr)
}
