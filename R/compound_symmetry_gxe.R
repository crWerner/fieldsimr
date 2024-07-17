#' Simulate genetic values based on a compound symmetry model for GxE interaction - `AlphaSimR`
#' input parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} to simulate
#' genetic values in multiple environments for one or more traits based on a compound symmetry
#' model for genotype-by-environment (GxE) interaction. \cr
#' This function utilises the ability of `AlphaSimR` to simulate correlated traits.
#' The wrapper function \code{compsym_asr_input()} is used to specify the input parameters required in `AlphaSimR`.
#' After simulating the genetic values, the wrapper function \link[FieldSimR]{compsym_asr_output} can be used to
#' generate a data frame with output values.
#'
#' The compound symmetry model assumes the same genetic variance for each environment
#' and the same genetic covariance between each pair of environments. New functionality
#' is being implemented which relaxes the former assumption
#' (also see \link[FieldSimR]{unstr_asr_output}).
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
#' @param var A vector of genetic variances for each trait. \cr
#'   \strong{Note:} When \code{useVarA = TRUE} is specified in `AlphaSimR` (default), the values in
#'   \code{var} represent the additive genetic variances, otherwise they represent the
#'   total (additive + non-additive) genetic variances.
#' @param prop.main  A vector defining the proportion of main effect variance for each trait.
#'   If only one value is specified, all traits will be assigned the same proportion. \cr
#'   \strong{Note:} \code{0 < prop.main < 1}.
#' @param corA A matrix of additive genetic correlations between traits. By default, a diagonal
#'   matrix is constructed.
#' @param meanDD A vector of mean dominance degrees for each environment-within-trait combination
#'   (similar to \code{mean}). If only one value is specified, all combinations will be assigned
#'   the same mean. By default, \code{meanDD = NULL} and dominance is not simulated.
#' @param varDD A vector of dominance degree variances for each trait.
#' @param prop.mainDD A vector defining the proportion of dominance degree main effect
#'   variance for each trait (similar to \code{prop.main}).
#'   If only one value is specified, all traits will be assigned the same proportion. \cr
#'   \strong{Note:} \code{0 < prop.mainDD < 1}.
#' @param corDD A matrix of dominance degree correlations between traits (similar
#'   to \code{corA}). If not specified and dominance is simulated, a diagonal matrix is constructed.
#' @param relAA A vector defining the relative magnitude of additive-by-additive (epistatic) variance
#'   to additive genetic variance for each trait, that is in a diploid organism with
#'   allele frequency of 0.5.
#'   If only one value is specified, all traits will be assigned the same relative magnitude.
#' @param prop.mainAA A vector defining the proportion of epistatic main effect variance for each
#'   trait (similar to \code{prop.main}). If only one value is specified, all traits will be assigned the
#'   same proportion. \cr
#'   \strong{Note:} \code{0 < prop.mainAA < 1}.
#' @param corAA A matrix of epistatic correlations between traits (similar to
#'   \code{corA}). If not specified and epistasis is simulated, a diagonal matrix is constructed.
#'
#' @return A list with input parameters for `AlphaSimR`, which are used to simulate
#'   correlated genetic values based on a compound symmetry model for GxE interaction.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
#' # in two environments based on a compound symmetry model.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.08, 13) # Different values for Traits 1 and 2
#' varDD <- 0.2 # Same value for Traits 1 and 2
#'
#' # Proportion of additive and dominance degree main effect variances.
#' prop.main <- c(0.4, 0.6) # Different values for Traits 1 and 2
#' prop.mainDD <- 0.4 # Same value for Traits 1 and 2
#'
#' # Additive and dominance degree correlations between the two simulated traits.
#' corA <- matrix(c(
#'   1.0, 0.5,
#'   0.5, 1.0
#' ), ncol = 2)
#' corDD <- diag(2) # Assuming independence
#'
#' input_asr <- compsym_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
#'   mean = mean,
#'   var = var,
#'   prop.main = prop.main,
#'   corA = corA,
#'   meanDD = meanDD,
#'   varDD = varDD,
#'   prop.mainDD = prop.mainDD,
#'   corDD = corDD
#' )
#'
#' @export
compsym_asr_input <- function(ntraits = 1,
                              nenvs = 2,
                              mean = 0,
                              var = 1,
                              prop.main = 0.5,
                              corA = NULL,
                              meanDD = NULL,
                              varDD = NULL,
                              prop.mainDD = NULL,
                              corDD = NULL,
                              relAA = NULL,
                              prop.mainAA = NULL,
                              corAA = NULL) {
  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
  if (nenvs < 2 | nenvs %% 1 != 0) stop("'nenvs' must be an integer > 1")

  labels <- "A"
  if (!is.null(meanDD) | !is.null(varDD) | !is.null(prop.mainDD)) {
    labels <- c("A", "DD")
  }
  if (!is.null(relAA) | !is.null(prop.mainAA)) {
    labels <- c(labels, "AA")
  }

  for (i in labels) {
    if (i == "A") {
      if (length(mean) == 1) {
        mean_vals <- rep(mean, each = ntraits * nenvs)
      } else if (length(mean) == (ntraits * nenvs)) {
        mean_vals <- mean
      } else {
        stop("Number of values in 'mean' must be 1 or match
             number of environment-within-trait combinations")
      }

      if (length(var) == 1) {
        var <- rep(var, each = ntraits)
      } else if (length(var) != ntraits) {
        stop("Number of values in 'var' must be 1 or match number of traits")
      }

      if (length(prop.main) == ntraits) {
        prop_main <- prop.main
      } else if (length(prop.main) == 1) {
        prop_main <- rep(prop.main, ntraits)
      } else {
        stop("Number of values in 'prop.main' must be 1 or
               match number of traits")
      }

      if (any(prop.main < 0) | any(prop.main > 1)) {
        stop("'prop.main' must contain values between 0 and 1")
      }

      if (is.null(corA)) corA <- diag(ntraits)

      if (nrow(corA) != ntraits | ncol(corA) != ntraits) {
        stop("Dimensions of 'corA' must match number of traits")
      }

      if (any(unique(diag(corA)) != 1) | any(corA > 1) | any(corA < (-1)) | !isSymmetric(corA)) {
        stop("'corA' must be a symmetric correlation matrix")
      }

      main_mean <- colMeans(matrix(mean_vals, ncol = ntraits))
      vars <- var
      Tcor <- as.matrix(corA)
    }

    if (i == "DD") {
      if (is.null(meanDD)) stop("'meanDD' must be specified when simulating dominance")
      if (length(meanDD) == 1) {
        mean_vals <- rep(meanDD, each = ntraits * nenvs)
      } else if (length(meanDD) == (ntraits * nenvs)) {
        mean_vals <- meanDD
      } else {
        stop("Number of values in 'meanDD' must be 1 or match
              number of environment-within-trait combinations")
      }

      if (is.null(varDD)) stop("'varDD' must be specified when simulating dominance")
      if (length(varDD) == 1) {
        varDD <- rep(varDD, each = ntraits)
      } else if (length(varDD) != ntraits) {
        stop("Number of values in 'varDD' must be 1 or match number of traits")
      }

      if (is.null(prop.mainDD)) stop("'prop.mainDD' must be specified when simulating dominance")

      if (length(prop.mainDD) == ntraits) {
        prop_main <- prop.mainDD
      } else if (length(prop.mainDD) == 1) {
        prop_main <- rep(prop.mainDD, ntraits)
      } else {
        stop("Number of values in 'prop.mainDD' must be 1 or
               match number of traits")
      }

      if (any(prop.mainDD < 0) | any(prop.mainDD > 1)) {
        stop("'prop.mainDD' must contain values between 0 and 1")
      }

      if (is.null(corDD)) corDD <- diag(ntraits)

      if (nrow(corDD) != ntraits | ncol(corDD) != ntraits) {
        stop("Dimensions of 'corDD' must match number of traits")
      }

      if (any(unique(diag(corDD)) != 1) | any(corDD > 1) | any(corDD < (-1)) | !isSymmetric(corDD)) {
        stop("'corDD' must be a symmetric correlation matrix")
      }

      main_mean <- colMeans(matrix(mean_vals, ncol = ntraits))
      vars <- varDD
      Tcor <- as.matrix(corDD)
    }

    if (i == "AA") {
      if (is.null(relAA)) stop("'relAA' must be specified when simulating epistasis")
      if (length(relAA) == 1) {
        relAA <- rep(relAA, ntraits)
      } else if (length(relAA) != ntraits) {
        stop("Number of values in 'relAA' must be 1 or match number of traits")
      }

      if (is.null(prop.mainAA)) stop("'prop.mainAA' must be specified when simulating epistasis")

      if (length(prop.mainAA) == ntraits) {
        prop_main <- prop.mainAA
      } else if (length(prop.mainAA) == 1) {
        prop_main <- rep(prop.mainAA, ntraits)
      } else {
        stop("Number of values in 'prop.mainAA' must be 1 or
               match number of traits")
      }

      if (any(prop.mainAA < 0) | any(prop.mainAA > 1)) {
        stop("'prop.mainAA' must contain values between 0 and 1")
      }

      if (is.null(corAA)) corAA <- diag(ntraits)

      if (nrow(corAA) != ntraits | ncol(corAA) != ntraits) {
        stop("Dimensions of 'corAA' must match number of traits'")
      }

      if (any(unique(diag(corAA)) != 1) | any(corAA > 1) | any(corAA < (-1)) | !isSymmetric(corAA)) {
        stop("'corAA' must be a symmetric correlation matrix")
      }


      mean_vals <- rep(0, ntraits * nenvs)
      main_mean <- rep(0, ntraits)
      vars <- relAA
      Tcor <- as.matrix(corAA)
    }

    gxe_devs <- as.list(as.data.frame(matrix(mean_vals, ncol = ntraits) - rep(main_mean, each = nenvs)))
    main_var <- vars * prop_main
    gxe_var <- vars - main_var

    mean_pseudo <- c(mapply(function(x, y) c(x, y), x = main_mean, y = gxe_devs))
    var_pseudo <- c(mapply(function(x, y) c(x, rep(y, nenvs)), x = main_var, y = gxe_var))

    cor_pseudo <- kronecker(Tcor, (diag(nenvs + 1)))

    if (i == "A") {
      input_asr <- list(
        mean = mean_pseudo,
        var = var_pseudo,
        corA = cor_pseudo
      )
    }

    if (i == "DD") {
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
      input_asr <- c(
        input_asr,
        list(
          relAA = var_pseudo,
          corAA = cor_pseudo
        )
      )
    }
  }

  return(input_asr)
}

#' Simulate genetic values based on a compound symmetry model for GxE interaction -
#' Simulation with `AlphaSimR`
#'
#' Creates a data frame of simulated genetic values in multiple environments for one or more traits
#' based on a compound symmetry model for genotype-by-environment (GxE) interaction. The wrapper function
#' \code{compsym_asr_output()} requires an \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`}
#' population object generated with \link[FieldSimR]{compsym_asr_input}.
#'
#' @param pop An \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} population object
#'   (\href{https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html}{Pop-class} or
#'   \href{https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html}{HybridPop-class})
#'   generated with \link[FieldSimR]{compsym_asr_input}.
#' @param ntraits Number of traits specified in \link[FieldSimR]{compsym_asr_input}.
#' @param nenvs Number of environments specified in \link[FieldSimR]{compsym_asr_input}.
#' @param nreps A vector defining the number of replicates in each environment. If only one value
#'   is specified, all environments will be assigned the same number.
#' @param return.effects When \code{TRUE} (default is \code{FALSE}), a list is returned with additional
#'   entries containing the genotype main effects and GxE interaction effects for each trait.
#'
#' @return A data frame with columns 'env', genotype 'id', and 'rep', followed by the
#'   simulated genetic values for each trait. When \code{return.effects = TRUE}, a list is returned with
#'   additional entries containing the genotype main effects and GxE interaction effects for each trait.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
#' # in two environments based on a compound symmetry model.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values and mean dominance degrees.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#' meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments
#'
#' # Additive genetic variances and dominance degree variances.
#' var <- c(0.08, 13) # Different values for Traits 1 and 2
#' varDD <- 0.2 # Same value for Traits 1 and 2
#'
#' # Proportion of additive and dominance degree main effect variances.
#' prop.main <- c(0.4, 0.6) # Different values for Traits 1 and 2
#' prop.mainDD <- 0.4 # Same value for Traits 1 and 2
#'
#' # Additive and dominance degree correlations between the two simulated traits.
#' corA <- matrix(c(
#'   1.0, 0.5,
#'   0.5, 1.0
#' ), ncol = 2)
#' corDD <- diag(2) # Assuming independence
#'
#' input_asr <- compsym_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
#'   mean = mean,
#'   var = var,
#'   prop.main = prop.main,
#'   corA = corA,
#'   meanDD = meanDD,
#'   varDD = varDD,
#'   prop.mainDD = prop.mainDD,
#'   corDD = corDD
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values with 'AlphaSimR' based on a
#' # compound symmetry model.
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
#' gv_ls <- compsym_asr_output(
#'   pop = pop,
#'   ntraits = 2,
#'   nenvs = 2,
#'   nreps = 2,
#'   return.effects = TRUE
#' )
#'
#' @export
compsym_asr_output <- function(pop,
                               ntraits = 1,
                               nenvs,
                               nreps = 1,
                               return.effects = FALSE) {
  if (ntraits < 1 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
  if (nenvs < 2 | nenvs %% 1 != 0) stop("'nenvs' must be an integer > 1")

  if ((sum(nreps < 1) > 0) | (sum(nreps %% 1 != 0) > 0)) {
    stop("'nreps' must contain positive integers")
  }
  if (length(nreps) == 1) nreps <- rep(nreps, nenvs)

  envs <- factor(rep(1:nenvs, times = length(pop@id) * nreps))
  reps <- factor(unlist(lapply(nreps, function(x) rep(1:x, each = length(pop@id)))))
  if (all(!grepl('\\D', pop@id))) {
    ids <- factor(as.numeric(as.character(pop@id)))
  } else {ids <- factor(as.character(pop@id))}

  main <- as.list(as.data.frame(pop@gv[, seq(1, (ntraits + ntraits * nenvs), (nenvs + 1))]))
  int <- pop@gv[, -seq(1, (ntraits + ntraits * nenvs), nenvs + 1)]
  index <- as.list(as.data.frame(t(matrix(1:(ntraits * nenvs), ncol = ntraits))))
  int <- lapply(index, function(x) int[, x])

  gv <- lapply(int, function(x) x + do.call(cbind, main))
  gv <- do.call(rbind, mapply(function(x, y) cbind(x[rep(1:nrow(x), y), ]), x = gv, y = as.list(nreps), SIMPLIFY = F))
  colnames(gv) <- paste0("gv.Trait", 1:ntraits)

  output_asr <- data.frame(
    env = envs,
    id = ids,
    rep = reps,
    gv
  )
  output_asr <- output_asr[order(output_asr$env, output_asr$rep, output_asr$id), ]
  rownames(output_asr) <- NULL

  if (return.effects) {
    main <- mapply(cbind, list(data.frame(id = pop@id)), main = main, SIMPLIFY = F)

    int <- pop@gv[, -seq(1, (ntraits + ntraits * nenvs), nenvs + 1)]
    index_eff <- as.list(as.data.frame(t(matrix(1:(ntraits * nenvs), nrow = ntraits, byrow = TRUE))))
    int <- lapply(index_eff, function(x) int[, x])
    int <- lapply(int, function(x) {
      colnames(x) <- paste0("Env", 1:nenvs)
      x
    })

    effects_df <- mapply(cbind, main, int = int, SIMPLIFY = F)
    effects_df <- lapply(effects_df, function(x) {
      x[[1]] <- factor(as.numeric(as.character(x[[1]])))
      x <- x[order(x[[1]]), ]
    })

    list_names <- c("gv.df", paste0("Trait", 1:ntraits))
    output_asr <- c(list(output_asr), effects_df)
    names(output_asr) <- list_names
  }

  return(output_asr)
}
