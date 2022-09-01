#' Compound symmetry model for genotype-by-environment interaction using AlphaSimR -
#' input parameters
#'
#' Creates a list of simulation parameters for use with AlphaSimR to simulate
#' genetic values across multiple environments and traits based on a compound
#' symmetry model for genotype-by-environment (GxE) interaction. By default,
#' AlphaSimR does not support complex models for GxE interaction. However, its
#' functionality to simulate correlated genetic values can be utilised for this
#' purpose by providing the required variance structure. \code{compsym_asr_input}
#' constructs the required variance structure for a compound symmetry model.
#' It assumes a separable structure between traits and
#' environments.
#'
#' \strong{Note:} for additive traits use AlphaSimR's \code{addTraitA()};
#' for additive + dominance traits use \code{addTraitAD()}; for additive +
#' epistatic traits use \code{addTraitAE()}; and for additive + dominance +
#' epistatic traits use \code{addTraitADE()}. If non-additive effects are to be
#' simulated, check the \code{useVarA} argument of AlphaSimR's \code{addTrait}
#' methods.
#'
#' @param nEnvs Number of environments to be simulated. A minimum of two
#'   environments is required.
#' @param nTraits Number of traits to be simulated.
#' @param mean A vector of desired mean genetic values for each trait-by-environment
#'   combination (ordered as environments within traits). Simulated traits can
#'   have a different mean for each environment. If the length of \code{mean}
#'   corresponds to \code{nTraits}, however, then traits will be assigned the
#'   same mean for each environment.
#' @param var A vector of desired genetic variances for each trait.
#'   Simulated traits are restricted by the compound symmetry model to having
#'   the same variance for each environment (i.e. main effect + GxE interaction)
#'   and the same covariance between each pair of environments (main effect).
#'   \strong{Note:} when \code{useVarA = TRUE} is set in AlphaSimR (default) the
#'   values in \code{var} represent the desired \code{additive} genetic variances,
#'   otherwise they will represent the \code{total} (additive + dominance +
#'   epistatic) genetic variances.
#' @param relMainEffA  A vector with the desired magnitude of the additive main
#'   effect variance relative to the additive main effect + GxE interaction
#'   variance for each trait. If only one value is provided and \code{nTraits > 1},
#'   then all traits will be assigned the same value. \strong{Note:} \code{0 <
#'   relMainEffA < 1}.
#' @param corA A matrix of additive genetic correlations between more than one
#'   traits. If not defined and \code{nTraits > 1}, a diagonal matrix is assigned.
#' @param meanDD A vector of mean dominance degrees for each environment and
#'   trait combination (ordered as environments within traits), similar to
#'   \code{mean}. By default, \code{meanDD = NULL} and dominance is not simulated.
#' @param varDD A vector of dominance degree variances for each trait. Simulated
#'   traits have the same dominance degree variance for each environment and the
#'   same dominance degree covariance between each pair of environments, similar
#'   to \code{var}. By default, \code{varDD = NULL}.
#' @param relMainEffDD A vector with the desired magnitude of the dominance
#'   degree main effect variance relative to the main effect + GxE interaction
#'   variance for each trait, similar to \code{relMainEffA}. \strong{Note:}
#'   \code{0 < relMainEffDD < 1}. By default, \code{relMaiEffDD = NULL}.
#' @param corDD A matrix of dominance degree correlations between more than one
#'   traits, similar to \code{corA}. If not defined and \code{nTraits > 1}, a
#'   diagonal matrix is assigned. By default, \code{corDD = NULL}.
#' @param relAA  A vector with the magnitude of additive-by-additive (epistatic)
#'   variance relative to additive genetic variance for each trait, that is in a
#'   diploid organism with allele frequency 0.5. Simulated traits have the same
#'   epistatic variance for each environment and same epistatic covariance
#'   between each pair of environments, similar to \code{var}. By default,
#'   \code{relAA = NULL} and epistasis is not simulated.
#' @param relMainEffAA A vector with the magnitude of the epistatic main effect
#'   variance relative to the main effect + GxE interaction variance for each
#'   trait, similar to \code{relMainEffA}. \strong{Note:} \code{0 < relMainEffAA < 1}.
#'   By default, \code{relMainEffAA = NULL}.
#' @param corAA A matrix of epistatic correlations between  more than one traits,
#'   similar to \code{corA}. If not defined and \code{nTraits > 1}, a diagonal
#'   matrix is assigned. By default, \code{corAA = NULL}.
#'
#' @return A list containing input parameters for AlphaSimR, which is then used
#'   to simulate correlated genetic effects based on a compound symmetry model for
#'   GxE interaction.
#'
#' @examples
#' # Simulation of genetic values in AlphaSimR for two additive + dominance traits
#' # and three environments based on a compound symmetry model for GxE interaction.
#'
#' # 1. Assign genetic architecture of traits
#' # Mean genetic values and mean dominance degrees for trait 1 in 3 environments
#' # and trait 2 in 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # trait 1 by 3 envs, trait 2 by 3 envs.
#' meanDD <- c(0.1, 0.4) # trait 1 and trait 2, same values assigned to all 3 envs for each trait
#'
#' # Additive genetic variances (set useVarA=TRUE) and dominance degree variances for traits 1 and 2.
#' var <- c(0.2, 10)
#' varDD <- c(0.1, 0.2)
#'
#' # Relative magnitude of additive and dominance degree main effect variance for traits 1 and 2.
#' relMainEffA <- c(0.4, 0.6) # different values for traits 1 and 2
#' relMainEffDD <- 0.8 # same values used for traits 1 and 2
#'
#' # Additive and dominance degree correlations between traits 1 and 2.
#' corA <- matrix(c(1.0, 0.3, 0.3, 1.0), ncol = 2)
#' corDD <- diag(2) # assuming independence between traits
#'
#' input_asr <- compsym_asr_input(
#'   nEnvs = 3, nTraits = 2,
#'   mean = mean, var = var,
#'   relMainEffA = relMainEffA, corA = corA,
#'   meanDD = meanDD, varDD = varDD,
#'   relMainEffDD = relMainEffDD, corDD = corDD
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values in AlphaSimR based on a compound symmetry model for
#' #    GxE interaction.
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
#'   meanDD = input_asr$meanDD,
#'   varDD = input_asr$varDD,
#'   corA = input_asr$corA,
#'   corDD = input_asr$corDD,
#'   useVarA = TRUE
#' ) # Variance in var is used as additive variance.
#' # If FALSE, var = total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for each of the two traits
#' #    and three environments.
#'
#' nReps <- c(2, 3, 2) # Vector with the number of complete replicates in each environment
#'
#' trial_df <- compsym_asr_output(
#'   pop = pop, nEnvs = 3, nReps = nReps,
#'   nTraits = 2, effects = TRUE
#' )
#'
#' @export
compsym_asr_input <- function(nEnvs,
                              nTraits,
                              mean,
                              var,
                              relMainEffA,
                              corA = NULL,
                              meanDD = NULL,
                              varDD = NULL,
                              relMainEffDD = NULL,
                              corDD = NULL,
                              relAA = NULL,
                              relMainEffAA = NULL,
                              corAA = NULL) {
  if (nEnvs < 2) stop("'nEnvs' must be > 1")
  if (nEnvs %% 1 != 0) stop("'nEnvs' must be an integer")
  if (nTraits < 1 | nTraits %% 1 != 0) stop("'nTraits' must be an integer > 0")


  if (is.null(meanDD) & is.null(varDD) & is.null(relAA)) {
    labels <- "A"
  } else if (is.null(meanDD) & is.null(varDD)) {
    labels <- c("A", "AA")
  } else if (is.null(relAA)) {
    labels <- c("A", "DD")
  } else {
    labels <- c("A", "DD", "AA")
  }


  for (i in labels) {
    if (i == "A") {
      if (length(mean) == nTraits) {
        meanVals <- rep(mean, each = nEnvs)
      } else if (length(mean) == (nTraits * nEnvs)) {
        meanVals <- mean
      } else {
        stop("Number of values in argument 'mean' must either match number of
               traits or number of trait x environment combinations")
      }

      if (length(var) != nTraits) {
        stop("Number of values in argument 'var' must match number of traits")
      }

      if (length(relMainEffA) == nTraits) {
        relMainEff <- relMainEffA
      } else if (length(relMainEffA) == 1 & nTraits > 1) {
        relMainEff <- rep(relMainEffA, nTraits)
      } else {
        stop("Number of values in 'relMainEffA' has must either be 1 or must
                 match number of traits")
      }

      if (any(relMainEffA <= 0) | any(relMainEffA >= 1)) {
        stop("'relMainEffA' must contain values between 0 and 1")
      }

      if (is.null(corA)) corA <- diag(nTraits)

      if (nrow(corA) != nTraits | ncol(corA) != nTraits) {
        stop("Dimension of 'corA' does not match number of traits")
      }

      if (any(unique(diag(corA)) != 1) | any(corA > 1) | any(corA < (-1))) {
        stop("'corA' is not a correlation matrix")
      }

      if (!isSymmetric(corA)) stop("'corA' is not symmetric")


      mainMean <- colMeans(matrix(meanVals, ncol = nTraits))
      vars <- var
      Tcor <- as.matrix(corA)
    }

    if (i == "DD") {
      if (length(meanDD) == nTraits) {
        meanVals <- rep(meanDD, each = nEnvs)
      } else if (length(meanDD) == (nTraits * nEnvs)) {
        meanVals <- meanDD
      } else {
        stop("Number of values in argument 'meanDD' must either match number
            of traits or number of trait x environment combinations")
      }

      if (length(varDD) != nTraits) {
        stop("Number of values in argument 'varDD' must match number of traits")
      }

      if (is.null(relMainEffDD)) stop("'relMainEffDD' is not defined")

      if (length(relMainEffDD) == nTraits) {
        relMainEff <- relMainEffDD
      } else if (length(relMainEffDD) == 1 & nTraits > 1) {
        relMainEff <- rep(relMainEffDD, nTraits)
      } else {
        stop("Number of values in 'relMainEffDD' has must either be 1 or must
               match number of traits")
      }

      if (any(relMainEffDD <= 0) | any(relMainEffDD >= 1)) {
        stop("'relMainEffDD' must contain values between 0 and 1")
      }

      if (is.null(corDD)) corDD <- diag(nTraits)

      if (nrow(corDD) != nTraits | ncol(corDD) != nTraits) {
        stop("Dimension of 'corDD' does not match number of traits")
      }

      if (any(unique(diag(corA)) != 1) | any(corDD > 1) | any(corDD < (-1))) {
        stop("'corDD' is not a correlation matrix")
      }

      if (!isSymmetric(corDD)) stop("'corDD' is not symmetric")


      mainMean <- colMeans(matrix(meanVals, ncol = nTraits))
      vars <- varDD
      Tcor <- as.matrix(corDD)
    }

    if (i == "AA") {
      if (length(relAA) != nTraits) {
        stop("Number of values in argument 'relAA' must match number of traits")
      }

      if (is.null(relMainEffAA)) stop("'relMainEffAA' is not defined")

      if (length(relMainEffAA) == nTraits) {
        relMainEff <- relMainEffAA
      } else if (length(relMainEffAA) == 1 & nTraits > 1) {
        relMainEff <- rep(relMainEffAA, nTraits)
      } else {
        stop("Number of values in 'relMainEffAA' has must either be 1 or must
               match number of traits")
      }

      if (any(relMainEffAA <= 0) | any(relMainEffAA >= 1)) {
        stop("'relMainEffAA' must contain values between 0 and 1")
      }

      if (is.null(corAA)) corAA <- diag(nTraits)

      if (nrow(corAA) != nTraits | ncol(corAA) != nTraits) {
        stop("Dimension of 'corAA' does not match number of traits")
      }

      if (any(unique(diag(corA)) != 1) | any(corAA > 1) | any(corAA < (-1))) {
        stop("'corAA' is not a correlation matrix")
      }

      if (!isSymmetric(corAA)) stop("'corAA' is not symmetric")


      meanVals <- rep(0, nTraits * nEnvs)
      mainMean <- rep(0, nTraits)
      vars <- relAA
      Tcor <- as.matrix(corAA)
    }

    GEdevs <- as.list(as.data.frame(matrix(meanVals, ncol = nTraits) - rep(mainMean, each = nEnvs)))
    mainVar <- vars * relMainEff
    GEvar <- vars - mainVar

    meanPseudo <- c(mapply(function(x, y) c(x, y), x = mainMean, y = GEdevs))
    varPseudo <- c(mapply(function(x, y) c(x, rep(y, nEnvs)), x = mainVar, y = GEvar))

    corPseudo <- kronecker(Tcor, (diag(nEnvs + 1)))

    if (i == "A") {
      inputAsr <- list(
        mean = meanPseudo,
        var = varPseudo,
        corA = corPseudo
      )
    }

    if (i == "DD") {
      inputAsr <- c(
        inputAsr,
        list(
          meanDD = meanPseudo,
          varDD = varPseudo,
          corDD = corPseudo
        )
      )
    }

    if (i == "AA") {
      inputAsr <- c(
        inputAsr,
        list(
          meanAA = meanPseudo,
          relAA = varPseudo,
          corAA = corPseudo
        )
      )
    }
  }

  return(inputAsr)
}


#' Compound symmetry model for genotype-by-environment interaction using AlphaSimR -
#' simulated genetic values
#'
#' Creates a data frame of correlated genetic values across multiple traits
#' and multiple environments based on a compound symmetry model for
#' genotype-by-environment (GxE) interaction. This function requires an AlphaSimR
#' population object that was generated using input parameters from
#' \link[FieldSimR]{compsym_asr_input}.
#'
#' @param pop An AlphaSimR population object
#'   (\code{\link[AlphaSimR]{Pop-class}} or \code{\link[AlphaSimR]{HybridPop-class}})
#'   generated using \link[FieldSimR]{compsym_asr_input}.
#' @param nEnvs Number of simulated environments (same as used in
#'   \link[FieldSimR]{compsym_asr_input}).
#' @param nReps A vector with the number of complete replicates in each
#'   environment. If only one value is provided, then all environments will be
#'   assigned the same value.
#' @param nTraits Number of simulated traits (same as used in
#'   \link[FieldSimR]{compsym_asr_input}).
#' @param effects When TRUE, a list is returned  with additional columns
#'   containing the total (additive + dominance + epistatic) main effects and GxE
#'   interaction effects for each trait.
#'
#' @return A data-frame containing environment number, replicate number,
#'   genotype ID and simulated genetic values for each trait. When
#'   \code{effects = TRUE}, a list is returned with additional columns containing
#'   the total (additive + dominance + epistatic) main effects and GxE
#'   interaction effects for each trait.
#'
#' @examples
#' # Simulation of genetic values in AlphaSimR for two additive + dominance traits
#' # and three environments based on a compound symmetry model for GxE interaction.
#'
#' # 1. Assign genetic architecture of traits
#' # Mean genetic values and mean dominance degrees for trait 1 in 3 environments
#' # and trait 2 in 3 environments.
#' mean <- c(1, 3, 2, 80, 70, 100) # trait 1 by 3 envs, trait 2 by 3 envs.
#' meanDD <- c(0.1, 0.4) # trait 1 and trait 2, same values assigned to all 3 envs for each trait
#'
#' # Additive genetic variances (set useVarA=TRUE) and dominance degree variances for traits 1 and 2.
#' var <- c(0.2, 10)
#' varDD <- c(0.1, 0.2)
#'
#' # Relative magnitude of additive and dominance degree main effect variance for traits 1 and 2.
#' relMainEffA <- c(0.4, 0.6) # different values for traits 1 and 2
#' relMainEffDD <- 0.8 # same values used for traits 1 and 2
#'
#' # Additive and dominance degree correlations between traits 1 and 2.
#' corA <- matrix(c(1.0, 0.3, 0.3, 1.0), ncol = 2)
#' corDD <- diag(2) # assuming independence between traits
#'
#' input_asr <- compsym_asr_input(
#'   nEnvs = 3, nTraits = 2,
#'   mean = mean, var = var,
#'   relMainEffA = relMainEffA, corA = corA,
#'   meanDD = meanDD, varDD = varDD,
#'   relMainEffDD = relMainEffDD, corDD = corDD
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values in AlphaSimR based on a compound symmetry model for
#' #    GxE interaction.
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
#'   meanDD = input_asr$meanDD,
#'   varDD = input_asr$varDD,
#'   corA = input_asr$corA,
#'   corDD = input_asr$corDD,
#'   useVarA = TRUE
#' ) # Variance in var is used as additive variance.
#' # If FALSE, var = total genetic variance.
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame containing the simulated genetic values for each of the two traits
#' #    and three environments.
#'
#' nReps <- c(2, 3, 2) # Vector with the number of complete replicates in each environment
#'
#' trial_df <- compsym_asr_output(
#'   pop = pop, nEnvs = 3, nReps = nReps,
#'   nTraits = 2, effects = TRUE
#' )
#'
#' @export
compsym_asr_output <- function(pop, nEnvs, nReps, nTraits, effects = FALSE) {
  if (nEnvs < 2) stop("'nEnvs' must be > 1")
  if (nEnvs %% 1 != 0) stop("'nEnvs' must be an integer")

  if ((sum(nReps < 1) > 0) | (sum(nReps %% 1 != 0) > 0)) {
    stop("'nReps' must contain positive integer values")
  }
  if (length(nReps) == 1) nReps <- rep(nReps, nEnvs)

  if (nTraits < 1 | nTraits %% 1 != 0) stop("'nTraits' must be an integer > 0")

  envs <- rep(1:nEnvs, times = length(pop@id) * nReps)
  reps <- unlist(lapply(nReps, function(x) rep(1:x, each = length(pop@id))))

  g.Main <- as.list(as.data.frame(pop@gv[, seq(1, (nTraits + nTraits * nEnvs), (nEnvs + 1))]))
  gxe.Env <- pop@gv[, -seq(1, (nTraits + nTraits * nEnvs), nEnvs + 1)]
  index <- as.list(as.data.frame(t(matrix(1:(nTraits * nEnvs), ncol = nTraits))))
  gxe.Env <- lapply(index, function(x) gxe.Env[, x])

  gv <- lapply(gxe.Env, function(x) x + do.call(cbind, g.Main))
  gv <- do.call(rbind, mapply(function(x, y) x[rep(1:nrow(x), y), ], x = gv, y = as.list(nReps), SIMPLIFY = F))
  colnames(gv) <- paste0("Trait.", 1:nTraits)


  compsymAsr <- data.frame(
    env = envs,
    rep = reps,
    id = pop@id,
    gv = gv
  )


  if (effects) {
    g.Main <- mapply(cbind, list(data.frame(id = pop@id)), g.Main = g.Main, SIMPLIFY = F)

    gxe.Env <- pop@gv[, -seq(1, (nTraits + nTraits * nEnvs), nEnvs + 1)]
    indexEff <- as.list(as.data.frame(t(matrix(1:(nTraits * nEnvs), nrow = nTraits, byrow = TRUE))))
    gxe.Env <- lapply(indexEff, function(x) gxe.Env[, x])

    effComps <- mapply(cbind, g.Main, gxe.Env = gxe.Env, SIMPLIFY = F)

    listNames <- c("Trial.df", paste0("Trait.", 1:nTraits))
    compsymAsr <- list(compsymAsr)
    compsymAsr <- c(compsymAsr, effComps)
    names(compsymAsr) <- listNames
  }

  return(compsymAsr)
}
