#' Simulate genetic values based on a multiplicative model for GxE interaction - `AlphaSimR`
#' input parameters
#'
#' Creates a list of input parameters for
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} to simulate
#' genetic values in multiple environments for one or more traits based on a (reduced rank)
#' multiplicative model for genotype-by-environment (GxE) interaction. \cr
#' This function utilises the ability of `AlphaSimR` to simulate correlated traits.
#' The wrapper function \code{multi_asr_input()} is used to specify the input parameters required in `AlphaSimR`.
#' After simulating the genetic values, the wrapper function \link[FieldSimR]{multi_asr_output} can be used to
#' generate a data frame with output values.
#'
#' Currently supports additive traits only, but other (non-additive) traits are being implemented.
#'
#' @param ntraits Number of traits to be simulated.
#' @param nenvs Number of environments to be simulated (minimum of two).
#' @param mean A vector of mean genetic values for each environment-within-trait combination.
#'   If only one value is specified, all combinations will be assigned the same mean.
#' @param var A vector of additive genetic variances for each environment-within-trait combination.
#'   If only one value is specified, all combinations will be assigned the same variance.
#' @param corA A matrix of additive genetic correlations between environment-within-trait
#'   combinations. By default, a diagonal matrix is constructed.
#' @param nterms A scalar defining the number of multiplicative terms to be simulated. By default,
#'   the number of terms is set to the number of environment-within-trait combinations.
#'   \strong{Note:} when \code{nterms} is less than the number of environment-within-trait combinations,
#'   the values in \code{mean} will be approximated.
#'
#' @return A list with input parameters for `AlphaSimR`, which are used to simulate
#'   correlated genetic values based on a multiplicative model for GxE interaction.
#'   Covariates are also supplied for use in \link[FieldSimR]{multi_asr_output}.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive traits in two
#' # environments based on a multiplicative model with three terms.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#'
#' # Additive genetic variances.
#' var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
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
#' # Construct separable additive genetic correlation matrix.
#' corA <- kronecker(TcorA, EcorA)
#'
#' input_asr <- multi_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
#'   mean = mean,
#'   var = var,
#'   corA = corA,
#'   nterms = 3
#' )
#'
#' @export
multi_asr_input <- function(ntraits = 1,
                            nenvs = 2,
                            mean = 0,
                            var = 1,
                            corA = NULL,
                            nterms = NULL) {
  if (!ntraits > 0 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
  if (!nenvs > 1 | nenvs %% 1 != 0) stop("'nenvs' must be an integer > 1")

  if (is.null(nterms)) {
    nterms <- ceiling(ntraits * nenvs)
  }
  if (!nterms > 0 | nterms %% 1 != 0) stop("'nterms' must be a positive integer")
  if (nterms > ntraits * nenvs) stop("'nterms' must be less than or equal to the number of environment-within-trait combinations")

  if (length(mean) == 1) {
    mean <- rep(mean, each = ntraits * nenvs)
  }
  if (length(mean) != (ntraits * nenvs)) {
    stop("Number of values in 'mean' must be 1 or match number of environment-within-trait combinations")
  }

  if (length(var) == 1) {
    var <- rep(var, each = ntraits * nenvs)
  }
  if (length(var) != (ntraits * nenvs)) {
    stop("Number of values in 'var' must be 1 or match number of environment-within-trait combinations")
  }
  if(any(var < 0)) stop("All values in 'var' must be greater than or equal to 0")

  if (is.null(corA)) {
    corA <- diag(ntraits * nenvs)
  }
  if (nrow(corA) != length(mean)) {
    stop("Dimensions of 'corA' must match number of environment-within-trait combinations")
  }

  if (any(unique(diag(corA)) != 1) | any(corA > 1) | any(corA < -1) | !isSymmetric(corA)) {
    stop("'corA' must be a symmetric correlation matrix")
  }

  covA <- diag(sqrt(var)) %*% corA %*% diag(sqrt(var))
  eigen_decom <- eigen(covA)
  if (any(eigen_decom$values[1:nterms] < 0)) {
    stop("'corA' must be positive (semi)-definite with regards to 'nterms'")
  }

  rank <- sum(eigen_decom$values > 1e-8)
  if (nterms == rank) {
    covariates <- cbind(eigen_decom$vectors[, 1:nterms])
    var_pseudo <- eigen_decom$values[1:nterms]
  } else if (nterms < rank) {
    term_char <- "terms"
    if (nterms == 1) {
      term_char <- "term"
    }
    message(paste0(
      "Warning message: \n 'nterms' is less than rank of 'corA', ",
      round(100 * sum(eigen_decom$values[1:nterms]) / sum(eigen_decom$values), 2), "% of variation captured with ", nterms, " ", term_char
    ))

    covariates <- cbind(eigen_decom$vectors[, 1:nterms])
    var_pseudo <- eigen_decom$values[1:nterms]
  } else if (nterms > rank) {
    message("Warning message: \n 'nterms' is greater than rank of 'corA', some terms added")
    covariates <- cbind(eigen_decom$vectors[, 1:rank])
    covariates <- cbind(covariates, matrix(0, ncol = (nterms - rank), nrow = ntraits * nenvs))
    var_pseudo <- eigen_decom$values[1:rank]
    var_pseudo <- c(var_pseudo, rep(0, nterms - rank))
  }
  if (nterms < (ntraits * nenvs) | rank < (ntraits * nenvs)) {
    message("Warning message: \n 'nterms' and/or rank of 'corA' are less than number of environment-within-trait combinations, values in 'mean' will be approximated")
  }

  mean_pseudo <- c(solve(t(covariates) %*% covariates) %*% t(covariates) %*% mean)
  cor_pseudo <- diag(nterms)
  colnames(covariates) <- paste0("cov.Term", 1:nterms)

  input_asr <- list(
    mean = mean_pseudo,
    var = var_pseudo,
    corA = cor_pseudo,
    cov.mat = covariates
  )

  return(input_asr)
}

#' Simulate genetic values based on a multiplicative model for GxE interaction -
#' Simulation with `AlphaSimR`
#'
#' Creates a data frame of simulated genetic values in multiple environments for one or more traits
#' based on a (reduced rank) multiplicative model for genotype-by-environment (GxE) interaction.
#' This function requires an \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`}
#' population object generated with \link[FieldSimR]{multi_asr_input}.
#'
#' @param pop An \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} population object
#'   (\href{https://gaynorr.github.io/AlphaSimR/reference/Pop-class.html}{Pop-class} or
#'   \href{https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.html}{HybridPop-class})
#'   generated with \link[FieldSimR]{multi_asr_input}.
#' @param ntraits Number of traits specified in \link[FieldSimR]{multi_asr_input}.
#' @param nenvs Number of environments specified in \link[FieldSimR]{multi_asr_input}.
#' @param nreps A vector defining the number of replicates in each environment. If only one value
#'   is specified, all environments will be assigned the same number.
#' @param cov.mat A matrix of covariates that will be used to construct the genetic values, typically generated
#'   with \link[FieldSimR]{multi_asr_input}.
#' @param return.effects When \code{TRUE} (default is \code{FALSE}), a list is returned with additional
#'   entries containing the genotype slopes for each multiplicative term.
#'
#' @return A data frame with columns 'env', 'rep', and genotype 'id', followed by the
#'   simulated genetic values for each trait. When \code{return.effects = TRUE}, a list is returned with
#'   additional entries containing the genotype slopes for each multiplicative term.
#'
#' @examples
#' # Simulate genetic values with 'AlphaSimR' for two additive traits in two
#' # environments based on a multiplicative model with three terms.
#'
#' # 1. Define the genetic architecture of the simulated traits.
#' # Mean genetic values.
#' mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
#'
#' # Additive genetic variances.
#' var <- c(0.086, 0.12, 15.1, 8.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
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
#' # Construct separable additive genetic correlation matrix
#' corA <- kronecker(TcorA, EcorA)
#'
#' input_asr <- multi_asr_input(
#'   ntraits = 2,
#'   nenvs = 2,
#'   mean = mean,
#'   var = var,
#'   corA = corA,
#'   nterms = 3
#' )
#'
#'
#' # 2. Use input_asr to simulate genetic values in 'AlphaSimR' based on a
#' # multiplicative model with three terms.
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
#' SP$addTraitA(
#'   nQtlPerChr = 20,
#'   mean = input_asr$mean,
#'   var = input_asr$var,
#'   corA = input_asr$corA
#' )
#'
#' pop <- newPop(FOUNDERPOP)
#'
#'
#' # 3. Create a data frame with simulated genetic values for the two traits in the two
#' # environments, with two replicates of each genotype.
#'
#' # The covariates are obtained from input_asr.
#'
#' gv_ls <- multi_asr_output(
#'   pop = pop,
#'   ntraits = 2,
#'   nenvs = 2,
#'   nreps = 2,
#'   cov.mat = input_asr$cov.mat,
#'   return.effects = TRUE
#' )
#'
#' @export
multi_asr_output <- function(pop,
                             ntraits,
                             nenvs,
                             nreps = 1,
                             cov.mat = NULL,
                             return.effects = FALSE) {
  if (!nenvs > 1 | nenvs %% 1 != 0) stop("'nenvs' must be an integer > 1")
  if (!ntraits > 0 | ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")

  if ((sum(nreps < 1) > 0) | (sum(nreps %% 1 != 0) > 0)) {
    stop("'nreps' must contain positive integers")
  }
  if (length(nreps) == 1) {
    nreps <- rep(nreps, nenvs)
  }

  envs <- factor(rep(1:nenvs, times = length(pop@id) * nreps))
  reps <- factor(unlist(lapply(nreps, function(x) rep(1:x, each = length(pop@id)))))
  ids <- factor(as.numeric(as.character(pop@id)))

  slopes <- pop@gv
  nterms <- ncol(slopes)
  if (is.null(cov.mat)) stop("'cov.mat' is not specified")
  rank <- ncol(cov.mat)
  if (nterms != rank) stop("Number of columns in 'cov.mat' does not match number of additive terms simulated")

  gv <- slopes %*% t(cov.mat)
  index <- as.list(as.data.frame(t(matrix(1:(ntraits * nenvs), ncol = ntraits))))
  gv <- lapply(index, function(x) cbind(gv[, x]))
  gv <- do.call(rbind, mapply(function(x, y) cbind(x[rep(1:nrow(x), y), ]), x = gv, y = as.list(nreps), SIMPLIFY = F))
  colnames(gv) <- paste0("gv.Trait", 1:ntraits)

  output_asr <- data.frame(
    env = envs,
    rep = reps,
    id = ids,
    gv
  )
  output_asr <- output_asr[order(output_asr$env, output_asr$rep, output_asr$id), ]

  if (return.effects) {
    colnames(slopes) <- paste0("slope.Term", 1:nterms)
    slopes <- data.frame(id = pop@id, slopes)

    listNames <- c("gv.df", "Terms")
    output_asr <- c(list(output_asr), list(slopes))
    names(output_asr) <- listNames
  }

  return(output_asr)
}
