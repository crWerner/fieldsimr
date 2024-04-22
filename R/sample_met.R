#' Sample environments from a target population
#'
#' Creates a list of environments sampled from a population with user-defined sample size.
#'
#' @param ntraits A scalar defining the number of traits.
#' @param nenvs A scalar defining the number of environments in the target population.
#' @param nsamples A scalar defining the number of samples to be taken.
#' @param sample.size A vector defining the number of environments in each sample.
#'   When only one value is specified, all samples will be assigned the same number.
#' @param replace When \code{TRUE} (default), samples are taken with replacement.
#'   Ignored when \code{nsamples = 1}.
#' @param cov.mat An optional matrix of environmental covariates for one or more traits.
#'   When supplied, the covariates are sampled and printed.
#'
#' @return A list with elements given by the sample of environments taken from the target population.
#'   When \code{cov.mat} is supplied, additional entries are given containing the sampled environmental
#'   covariates for each trait.
#'
#' @examples
#' # Sample environments from a target population of 1000, with each sample containing 20 environments.
#' cov_ls <- sample_met(
#'   nenvs = 1000,
#'   nsamples = 10,
#'   sample.size = 20,
#'   replace = TRUE
#' )
#'
#' @export
sample_met <- function(ntraits = 1,
                       nenvs = 1000,
                       nsamples = 10,
                       sample.size = 20,
                       replace = TRUE,
                       cov.mat = NULL) {
  if (!(is.atomic(nsamples) && length(nsamples) == 1L)) stop("'nsamples' must be a scalar")
  if (nsamples < 1 || nsamples %% 1 != 0) stop("'nsamples' must be a positive integer")

  if(!is.vector(sample.size)) stop("'sample.size' must be a vector or scalar")
  if (length(sample.size) == 1) {
    sample.size <- rep(sample.size, nsamples)
  }
  if (nsamples != length(sample.size)) stop("'nsamples' must match length of 'sample.size'")
  if (any(sample.size < 1) || any(sample.size %% 1 != 0)) stop("All values in 'sample.size' must be positive integers")

  if (!(is.atomic(nenvs) && length(nenvs) == 1L)) stop("'nenvs' must be a scalar")
  if (!nenvs > 1 | nenvs %% 1 != 0) stop("'nenvs' must be an integer > 1")

  if (replace) {
    if (any(sample.size > nenvs)) stop("All values in 'sample.size' must be less than or equal to 'nenvs'")
    sample_ls <- lapply(sample.size, function(x) sample(1:nenvs, x))
  } else if (!replace) {
    if (sum(sample.size) > nenvs) stop("Sum of 'sample.size' across all samples must be less than or equal to 'nenvs' when 'replace = FALSE'")
    cumsum_sample_size <- cumsum(sample.size)
    if (nsamples == 1) {
      sample_table <- matrix(c(1, cumsum_sample_size), ncol = 2)
    } else if (nsamples > 1) {
      sample_table <- matrix(c(1, cumsum_sample_size[1:(nsamples-1)] + 1, cumsum_sample_size), ncol = 2)
    }
    sample_ls <- sample(1:nenvs)
    sample_ls <- mapply(function(x,y) sample_ls[x:y], x = sample_table[,1], y = sample_table[,2], SIMPLIFY = FALSE)
  }

  if (!is.null(cov.mat)) {
    cov.mat <- cbind(cov.mat)
    if (!is.matrix(cov.mat)) stop("'cov.mat' must be a matrix")
    if (!(is.atomic(ntraits) && length(ntraits) == 1L)) stop("'ntraits' must be a scalar")
    if (ntraits < 1 || ntraits %% 1 != 0) stop("'ntraits' must be a positive integer")
    if (nrow(cov.mat) != (ntraits * nenvs)) {
      stop("Number of rows in 'cov.mat' must match number of environment-within-trait combinations")
    }
    nterms <- ncol(cov.mat)
    trait_grid <- lapply(sample.size, function(x) rep(nenvs*(seq(0,ntraits-1)), each = x))
    cov_ls <- mapply(function(x,y) cbind(cov.mat[x + y,]), x = sample_ls, y = trait_grid, SIMPLIFY = FALSE)
    cov_ls <- lapply(cov_ls, function(x) {
      colnames(x) <- paste0("cov.Term", 1:nterms)
      x
    })
    cov_ls <- list(sample = sample_ls, cov.mat = cov_ls)
  } else if (is.null(cov.mat)) {
    cov_ls <- sample_ls
  }
return(cov_ls)
}
