#' Phenotype simulation through combination of genetic values and plot-level errors
#'
#' Creates a data frame of simulated field trial phenotypes through combination of the genetic
#' values and the plot-level errors generated for one or more traits with 'FieldSimR'. If the
#' genetic values were obtained externally, they have to be arranged in a data frame with columns
#' "env", "rep", and "id" additional to the genetic values for each trait.
#'
#' @param gv_df A data frame of genetic values. Must contain the columns "env", "rep", and "id"
#'   additional to the genetic values for each trait.
#' @param error_df A data frame of plot-level errors. Must contain the columns "env", "block",
#'   "col", and "row" additional to the error values for each trait.
#' @param randomise When TRUE, genotypes are randomly allocated to plots within blocks to simulate
#'   a randomized complete block design (RCBD).\cr
#'   \strong{Note:} other experimental designs must be generated externally.
#'
#' @return A data-frame containing the environment id, block id, column id, row id, genotype id,
#'   and the phenotypic values for each trait.
#'
#' @examples
#' # TBA
#'
#' @export
make_phenotypes <- function(gv_df, error_df, randomise = FALSE) {
  if (inherits(gv_df, "list")) gv_df <- gv_df[[1]]
  if (inherits(error_df, "list")) error_df <- error_df[[1]]

  colnames(gv_df) <- toupper(colnames(gv_df))
  colnames(error_df) <- toupper(colnames(error_df))

  if (any(!c("ENV", "REP", "ID") %in% colnames(gv_df))) {
    stop("'gv_df' must contain the columns 'env', 'rep', 'id', and a column with genetic values for each simulated trait.")
  }

  if (any(!c("ENV", "BLOCK", "COL", "ROW") %in% colnames(error_df))) {
    stop("'error_df' must contain the columns 'env', 'block', 'col', 'row', and a column with error values for each simulated trait.")
  }

  n_traits <- ncol(gv_df) - 3
  if (ncol(error_df) - 4 != n_traits) {
    stop("'gv_df' must contain the columns 'env', 'rep', 'id', and a column with genetic values for each simulated trait.
         'error_df' must contain the columns 'env', 'block', 'col', 'row', and a column with error values for each simulated trait.")
  }

  error_df <- error_df[order(error_df$ENV, error_df$BLOCK), ]
  gv_df <- gv_df[order(gv_df$ENV, gv_df$REP), ]

  phe <- error_df[, !(colnames(error_df) %in% c("ENV", "BLOCK", "COL", "ROW"))] +
    gv_df[, !(colnames(gv_df) %in% c("ENV", "REP", "ID"))]

  error_df <- error_df[, (colnames(error_df) %in% c("ENV", "BLOCK", "COL", "ROW"))]
  error_df <- error_df[, c("ENV", "BLOCK", "COL", "ROW")]

  if (randomise) {
    id <- unique(gv_df$ID)
    n_blocks_total <- nrow(gv_df) / length(id)

    id_rand <- sample(id)
    for (i in 2:n_blocks_total) {
      id_rand <- c(id_rand, sample(id))
    }

    pheno_df <- cbind(
      error_df,
      id_rand,
      phe
    )
  } else {
    pheno_df <- cbind(
      error_df,
      gv_df$ID,
      phe
    )
  }

  colnames(pheno_df) <- c("env", "block", "col", "row", "id", paste0("phe.Trait.", 1:n_traits))
  return(pheno_df)
}
