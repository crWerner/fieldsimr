#' Simulate phenotypes - Combine simulated genetic values and plot errors
#'
#' Creates a data frame of simulated phenotypes for one or more traits by
#' combining simulated plot errors with genetic values (e.g. true, simulated or predicted).
#' The genetic values can be generated externally, but note that they must be stored in a
#' data frame as described below.
#'
#' @param gv_df A data frame of genetic values. Must contain the columns 'env', 'rep', and 'id',
#'   followed by the genetic values for each trait.
#' @param error_df A data frame of plot errors. Must contain the columns 'env', 'block',
#'   'col', and 'row', followed by the plot errors for each trait.
#' @param randomise When TRUE, genotypes are randomly allocated to plots within blocks to generate
#'   a randomized complete block design (RCBD).\cr
#'   \strong{Note:} Other experimental designs must be generated externally.
#'
#' @return A data frame with columns 'env', 'block', 'column', 'row' and 'genotype', followed by
#'   the phenotypes for each trait.
#'
#' @examples
#' # Simulate phenotypes by combining the genetic values and plot errors provided in
#' # the two example data frames 'df_gv_unstr' and 'df_error_bivar'.
#'
#' gv_df <- df_gv_unstr
#' error_df <- df_error_bivar
#'
#' pheno_df <- make_phenotypes(
#'   gv_df,
#'   error_df,
#'   randomise = TRUE
#' )
#' @export
make_phenotypes <- function(gv_df,
                            error_df,
                            randomise = FALSE) {
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

  if (randomise) {
    id <- unique(gv_df$ID)
    n_blocks_total <- nrow(gv_df) / length(id)

    gv_df$ORD <- unlist((lapply(seq_len(n_blocks_total), function(x) sample(id))))
    gv_df <- gv_df[order(gv_df$ENV, gv_df$REP, gv_df$ORD), ]
  }

  phe <- error_df[, !(colnames(error_df) %in% c("ENV", "BLOCK", "COL", "ROW"))] +
    gv_df[, !(colnames(gv_df) %in% c("ENV", "REP", "ID", "ORD"))]
  error_df <- error_df[, c("ENV", "BLOCK", "COL", "ROW")]
  ids <- factor(as.numeric(as.character(gv_df$ID)))

  pheno_df <- cbind(
    error_df,
    id = ids,
    phe
  )

  pheno_df <- pheno_df[order(pheno_df$ENV, pheno_df$COL, pheno_df$ROW), ]
  colnames(pheno_df) <- c("env", "block", "col", "row", "id", paste0("phe.Trait", 1:n_traits))
  return(pheno_df)
}
