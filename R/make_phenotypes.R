#' Generate phenotypes - Combine genetic values and plot errors
#'
#' Creates a data frame of phenotypes by combining genetic values with plot errors
#' generated with the function \link[FieldSimR]{field_trial_error}.
#' Requires genetic values generated with the functions \link[FieldSimR]{compsym_asr_output}
#' or \link[FieldSimR]{unstr_asr_output},
#' or any data frame matching the description below.
#'
#' @param gv.df A data frame of genetic values. Must contain the columns 'env', 'rep', genotype 'id',
#'   and the genetic values for each trait.
#' @param error.df A data frame of plot errors. Must contain the columns 'env', 'block',
#'   'col', 'row', and the plot errors for each trait.
#' @param randomise When \code{TRUE} (default is \code{FALSE}), genotypes are randomly allocated to plots according to
#'   a randomized complete block design (RCBD).\cr
#'   \strong{Note:} Other experimental designs are being implemented and should be generated externally.
#' @param return.effects When \code{TRUE} (default is \code{FALSE}), a list is returned with additional
#'   entries containing the genetic values and plot errors for each trait.
#'
#' @return A data frame with columns 'env', 'block', 'column', 'row', genotype 'id', and
#'   the phenotypes for each trait. When \code{return.effects = TRUE}, a list is returned with additional
#'   entries containing the genetic values and plot errors for each trait.
#'
#' @examples
#' # Generate phenotypes by combining the genetic values and plot errors provided in
#' # the two example data frames gv_df_unstr and error_df_bivar.
#'
#' pheno_df <- make_phenotypes(
#'   gv.df = gv_df_unstr,
#'   error.df = error_df_bivar,
#'   randomise = TRUE,
#'   return.effects = TRUE
#' )
#'
#' @export
make_phenotypes <- function(gv.df,
                            error.df,
                            randomise = FALSE,
                            return.effects = FALSE) {
  if (inherits(gv.df, "list")) gv.df <- gv.df[[1]]
  if (inherits(error.df, "list")) error.df <- error.df[[1]]

  if (!is.data.frame(gv.df)) {
    stop("'gv.df' must be a data frame")
  }
  if (!is.data.frame(error.df)) {
    stop("'error.df' must be a data frame")
  }

  colnames(gv.df)[grep("env|rep|id", tolower(colnames(gv.df)))] <- tolower(colnames(gv.df))[grep("env|rep|id", tolower(colnames(gv.df)))]
  colnames(error.df)[grep("env|block|col|row", tolower(colnames(error.df)))] <- tolower(colnames(error.df))[grep("env|block|col|row", tolower(colnames(error.df)))]

  if (any(!c("env", "rep", "id") %in% colnames(gv.df))) {
    stop("'gv.df' must contain the columns 'env', 'rep', 'id', and the genetic values for each trait")
  }

  colnames(error.df)[grep("column", colnames(error.df))] <- "col"
  if (any(!c("env", "block", "col", "row") %in% colnames(error.df))) {
    stop("'error.df' must contain the columns 'env', 'block', 'col', 'row', and the plot errors for each trait")
  }

  ntraits <- ncol(gv.df) - 3
  if (ncol(error.df) - 4 != ntraits) {
    stop("'gv.df' must contain the columns 'env', 'rep', 'id', and the genetic values for each trait.
         'error.df' must contain the columns 'env', 'block', 'col', 'row', and the plot errors for each trait")
  }

  error.df$env <- factor(as.numeric(as.character(error.df$env)))
  error.df$block <- factor(as.numeric(as.character(error.df$block)))
  error.df$col <- factor(as.numeric(as.character(error.df$col)))
  error.df$row <- factor(as.numeric(as.character(error.df$row)))
  error.df <- error.df[order(error.df$env, error.df$block), ]
  error.df <- unique(error.df)

  gv.df$env <- factor(as.numeric(as.character(gv.df$env)))
  gv.df$rep <- factor(as.numeric(as.character(gv.df$rep)))
  gv.df$id <- factor(as.numeric(as.character(gv.df$id)))
  gv.df <- gv.df[order(gv.df$env, gv.df$rep), ]
  gv.df <- unique(gv.df)

  if (nrow(gv.df) != nrow(error.df)) stop("number of rows in 'gv.df' and 'error.df' must match")

  if (randomise) {
    id <- unique(gv.df$id)
    nblocks_total <- nrow(gv.df) / length(id)
    gv.df$ord <- unlist((lapply(seq_len(nblocks_total), function(x) sample(id))))
    gv.df <- gv.df[order(gv.df$env, gv.df$rep, gv.df$ord), ]
  }

  y <- error.df[, !(colnames(error.df) %in% c("env", "block", "col", "row"))] +
    gv.df[, !(colnames(gv.df) %in% c("env", "rep", "id", "ord"))]
  error.df.names <- error.df[, c("env", "block", "col", "row")]
  ids <- factor(as.numeric(as.character(gv.df$id)))

  pheno_df <- cbind(
    error.df.names,
    id = ids,
    y
  )

  pheno_df <- pheno_df[order(pheno_df$env, pheno_df$col, pheno_df$row), ]
  rownames(pheno_df) <- NULL
  colnames(pheno_df) <- c("env", "block", "col", "row", "id", paste0("y.Trait", 1:ntraits))

  if (return.effects) {
    gv_ls <- as.list(as.data.frame(gv.df[, 4:(ntraits + 3)]))
    error_ls <- as.list(as.data.frame(error.df[, 5:(ntraits + 4)]))
    effects_df <- mapply(function(x,y) cbind(pheno_df[,1:5], x, y), x = gv_ls, y = error_ls, SIMPLIFY = FALSE)

    effects_df <-  mapply(function(x,y) {
      colnames(x)[6:(ntraits + 5)] <- c(paste0("gv.Trait", y), paste0("e.Trait", y))
      x
    }, x = effects_df, y = seq_len(ntraits), SIMPLIFY = FALSE)

    list_names <- c("pheno.df", paste0("Trait", 1:ntraits))
    pheno_df <- c(list(pheno_df), effects_df)
    names(pheno_df) <- list_names
  }

  return(pheno_df)
}
