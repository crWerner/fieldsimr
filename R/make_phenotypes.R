#' Generate phenotypes - Combine genetic values and plot errors
#'
#' Creates a data frame of phenotypes by combining genetic values with plot errors
#' generated with the function \link[FieldSimR]{field_trial_error}.
#' Requires genetic values generated with the functions \link[FieldSimR]{compsym_asr_output}
#' or \link[FieldSimR]{unstr_asr_output},
#' or any data frame matching the description below.
#'
#' @param gv.df A data frame of genetic values. Must contain the columns 'env', genotype 'id', 'rep',
#'   and the genetic values for each trait.
#' @param error.df A data frame of plot errors. Must contain the columns 'env', 'block',
#'   'col', 'row', and the plot errors for each trait.
#' @param design.df A optional data frame of frequencies for generating incomplete block designs.
#'   Must contain the columns 'env', 'id', and 'nreps' indicating the number of replicates per individual for each environment.
#' @param randomise When \code{TRUE} (default), genotypes are randomly allocated to plots according to
#'   a randomized complete (or incomplete) block design.\cr
#'   \strong{Note:} Other experimental designs are being implemented and should be generated externally.
#' @param return.effects When \code{TRUE} (default is \code{FALSE}), a list is returned with additional
#'   entries containing the genetic values and plot errors for each trait.
#'
#' @return A data frame with columns 'env', 'block', 'column', 'row', genotype 'id', 'rep', and
#'   the phenotypes for each trait. When \code{return.effects = TRUE}, a list is returned with additional
#'   entries containing the genetic values and plot errors for each trait.
#'
#' @examples
#' # Generate phenotypes by combining the genetic values and plot errors provided
#' # in the two example data frames gv_df_unstr and error_df_bivar.
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
                            design.df = NULL,
                            randomise = TRUE,
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
    stop("'gv.df' must contain the columns 'env', 'id', 'rep', and the genetic values for each trait")
  }

  colnames(error.df)[grep("column", colnames(error.df))] <- "col"
  if (any(!c("env", "block", "col", "row") %in% colnames(error.df))) {
    stop("'error.df' must contain the columns 'env', 'block', 'col', 'row', and the plot errors for each trait")
  }

  ntraits <- ncol(gv.df) - 3
  if (ncol(error.df) - 4 != ntraits) {
    stop("Number of traits in 'gv.df' must match number of traits in 'error.df'")
  }

  design <- FALSE
  if (!is.null(design.df)) {
    design <- TRUE
    if (inherits(design.df, "table")) {
      design.df <- as.data.frame(t(design.df))
      colnames(design.df) <- c("env", "id", "nreps")
    }
    if (!is.data.frame(design.df)) {
      stop("'design.df' must be a data frame")
    }
    colnames(design.df)[grep("env|id|nreps", tolower(colnames(design.df)))] <- tolower(colnames(design.df))[grep("env|id|nreps", tolower(colnames(design.df)))]

    if (any(!c("env", "id", "nreps") %in% colnames(design.df))) {
      stop("'design.df' must contain the columns 'env', 'id' and 'nreps'")
    }
  }

  if (all(!grepl('\\D', error.df$env))) {
    error.df$env <- factor(as.numeric(as.character(error.df$env)))
  } else {error.df$env <- factor(as.character(error.df$env))}
  error.df$block <- factor(as.numeric(as.character(error.df$block)))
  error.df$col <- factor(as.numeric(as.character(error.df$col)))
  error.df$row <- factor(as.numeric(as.character(error.df$row)))
  error.df <- error.df[order(error.df$env, error.df$block), ]
  error.df <- unique(error.df)
  rownames(error.df) <- NULL

  if (all(!grepl('\\D', gv.df$env))) {
    gv.df$env <- factor(as.numeric(as.character(gv.df$env)))
  } else {gv.df$env <- factor(as.character(gv.df$env))}
  if (all(!grepl('\\D', gv.df$id))) {
    gv.df$id <- factor(as.numeric(as.character(gv.df$id)))
  } else {gv.df$id <- factor(as.character(gv.df$id))}
  gv.df$block <- gv.df$rep <- factor(as.numeric(as.character(gv.df$rep)))
  gv.df <- gv.df[order(gv.df$env, gv.df$block), ]
  gv.df <- unique(gv.df)
  rownames(gv.df) <- NULL

  if (design) {
    design.df <- droplevels(design.df[design.df$nreps != 0, ])
    if (all(!grepl('\\D', design.df$env))) {
      design.df$env <- factor(as.numeric(as.character(design.df$env)))
    } else {design.df$env <- factor(as.character(design.df$env))}
    if (all(!grepl('\\D', design.df$id))) {
      design.df$id <- factor(as.numeric(as.character(design.df$id)))
    } else {design.df$id <- factor(as.character(design.df$id))}
    design.df$nreps <- factor(as.numeric(as.character(design.df$nreps)))
    design.df <- design.df[order(design.df$env, design.df$id), ]
    design.df <- unique(design.df)
    rownames(design.df) <- NULL

    if (any(duplicated(design.df[, c("env", "id")]))) {
      stop("Individuals must not have multiple entries in 'design.df' for each environment")
    }
    design_reps <- with(design.df, tapply(nreps, env, unique))
    error_blocks <- with(error.df, tapply(block, env, function(x) c(0, 1, max(as.numeric(as.character(x))))))
    reps_conform <- mapply(function(x, y) any(!x %in% y), x = design_reps, y = error_blocks)
    if (any(reps_conform)) {
      stop("Number of replicates in 'design.df' must be 0, 1 or match number of blocks in 'error.df' for each environment")
    }
    design_plots <- with(design.df, tapply(nreps, env, function(x) sum(as.numeric(as.character(x)))))
    error_plots <- with(error.df, table(env))
    if (length(error_plots) != length(design_plots)) {
      stop("Number of environments in 'design.df' and 'error.df' must match")
    }
    if (any(design_plots - error_plots != 0)) {
      stop("Number of plots dictated by 'design.df' must match number of plots in 'error.df' for each environment")
    }

    gv_env_names <- levels(gv.df[["env"]])
    design_env_names <- levels(design.df[["env"]])
    if (length(gv_env_names) != length(design_env_names)) {
      stop("Number of environments in 'gv.df' and 'design.df' must match")
    }
    if (any(gv_env_names != design_env_names)) {
      stop("'env' names in 'gv.df' and 'design.df' must match")
    }
    gv_id_names <- with(gv.df, tapply(id, env, unique))
    design_id_names <- with(design.df, tapply(id, env, unique))
    missing_names <- mapply(function(x, y) any(!x %in% y), x = design_id_names, y = gv_id_names)
    if (any(missing_names)) {
      stop("All 'id' names in 'design.df' must be present in 'gv.df' for relevant environments")
    }

    gv.df$env.id <- paste0(gv.df$env, ":", gv.df$id)
    design.df$env.id <- paste0(design.df$env, ":", design.df$id)
    gv.df <- gv.df[gv.df$env.id %in% unique(design.df$env.id), ]
    gv.df <- gv.df[order(gv.df$env, gv.df$id), ]
    rownames(gv.df) <- NULL

    gv_reps <- as.data.frame(with(gv.df, table(id, env)))
    nreps <- as.numeric(as.character(design.df$nreps))
    if (any(nreps != gv_reps$Freq[gv_reps$Freq != 0])) {
      gv.df <- gv.df[!duplicated(gv.df$env.id), ]
      gv.df <- gv.df[rep(rownames(gv.df), times = nreps), ]
    }
    gv.df$block <- gv.df$rep <- factor(unlist(lapply(nreps, function(x) 1:x)))
    gv_blocks <- with(gv.df, tapply(block, env, function(x) 1:max(as.numeric(as.character(x))), simplify = FALSE))
    nleftovers <- table(design.df$env[design.df$nreps == 1]) / unlist(lapply(gv_blocks, function(x) max(x)))
    gv.df$block[gv.df$env.id %in% design.df$env.id[design.df$nreps == 1]] <- unlist(mapply(function(x, y) sample(rep(x, each = y)), x = gv_blocks, y = nleftovers))

    gv.df <- gv.df[order(gv.df$env, gv.df$block, gv.df$id), ]
    gv.df$env.id <- NULL
    gv.df <- unique(gv.df)
    rownames(gv.df) <- NULL
  }

  nenvs <- nlevels(gv.df$env)
  if (nlevels(error.df$env) != nenvs) {
    stop("Number of environments in 'gv.df' and 'error.df' must match")
  }

  if (nrow(gv.df) != nrow(error.df)) {
    gv.df$env.id <- paste0(gv.df$env, ":", gv.df$id)
    gv.df <- gv.df[!duplicated(gv.df$env.id), ]
    gv_id_names <- with(gv.df, tapply(id, env, function(x) length(unique(x))))
    block_plots <- with(error.df[error.df$block == 1, ], tapply(block, env, length))

    if (all(gv_id_names - block_plots == 0)) {
      error_blocks <- with(error.df, tapply(block, env, function(x) max(as.numeric(as.character(x)))))
      nreps <- rep(error_blocks, times = gv_id_names)
      gv.df <- gv.df[rep(rownames(gv.df), times = nreps), ]
      gv.df$block <- gv.df$rep <- factor(unlist(lapply(nreps, function(x) 1:x)))
      gv.df <- gv.df[order(gv.df$env, gv.df$block, gv.df$id), ]
      rownames(gv.df) <- NULL
      gv.df$env.id <- NULL
    } else {
      stop("Number of rows in 'gv.df' and 'error.df' do not match (after removing duplicated entries)")
    }
  }

  gv_reps <- with(gv.df, tapply(block, list(block, env), length))
  error_blocks <- with(error.df, tapply(block, list(block, env), length))
  if (any(dim(gv_reps) - dim(error_blocks) != 0)) {
    stop("Number of replicates in 'gv.df' must match number of blocks in 'error.df' for each environment")
  }
  if (any(gv_reps - error_blocks != 0, na.rm = TRUE)) {
    stop("Number of individuals per replicate in 'gv.df' must match number of plots per block in 'error.df' for each environment")
  }

  if (randomise) {
    gv.df$ord <- sample(1:nrow(gv.df))
    gv.df <- gv.df[order(gv.df$env, gv.df$block, gv.df$ord), ]
  }

  y <- error.df[, !(colnames(error.df) %in% c("env", "block", "col", "row"))] +
    gv.df[, !(colnames(gv.df) %in% c("env", "id", "rep", "block", "ord"))]

  gv_env_names <- as.character(gv.df[["env"]])
  error_env_names <- as.character(error.df[["env"]])
  if (any(gv_env_names != error_env_names)) warning("'env' names in 'gv.df' and 'error.df' do not match, names in 'gv.df' will be used")

  gv_df_names <- data.frame(env = factor(gv_env_names))
  error_df_names <- error.df[, c("block", "col", "row")]
  ids <- gv.df$id
  reps <- gv.df$rep

  pheno_df <- cbind(
    gv_df_names,
    error_df_names,
    id = ids,
    rep = reps,
    y
  )

  pheno_df <- pheno_df[order(pheno_df$env, pheno_df$col, pheno_df$row), ]
  rownames(pheno_df) <- NULL
  colnames(pheno_df) <- c("env", "block", "col", "row", "id", "rep", paste0("y.Trait", 1:ntraits))

  if (return.effects) {
    gv_ls <- as.list(as.data.frame(gv.df[, 4:(ntraits + 3)]))
    error_ls <- as.list(as.data.frame(error.df[, 5:(ntraits + 4)]))
    effects_df <- mapply(function(x, y) cbind(pheno_df[, 1:6], x, y), x = gv_ls, y = error_ls, SIMPLIFY = FALSE)

    effects_df <- mapply(function(x, y) {
      colnames(x)[7:8] <- c(paste0("gv.Trait", y), paste0("e.Trait", y))
      x
    }, x = effects_df, y = seq_len(ntraits), SIMPLIFY = FALSE)

    list_names <- c("pheno.df", paste0("Trait", 1:ntraits))
    pheno_df <- c(list(pheno_df), effects_df)
    names(pheno_df) <- list_names
  }

  return(pheno_df)
}
