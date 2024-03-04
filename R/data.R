#' Plot errors - Example data frame
#'
#' An example data frame with simulated plot errors for two traits in three environments.
#' Environments 1 and 2 comprise two blocks, while Environment 3 comprises three blocks.
#' The blocks are aligned in the column direction (side-by-side) and comprise 5 columns and 20 rows.
#' The data frame was generated using the function \link[FieldSimR]{field_trial_error} with
#' bivariate interpolation. The simulation is demonstrated in the vignette
#' \href{https://crwerner.github.io/fieldsimr/articles/spatial_error_demo.html}{Simulation of plot errors and phenotypes in multi-environment field trials}.
#'
#' @format A data frame with 700 rows and 6 columns:
#' \describe{
#'   \item{env}{Environment number}
#'   \item{block}{Block number}
#'   \item{col}{Column number}
#'   \item{row}{Row number}
#'   \item{e.Trait1}{Simulated plot errors for Trait 1}
#'   \item{e.Trait2}{Simulated plot errors for Trait 2}
#' }
"error_df_bivar"

#' Genetic values - Example data frame
#'
#' An example data frame with simulated genetic values of 100 genotypes for two traits in three environments.
#' Environments 1 and 2 comprise two replicates of each genotype, while Environment 3 comprises three replicates.
#' The data frame was generated using the wrapper functions
#' \link[FieldSimR]{unstr_asr_input} and \link[FieldSimR]{unstr_asr_output}, which simulate correlated
#' genetic values in \href{https://CRAN.R-project.org/package=AlphaSimR}{AlphaSimR}. The simulation is
#' demonstrated in the vignette
#' \href{https://crwerner.github.io/fieldsimr/articles/unstructured_GxE_demo.html}{Simulation of genetic values using an unstructured model for GxE interaction}.
#'
#' @format A data frame with 700 rows and 5 columns:
#' \describe{
#'   \item{env}{Environment number}
#'   \item{rep}{Replicate number}
#'   \item{id}{Genotype identifier}
#'   \item{gv.Trait1}{Simulated genetic values for Trait 1}
#'   \item{gv.Trait2}{Simulated genetic values for Trait 2}
#' }
"gv_df_unstr"
