#' Plot errors - example data frame
#'
#' An example data frame of simulated plot errors for two traits across three locations. The first
#' and the third location include two blocks, and the second location includes three blocks.
#' Each block comprises 20 rows and 5 columns. The data frame was generated using the
#' function \link[FieldSimR]{field_trial_error} with bivariate interpolation. The simulation of
#' the plot errors is shown in the vignette on the
#' \href{https://crwerner.github.io/fieldsimr/articles/spatial_error_demo.html}{Simulation of plot errors and phenotypes for a multi-environment plant breeding trial}.
#'
#' @format A data frame with 700 rows and 6 columns:
#' \describe{
#'   \item{env}{Environment id}
#'   \item{block}{Block id}
#'   \item{col}{Column id}
#'   \item{row}{Row id}
#'   \item{e.Trait1}{Simulated plot error for trait 1}
#'   \item{e.Trait2}{Simulated plot error for trait 2}
#' }
"df_error_bivar"

#' Genetic values - example data frame
#'
#' An example data frame of simulated genetic values for 100 genotypes with two traits across
#' three environments. The data frame was generated using the wrapper functions
#' \link[FieldSimR]{unstr_asr_input} and \link[FieldSimR]{unstr_asr_output} to simulate correlated
#' genetic values based on an unstructured model for genotype-by-environment (GxE) interaction
#' with \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR'}. The simulation of
#' the genetic values is shown in the vignette on the
#' \href{https://crwerner.github.io/fieldsimr/articles/unstructured_GxE_demo.html}{Simulation of genetic values using an unstructured model for GxE interaction}.
#'
#' @format A data frame with 700 rows and 5 columns:
#' \describe{
#'   \item{env}{Environment id}
#'   \item{rep}{Replicate id}
#'   \item{id}{Genotype id}
#'   \item{gv.Trait1}{Simulated genetic values for trait 1}
#'   \item{gv.Trait2}{Simulated genetic values for trait 2}
#' }
"gv_df_unstr"
