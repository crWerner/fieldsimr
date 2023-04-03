#' Plot errors - example data frame
#'
#' An example data frame of simulated plot errors for two traits across three locations. The first
#' two locations include three replicates (blocks), and the third location includes two replicates.
#' Each replicate comprises 10 rows and 10 columns. The data frame was generated using the
#' function \link[FieldSimR]{field_trial_error} with bivariate interpolation. The simulation of
#' the plot errors is shown in the vignette on the
#' \href{https://crwerner.github.io/fieldsimr/articles/spatial_error_demo.html}{Simulation of plot errors and phenotypes for a multi-environment plant breeding trial}.
#'
#' @format A data frame with 800 rows and 6 columns:
#' \describe{
#'   \item{env}{Environment number}
#'   \item{block}{Block (replicate) number}
#'   \item{col}{Column number}
#'   \item{row}{Row number}
#'   \item{e.Trait.1}{Simulated plot error for trait 1}
#'   \item{e.Trait.2}{Simulated plot error for trait 2}
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
#' \href{https://crwerner.github.io/fieldsimr/articles/unstructured_GxE_demo.html}{Simulation of genetic values using an unstructured GxE interaction model}.
#'
#' @format A data frame with 800 rows and 5 columns:
#' \describe{
#'   \item{env}{Environment number}
#'   \item{rep}{Replicate number}
#'   \item{id}{Genotype id}
#'   \item{gv.Trait.1}{Simulated genetic values for trait 1}
#'   \item{gv.Trait.2}{Simulated genetic values for trait 2}
#' }
"df_gv_unstr"
