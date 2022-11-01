#' Field trial error example data frame
#'
#' An example data frame of simulated field trial plot errors for two traits
#' tested in three environments. The data frame was generated using the function
#' \link[FieldSimR]{field_trial_error} with bivariate interpolation.
#'
#' @format A data frame with 700 rows and 5 columns:
#' \describe{
#'   \item{env}{Environment id}
#'   \item{block}{Block id}
#'   \item{col}{Column id}
#'   \item{row}{Row id}
#'   \item{e.Trait.1}{Simulated plot error for trait 1}
#'   \item{e.Trait.2}{Simulated plot error for trait 2}
#' }
"df_error_bivar"

#' Genetic values example data frame
#'
#' An example data frame of simulated genetic values for two traits tested in three environments.
#' The data frame was generated using the wrapper functions \link[FieldSimR]{unstr_asr_input} and
#' \link[FieldSimR]{unstr_asr_output} to simulate correlated genetic values based on an
#' unstructured model for genotype-by-environment (GxE) interaction with
#' \href{https://CRAN.R-project.org/package=AlphaSimR}{'AlphaSimR'}.
#'
#' @format A data frame with 700 rows and 5 columns:
#' \describe{
#'   \item{env}{Environment id}
#'   \item{rep}{Replicate number}
#'   \item{id}{Genotype id}
#'   \item{gv.Trait.1}{Simulated genetic values for trait 1}
#'   \item{gv.Trait.2}{Simulated genetic values for trait 2}
#' }
"df_gv_unstr"
