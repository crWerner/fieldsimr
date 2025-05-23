% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{gv_df_unstr}
\alias{gv_df_unstr}
\title{Genetic values - Example data frame}
\format{
A data frame with 700 rows and 5 columns:
\describe{
  \item{env}{Environment number}
  \item{id}{Genotype identifier}
  \item{rep}{Replicate number}
  \item{gv.Trait1}{Simulated genetic values for Trait 1}
  \item{gv.Trait2}{Simulated genetic values for Trait 2}
}
}
\usage{
gv_df_unstr
}
\description{
An example data frame with simulated genetic values of 100 genotypes for two traits in three environments.
Environments 1 and 2 comprise two replicates of each genotype, while Environment 3 comprises three replicates.
The data frame was generated using the wrapper functions
\link[FieldSimR]{unstr_asr_input} and \link[FieldSimR]{unstr_asr_output}, which simulate correlated
genetic values with \href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`}. The simulation is
demonstrated in the vignette
\href{https://crwerner.github.io/fieldsimr/articles/unstructured_GxE_demo.html}{Simulation of genetic values based on an unstructured model for GxE interaction}.
}
\keyword{datasets}
