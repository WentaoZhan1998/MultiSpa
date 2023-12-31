% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{interaction_gen_group}
\alias{interaction_gen_group}
\title{Sequence optimize}
\usage{
interaction_gen_group(
  neighbor_mat,
  seuInt,
  expr_mat,
  sampletable,
  process_fun,
  id = NULL,
  id_temp_fun = NULL,
  markers = NULL,
  markers2 = NULL,
  group1 = NULL,
  group2 = NULL,
  outnames = c("variable", "neighbor_variables"),
  expr_needed = F,
  zero_rm = F
)
}
\arguments{
\item{design}{matrix, the design matrix X for the predictors.}

\item{outcome}{vector, the output Y for the response.}

\item{model}{the model to be used. The default is set as 'linear' to use the linear model.}

\item{option}{a list specifies the options for model solving, currently supporting 'mle_solver' only.}
}
\value{
A list including the information of the fitted model
}
\description{
Do sequence optimize
}
