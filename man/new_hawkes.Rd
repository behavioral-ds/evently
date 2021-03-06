% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hawkes.R
\name{new_hawkes}
\alias{new_hawkes}
\title{Create a new hawkes model with given arguments}
\usage{
new_hawkes(
  model_type,
  par = NULL,
  data = NULL,
  init_par = NULL,
  observation_time = NULL,
  lower_bound = NULL,
  upper_bound = NULL,
  model_vars = NULL,
  limit_event = NULL
)
}
\arguments{
\item{model_type}{A string indicates the model tyep, e.g. EXP for a Hawkes process
with an exponential kernel}

\item{par}{A named vector denotes the model parameters where the names are model
parameters and the values are the corresponding parameter values}

\item{data}{A list of data.frame(s) where each data.frame is an event cascade with event
tims and event magnitudes (optional)}

\item{init_par}{Initial parameter values used in fitting}

\item{observation_time}{The event cascades observation time. It is assumed that all cascades in data
are observed until a common time.}

\item{lower_bound}{Model parameter lower bounds. A named vector where names are model parameters and
values are the lowest possible values.}

\item{upper_bound}{Model parameter upper bounds. A named vector where names are model parameters and
values are the largest possible values.}

\item{model_vars}{A named list of extra variables provided to hawkes objects}

\item{limit_event}{choose how to optimize the computation by reducing the number of events added in log-likelihood functions.}
}
\value{
A model object with class [hawkes] and [hawkes_`model_type`] where `model_type` is replaced
by the given model_type
}
\description{
Create a new hawkes model with given arguments
}
\examples{
data <- list(data.frame(time = c(0, 0.5, 1)))
new_hawkes(model_type = 'EXP', par = c(K = 0.9, theta = 1),
           data = data, observation_time = Inf)
}
