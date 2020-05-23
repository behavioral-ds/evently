# this script implements methods for the background rate type: constant, i.e. a homogeneous poisson

get_param_names.hawkes_CONST <- function(model) {
  'lambda'
}

get_ampl_likelihood.hawkes_CONST <- function(model) {
  stop('Not implemented yet!')
}

get_ampl_constraints.hawkes_CONST <- function(model) {
  ''
}

generate_random_points.hawkes_CONST <- function(model) {
  data.frame(lambda = c(stats::runif(8, min = .Machine$double.eps, max = 300),
                        NA, Inf))
}

get_lower_bound.hawkes_CONST <- function(model) {
  c(lambda = 1e-100)
}

get_upper_bound.hawkes_CONST <- function(model) {
  c(lambda = 1e6)
}
