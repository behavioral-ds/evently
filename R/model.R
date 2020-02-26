# this script define methods for the model class

# default model class methods

get_ampl_likelihood.default <- function(model) {
  stop('Unknown model type!')
}

get_ampl_constraints.default <- function(model) {
  stop('Unknown model type!')
}

#' @export
get_branching_factor.default <- function(model) {
  stop('Unknown model type!')
}

#' @export
get_viral_score.default <- function(model) {
  stop('Unknown model type!')
}

# function dispatchers
generate_random_points <- function(obj) {
  UseMethod('generate_random_points', obj)
}

get_param_names <- function(obj) {
  UseMethod('get_param_names')
}

get_lower_bound <- function(obj) {
  UseMethod('get_lower_bound')
}

get_upper_bound <- function(obj) {
  UseMethod('get_upper_bound')
}

get_ampl_likelihood <- function(obj) {
  UseMethod('get_ampl_likelihood')
}

get_ampl_constraints <- function(obj) {
  UseMethod('get_ampl_constraints')
}

#' Branching factor is the expected number of events generated
#' by a single event.
#' @param model a model object for computing the branching factor.
#' @export
get_branching_factor <- function(model) {
  UseMethod('get_branching_factor')
}

#' Viral score is the total reaction of the system to a single promotion,
#' i.e. the expected cascade size started by a single event of magnitude
#' @param model a model object for computing the branching factor.
#' @param mu the magnitude of the initial event
#' @export
get_viral_score <- function(model, mu) {
  UseMethod('get_viral_score')
}

get_ampl_data_output <- function(obj) {
  UseMethod('get_ampl_data_output')
}

get_ampl_model_output <- function(obj) {
  UseMethod('get_ampl_model_output')
}

get_ampl_execution_options <- function(obj) {
  UseMethod('get_ampl_execution_options')
}
