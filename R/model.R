# this script define methods for the model class

# default model class methods

# translate the given model types to the right S3 class name.
# The model types could be a combination of different models
interpret_model_type <- function(model_type) {
  IMMIGRANT_TYPES <- c('CONST')
  DECAY_TYPES <- c('EXPN', 'EXP', 'mEXP', 'mEXPN', 'PL', 'mPL', 'PLN', 'mPLN')
  if (class(model_type) == 'hawkes_model_type') {
    return(model_type)
  } else if (length(model_type) == 1 && is.character(model_type)) {
    if (model_type %in% IMMIGRANT_TYPES) return(hawkes_model_type(hawkes_immigrant_type = model_type))
    else return(hawkes_model_type(hawkes_decay_type = model_type))
  } else if (length(model_type) == 2) {
    if (is.null(names(model_type))) {
      return(hawkes_model_type(hawkes_decay_type = DECAY_TYPES[DECAY_TYPES %in% model_type],
                               hawkes_immigrant_type = IMMIGRANT_TYPES[IMMIGRANT_TYPES %in% model_type]))
    } else {
      stopifnot(all(c('hawkes_decay_type', 'hawkes_immigrant_type') %in% names(model_type)))
      return(hawkes_model_type(hawkes_decay_type = model_type[['hawkes_decay_type']],
                               hawkes_immigrant_type = model_type[['hawkes_immigrant_type']]))
    }
  } else {
    stop('Unknown model type!')
  }
}

hawkes_model_type <- function(hawkes_decay_type = NULL, hawkes_immigrant_type = NULL) {
  ret <- list()
  if (!is.null(hawkes_decay_type)) ret$hawkes_decay_type <- hawkes_decay_type
  if (!is.null(hawkes_immigrant_type)) ret$hawkes_immigrant_type <- hawkes_immigrant_type
  class(ret) <- 'hawkes_model_type'
  ret
}

#' @export
as.character.hawkes_model_type <- function(x, ..) {
  paste(unlist(x), collapse = '_')
}

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
predict_final_popularity.default <- function(model) {
  stop('Unknown model type!')
}

#' @export
get_a1.default <- function(model) {
  stop('Unknown model type!')
}

#' @export
get_viral_score.default <- function(model, mu) {
  stop('Unknown model type!')
}

#' @export
get_model_intensity_at.default <- function(model, t) {
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

#' Predict the final popularity (event count) of give histories and
#' its model parameters.
#' @param model a model object provides data, model_type, observation_time
#' and model parameters
#' @return a vector of predicted final popularities whose length is the same
#' as the number of cascades in the provided model object
#' @export
predict_final_popularity <- function(model) {
  UseMethod('predict_final_popularity')
}

#' Calculating the expected size of first level of descendants
#' @param model a model object provides data, model_type, observation_time
#' and model parameters
#' @return a vector of the expected sizes of first level of descendants of the
#' given cascades
#' @export
get_a1 <- function(model) {
  UseMethod('get_a1')
}

#' Viral score is the total reaction of the system to a single promotion,
#' i.e. the expected cascade size started by a single event of magnitude
#' @param model a model object for computing the branching factor.
#' @param mu the magnitude of the initial event
#' @export
get_viral_score <- function(model, mu) {
  UseMethod('get_viral_score')
}

#' Compute the intensity value of a given model at time t
#' @param model a model object for computing the intensity value
#' @param t the given time to compute the intensity
#' @param cascade_index determine which cascade in the list of cascades to compute, defaults to 1
#' @export
get_model_intensity_at <- function(model, t, cascade_index = 1) {
  UseMethod('get_model_intensity_at')
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
