# implementation of a common benchmark SEISMIC as a model

fit_series_by_model.hawkes_SEISMIC <- function(model, cores, init_pars, parallel_type, .init_no, ...) {
  return(model)
}

get_param_names.hawkes_SEISMIC <- function(model) {
  ''
}

#' @export
predict_final_popularity.hawkes_SEISMIC <- function(model, data = NULL, observation_time = NULL) {
  if (!is.null(data) && !is.null(observation_time)) {
    model$data <- if (is.data.frame(data)) list(data) else data
    model$observation_time <- observation_time
  }
  check_required_packages('seismic')
  library(seismic)
  sum(sapply(seq_along(model$data), function(i) {
    cascade <- model$data[[i]]
    observation_time <- if (length(model$observation_time) == length(model$data)) model$observation_time[[i]] else model$observation_time
    if (!('magnitude' %in% names(cascade))) stop('User magnitudes are required for SEISMIC!')
    infectiousness <- get.infectiousness(cascade$time, cascade$magnitude, observation_time)
    # add 1 here as seismic prediction doesn't take the initial tweet into account
    pred.cascade(observation_time, infectiousness$infectiousness, cascade$time, cascade$magnitude, n.star = 100)[1, 1] + 1
  }))
}
