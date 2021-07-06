# This script hosts functions that handle the data preprocessing and model selection.

#' Fit a Hawkes process or HawkesN process model on one or many event cascades
#' and learn model parameters.
#'
#' @param data A list of data.frame(s) where each data.frame is an event cascade with event
#' tims and event magnitudes (optional)
#' @param model_type A string representing the model type, e.g. EXP for Hawkes processes with
#' an exponential kernel function
#' @param cores The number of cores used for parallel fitting, defaults to 1 (non-parallel)
#' @param init_pars A data.frame of initial parameters passed to the fitting program. Parameters should be
#' aligned with required ones for the corresponding "model_type". The default initial parameters will
#' be used if not provided.
#' @param .init_no If initi_pars is not provided, currently 10 random starting parameters are generated
#' for fitting. This controls which random points are used. Defaults to NULL
#' @param observation_time The event cascades observation time(s). This can either be a single number indicating
#' a common observation time for all cascades or a vector of observation times which has the same length as
#' the number of cascades.
#' @param lower_bound Model parameter lower bounds. A named vector where names are model parameters and
#' values are the lowest possible values.
#' @param upper_bound Model parameter upper bounds. A named vector where names are model parameters and
#' values are the largest possible values.
#' @param limit_event Define the way to optimize the computation by reducing the number of events added in log-likelihood (LL) functions,
#' defaults to NULL, i.e., no optimization. To limit the number of events computed, a list with `type` and `value` shoud be provided.
#' For example, limit_event = list(type = "event", value = 10) limits the LL fitting to 10 events,
#' limit_event = list(type = "time", value = 10)  limits the LL fitting to the events within past 10 time units.
#' The best practice to trade-off the computation could be to limit to the largest number of events that one can afford.
#' @param model_vars A named list of extra variables provided to hawkes objects
#' @param parallel_type One of "PSOCK" or "FORK". Default to "PSOCK". See "Details" in makeCluster {parallel}.
#' @param ... Further arguments passed to ampl
#' @import parallel
#' @return A model object where the [par] is fitted on [data]. [convergence] indicates the fitting convergence
#' status and [value] is the negative log-likelihood value of the fitted model on [data].
#' @export
#' @examples
#' \dontrun{
#' data <- generate_series(model_type = 'EXP',
#'                                      par = c(K = 0.9, theta = 1),
#'                                      sim_no = 10, Tmax = Inf)
#' fitted <- fit_series(data, 'EXP', observation_time = Inf)
#' fitted$par # fitted parameters
#' fitted$convergence # convergence status
#' fitted$value # negative log-likelihood value
#' }
fit_series <- function(data, model_type, cores = 1, init_pars = NULL, .init_no = NULL, observation_time = NULL,
                       lower_bound = NULL, upper_bound = NULL, limit_event = NULL, model_vars = NULL, parallel_type = 'PSOCK', ...) {
  data <- preparation(data)
  model <- new_hawkes(data = data, model_type = model_type, observation_time = observation_time,
                      lower_bound = lower_bound, upper_bound = upper_bound, model_vars = model_vars, limit_event = limit_event)
  fit_series_by_model(model, cores = cores, init_pars = init_pars, .init_no = .init_no,
                      parallel_type = parallel_type, ...)
}

#' Compute the negative log-likelihood values of a given model on a list of given
#' event cascades.
#'
#' @param model An object of a specific model class where the `data` and the `par` fields
#' are required.
#' @param ... Further arguments passed to ampl
#' @param par Hawkes model parameters
#' @param data A list of data.frames of event cascades
#' @param model_type The Hawkes model type
#' @param observation_time The observation time of the given event cascades
#' @return A single number, the negative log-likelihood of the given model on data
#' @export
#' @examples
#' \dontrun{
#' data <- generate_series(model_type = 'EXP',
#'                                      par = c(K = 0.9, theta = 1),
#'                                      sim_no = 10, Tmax = Inf)
#' fitted <- fit_series(data, 'EXP', observation_time = Inf)
#' data_test <- generate_series(model_type = 'EXP',
#'                                           par = c(K = 0.9, theta = 1),
#'                                           sim_no = 10, Tmax = Inf)
#' get_hawkes_neg_likelihood_value(fitted, data = data_test)
#' }
get_hawkes_neg_likelihood_value <- function(model, ..., par, data, model_type, observation_time) {
  # par and data are required for computing log-likelihood values
  if (!missing(model) && missing(model_type)) {
    if (!missing(par)) model$par <- par
    if (!missing(data)) model$data <- convert_data_format(data)
    check_required_hawkes_fields(model, c('par', 'data'))
  } else if (missing(par) || missing(data) || missing(model_type)) {
    stop('Neither an model object nor par,data,model_type are provided!')
  } else {
    model <- new_hawkes(model_type = model_type, par = par, data = data)
  }
  if (!missing(observation_time)) model$observation_time <- observation_time

  # a trick to reuse existing functions
  # have made sure this won't affect the original model object
  model$init_par <- model$par
  ampl_run(model, goal = 'nll', ...)
}

model_selection <- function(models, cores, ...) {
  if (length(models) == 1) return(models[[1]])

  ## score each model -- don't trust the algorithms own value, redo my own.
  nLLs <- simplify2array(mclapply(models, function(model) {
    get_hawkes_neg_likelihood_value(model, ...)
  }, mc.cores = cores), higher = F)
  if (all(is.na(nLLs))) stop('something went wrong! All neg.likelihood values are missing')

  models[[which.min(nLLs)]]
}

preparation <- function(data) {
  # check if ampl is available
  if (any(Sys.which(c('ampl', 'ipopt')) == '') && Sys.getenv('AMPL_PATH') == '') {
    stop('Please set up ampl and ipopt before fitting!')
  }

  convert_data_format(data)
}
