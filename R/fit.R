# This script hosts functions that handle the data preprocessing and model selection.

#' Fit a Hawkes process or HawkesN process model on one or many event cascades
#' and learn model parameters.
#'
#' @param data a list of data.frame(s) where each data.frame is an event cascade with event
#' tims and event magnitudes (optional)
#' @param model_type a string representing the model type, e.g. EXP for Hawkes processes with
#' an exponential kernel function
#' @param cores the number of cores used for parallel fitting, defaults to 1 (non-parallel)
#' @param init_pars a data.frame of initial parameters passed to the fitting program. Parameters should be
#' aligned with required ones for the corresponding "model_type". The default initial parameters will
#' be used if not provided.
#' @param .init_no If initi_pars is not provided, currently 10 random starting parameters are generated
#' for fitting. This controls which random points are used. Defaults to NULL
#' @param observation_time the event cascades observation time(s). This can either be a single number indicating
#' a common observation time for all cascades or a vector of observation times which has the same length as
#' the number of cascades.
#' @param lower_bound model parameter lower bounds. A named vector where names are model parameters and
#' values are the lowest possible values.
#' @param upper_bound model parameter upper bounds. A named vector where names are model parameters and
#' values are the largest possible values.
#' @param model_vars a named list of extra variables provided to hawkes objects
#' @param parallel_type One of "PSOCK" or "FORK". Default to "PSOCK". See "Details" in makeCluster {parallel}.
#' @param ... further arguments passed to ampl
#' @import parallel
#' @export
fit_series <- function(data, model_type, cores = 1, init_pars, .init_no = NULL, observation_time = NULL,
                       lower_bound = NULL, upper_bound = NULL, model_vars = NULL, parallel_type = 'PSOCK', ...) {
  data <- preparation(data)
  model <- new_hawkes(data = data, model_type = model_type, observation_time = observation_time,
                            lower_bound = lower_bound, upper_bound = upper_bound, model_vars = model_vars)
  ampl_dat_file <- output_dat(model) # output ampl data file for being reused by fits with different initializations

  if (!missing(init_pars)) {
    ## use the provided init_pars
    if (all(get_param_names(model) %in% names(init_pars))) {
      points <- init_pars[names(init_pars) %in% get_param_names(model)]
    } else {
      stop('The provided initial parameters have a wrong parameter set!')
    }
  } else {
    ## get the initial points
    points <- generate_random_points(model)
  }
  models_with_initial_point <- lapply(seq(nrow(points)), function(i) {
    model$init_par <- unlist(points[i, , drop = F])
    model
  })

  ## if no .init_no, then do all and model selection at the end
  if (is.null(.init_no)) .init_no <- seq(models_with_initial_point)

  ## if we are asked for an init larger than our initial parameters, report errors
  if (sum(.init_no > nrow(models_with_initial_point)) > 0) stop('init_no is too large')

  inner_apply_func <- function(model, ...){
    # to make sure exact init_par is saved, mainly for lgo fitting
    saved_init_par <- model$init_par
    if (is.na(model$init_par[[1]])) {
      lgo_model <- ampl_run(model = model, solver = "lgo", dat_file = ampl_dat_file, ...)
      model[['init_par']] <- unlist(lgo_model$par)
    }
    model <- ampl_run(model = model, solver = "ipopt", dat_file = ampl_dat_file, ...)
    # to emphasize the init_par is found by lgo
    model$init_par <- saved_init_par

    return(model)
  }

  ## start fitting
  if (cores == 1 || length(models_with_initial_point) == 1) {
    fitted_models <- lapply(X = models_with_initial_point[.init_no], FUN = inner_apply_func, ...)
  } else {
    fitted_models <- switch (parallel_type,
      PSOCK = {
        cl <- makePSOCKcluster(cores)
        res <- parLapply(cl = cl, X = models_with_initial_point[.init_no], fun = inner_apply_func,
                         ampl_execution = .globals$execution, ...)
        stopCluster(cl)
        res
      },
      FORK = mclapply(X = models_with_initial_point[.init_no], FUN = inner_apply_func, ..., mc.cores = cores, mc.silent = F),
      stop('Unknown parallel type.')
    )
  }

  model_selection(models = fitted_models, cores = cores, dat_file = ampl_dat_file)
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
#' @export
get_hawkes_neg_likelihood_value <- function(model, ..., par, data, model_type, observation_time) {
  # par and data are required for computing log-likelihood values
  if (!missing(model)) {
    if (!missing(par)) model$par <- par
    if (!missing(data)) model$data <- data
    if (!missing(model_type)) model$model_type <- model_type
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
  # put data into a list if it is a data.frame
  if (is.data.frame(data)) {
    data <- list(data)
  }

  # validate data format
  if (!(is.list(data) && all(sapply(data, is.data.frame)))) {
    stop('The provided cascade(s) is in a wrong format.')
  }

  # check if ampl is available
  if (any(Sys.which(c('ampl', 'ipopt')) == '') && Sys.getenv('AMPL_PATH') == '') {
    stop('Please set up ampl and ipopt before fitting!')
  }

  data
}
