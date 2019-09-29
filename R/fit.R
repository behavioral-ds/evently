# This script hosts functions that handle the data preprocessing and model selection.

# for backward compactability
#' @export
fitSeries <- fitSeries <- function(history, parallel = T, kernel.type = 'EXPN', .init_no = NA, .additional.PLN.final.fit = F, .fit_engine = "AMPL", single_event_cascade = F, lowerBound = NULL, upperBound = NULL, ...) {
  cores <- 1
  if (parallel) cores <- detectCores()
  fit_series(history, model_type = kernel.type, cores = cores, lower_bound = lowerBound, upper_bound = upperBound, ...)
}


#' @import parallel
#' @export
fit_series <- function(history, model_type, cores = 1, .init_no = NULL, observation_time = NULL,
                       lower_bound = NULL, upper_bound = NULL, ...) {
  model <- new_hawkes_model(data = history, model_type = model_type, observation_time = observation_time,
                            lower_bound = lower_bound, upper_bound = upper_bound)

  ## get the initial points
  points <- generate_random_points(model)
  models_with_initial_point <- lapply(seq(nrow(points)), function(i) {
    model$init_par <- unlist(points[i, ])
    model
  })

  ## if no .init_no, then do all and model selection at the end
  if (is.null(.init_no)) .init_no <- seq(models_with_initial_point)

  ## if we are asked for an init larger than our initial params, report errors
  if (sum(.init_no > nrow(models_with_initial_point)) > 0) stop('init_no is too large')

  inner_apply_func <- function(i){
    model <- models_with_initial_point[[i]]

    tryCatch({
      lgo_model <- NA
      if (is.na(model$init_par[[1]])) {
        lgo_model <- ampl_run(model = model, solver = "lgo", ...)
        model[['init_par']] <- unlist(lgo_model$par)
      }
      model <- ampl_run(model = model, solver = "ipopt", ...)
    }, error = function(err) {
      print(paste("[fitSeries] Error in optim:  ", err))
    })

    return(model)
  }

  ## start fitting
  fitted_models <- mclapply(X = (1:nrow(initial_points))[.init_no], FUN = inner_apply_func, mc.cores = cores, mc.silent = F)

  model_selection(models = fitted_models, ...)
}

model_selection <- function(models, ...) {
  if (length(models) == 1) return(models[[1]])

  ## score each model -- don't trust the algorithms own value, redo my own.
  nLLs <- sapply(models, function(model) neg.log.likelihood(model, ...))
  if (all(is.na(nLLs))) stop('something went wrong! All neg.likelihood values are missing')

  models[[which.min(nLLs)]]
}
