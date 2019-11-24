# this script implements default methods for the hawkes_model class

# default model class methods

# Create a new hawkes model with parameters available
#' @param model_type A string indicates the model tyep, e.g. EXP for a Hawkes process
#' with an exponential kernel
#' @param par A named vector denotes the model parameters where the names are model
#' parameters and the values are the corresponding parameter values
#' @param data a list of data.frame(s) where each data.frame is an event cascade with event
#' tims and event magnitudes (optional)
#' @param observation_time the event cascades observation time. It is assumed that all cascades in data
#' are observed until a common time.
#' @param init_par Initial parameter values used in fitting
#' @param lower_bound model parameter lower bounds. A named vector where names are model parameters and
#' values are the lowest possible values.
#' @param uppper_bound model parameter upper bounds. A named vector where names are model parameters and
#' values are the largest possible values.
#' @export
new_hawkes_model <- function(model_type, par = NULL, data = NULL, init_par = NULL,
                             observation_time = NULL, lower_bound = NULL, upper_bound = NULL) {
  model <- list(model_type = model_type)
  class(model) <- c(paste0('hawkes_', model_type), 'hawkes_model')
  param_names <- get_param_names(model)

  if (!is.null(data)) model$data <- data

  # check if provided parameters are of the same length as required
  if (!is.null(par)) {
    stopifnot(length(param_names) == length(par))

    if (length(names(par)) == 0) {
      warning(paste0('Provided parameter vector is unnamed. Aussming the following order: ', paste(param_names, collapse = ', ')))
      names(par) <- param_names
    }

    par <- par[param_names]
  } else {
    par <- rep(NA, length(param_names))
    names(par) <- param_names
  }

  if (is.null(init_par)) {
    init_par <- rep(NA, length(param_names))
    names(init_par) <- param_names
  }

  model$init_par <- init_par
  model$par <- par
  model$value <- NA

  final_lower_bound <- get_lower_bound(model)
  final_upper_bound <- get_upper_bound(model)
  if (!is.null(lower_bound)) {
    if (length(lower_bound) <= length(final_lower_bound)) {
      final_lower_bound[names(lower_bound)] <- lower_bound
    } else if (length(lower_bound) > length(final_lower_bound)) {
      stop('Wrong lower bound provided!')
    }
  }

  if (!is.null(upper_bound)) {
    if (length(upper_bound) <= length(final_upper_bound)) {
      final_upper_bound[names(upper_bound)] <- upper_bound
    } else if (length(upper_bound) > length(final_upper_bound)) {
      stop('Wrong upper bound provided!')
    }
  }
  model$lower_bound <- final_lower_bound
  model$upper_bound <- final_upper_bound

  model$observation_time <- observation_time

  model
}

check_required_hawkes_model_fields <- function(model, fields) {
  for (f in fields) {
    stopifnot(hasName(model, f))
  }
}

#' @importFrom stats runif
generate_random_points.hawkes_model <- function(model) {
  init_k <- runif(5, min = .Machine$double.eps, max = 10)
  init_beta <- runif(5, min = .Machine$double.eps, max = 1)
  init_c <- runif(5, min = .Machine$double.eps, max = 300)
  init_theta <- runif(5, min = .Machine$double.eps, max = 3)
  init_n <- floor(runif(5, min = 50, max = 5000))
  ## 3 known good start points, all NA (which invokes the global optimizer) and all Inf (which lets IPOPT chose its own start point)
  init_k[6:10] <- c(0.1, 0.5, 0.8, NA, Inf)
  init_beta[6:10] <- c(0.001, 0.5, 0.8, NA, Inf)
  init_c[6:10] <- c(10, 100, 60, NA, Inf)
  init_theta[6:10] <- c(0.0001, 0.00001, 0.000001, NA, Inf)
  init_n[6:10] <- c(60, 200,600, NA, Inf)
  init <- data.frame('K' = init_k, 'beta' = init_beta,
                     'c' = init_c, 'theta' = init_theta, 'N' = init_n)

  init[names(init) %in% get_param_names(model)]
}

get_lower_bound.hawkes_model <- function(model) {
  lower_bound <- c(K = 1e-100, beta = 1e-100, c = 1e-100, theta = 1e-100, N = 1)
  lower_bound[names(lower_bound) %in% get_param_names(model)]
}

get_upper_bound.hawkes_model <- function(model) {
  upper_bound <- c(K = 10000, beta = 1.016 - 1e-100, c = 300, theta = 300, N = 1e7)
  upper_bound[names(upper_bound) %in% get_param_names(model)]
}

get_ampl_model_output.hawkes_model <- function(model) {
  ''
}

get_ampl_data_output.hawkes_model <- function(model) {
  ''
}

#' @export
print.hawkes_model <- function(x, ...) {
  for (n in names(x)) {
    if (n %in% 'data') {
      cat(paste('No. of cascades:', length(x$data), '\n'))
    } else if (n %in% c('model_type')) {
      cat(paste('Model:', x$model_type, '\n'))
    } else if (n %in% c('value')) {
      cat(paste('Neg Log Likelihood:', x$value, '\n'))
    } else if (n %in% c('convergence')) {
      cat(paste('convergence:', x$convergence, '\n'))
    } else if (n %in% c('init_par', 'par', 'lower_bound', 'upper_bound')) {
      cat(paste0(n, '\n'))
      cat('  ')
      cat(paste(names(x[[n]]), formatC(x[[n]], format = "e", digits = 2), collapse = '; '))
      cat('\n')
    }
  }
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

#' @export
get_branching_factor <- function(obj) {
  UseMethod('get_branching_factor')
}

get_ampl_data_output <- function(obj) {
  UseMethod('get_ampl_data_output')
}

get_ampl_model_output <- function(obj) {
  UseMethod('get_ampl_model_output')
}
