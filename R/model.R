# this script implements default methods for the hawkes_model class

#' @importFrom stats runif
preprocess_data <- function(data, observation_time) {
  offset_same_time <- function(history) {
    # check if two retweets happened at the same time. add a random number if so
    i <- 1
    while (i <= nrow(history) && history$time[i] == 0) i <- i + 1
    if ((nrow(history) + 1) == i) return(history)
    k <- history$time[i]
    i <- i + 1
    while (i <= nrow(history)) {
      if (history$time[i] == k) {
        history$time[i] <- history$time[i-1] + runif(1, 0, 1e-5)
      } else {
        k <- history$time[i]
      }
      i <- i + 1
    }
    history
  }

  if (length(data) == 1 && (is.null(observation_time) || max(data[[1]]$time) > observation_time)) {
    observation_time <- max(data[[1]]$time)
  }

  if (is.null(observation_time)) stop('Please specify an observation time when doing joint fitting!')
  data <- lapply(data, function(hist) {
    hist <- offset_same_time(hist)
    new_row <- data.frame(time = observation_time, magnitude = 0)
    hist <- rbind(hist, new_row)
    hist
  })

  # sanity check
  if (any(sapply(data, function(d) is.unsorted(d$time)))) stop('something went wrong..')

  data
}

# default model class methods

# create a new hawkes model class
new_hawkes_model <- function(data, model_type, init_par = NULL, observation_time = NULL,
                             lower_bound = NULL, upper_bound = NULL) {
  data <- preprocess_data(data, observation_time)
  model <- list(
    data = data,
    model_type = model_type
  )
  class(model) <- c(paste0('hawkes_', model_type), 'hawkes_model')
  param_names <- get_param_names(model)
  if (is.null(init_par)) {
    init_par <- rep(NA, length(param_names))
    names(init_par) <- param_names
  }

  model$init_par <- init_par
  model$par <- rep(NA, length(param_names))
  names(model$par) <- param_names
  model$value <- NA

  if (is.null(lower_bound)) lower_bound <- get_lower_bound(model)
  if (is.null(upper_bound)) upper_bound <- get_upper_bound(model)
  model$lower_bound <- lower_bound
  model$upper_bound <- upper_bound

  model$observation_time <- observation_time
  model
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

  init[, names(init) %in% get_param_names(model)]
}

get_lower_bound.hawkes_model <- function(model) {
  lower_bound <- c(K = 1e-100, beta = 1e-100, c = 1e-100, theta = 1e-100, N = 1)
  lower_bound[names(lower_bound) %in% get_param_names(model)]
}

get_upper_bound.hawkes_model <- function(model) {
  upper_bound <- c(K = 10000, beta = 1.016 - 1e-100, c = 300, theta = 300, N = 1e7)
  upper_bound[names(upper_bound) %in% get_param_names(model)]
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
