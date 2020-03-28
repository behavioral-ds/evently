# this script implements default methods for the hawkes class

#' Create a new hawkes model with parameters available
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
#' @param upper_bound model parameter upper bounds. A named vector where names are model parameters and
#' values are the largest possible values.
#' @param model_vars a named list of extra variables provided to hawkes objects
#' @export
new_hawkes <- function(model_type, par = NULL, data = NULL, init_par = NULL,
                       observation_time = NULL, lower_bound = NULL, upper_bound = NULL,
                       model_vars = NULL) {
  model <- list(model_type = model_type)
  class(model) <- c(paste0('hawkes_', model_type), 'hawkes')
  param_names <- get_param_names(model)

  if (!is.null(data)) model$data <- data

  # check if provided parameters are of the same length as required
  if (!is.null(par)) {
    if (length(param_names) != length(par)) stop('The provided parameters are not aligned with the required model parameters!')

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

  if (length(names(model_vars)) > 0) {
    for (var in names(model_vars)) {
      model[[var]] <- model_vars[[var]]
    }
  }

  model
}

#' @importFrom utils hasName
check_required_hawkes_fields <- function(model, fields) {
  stopifnot('hawkes' %in% class(model))
  for (f in fields) {
    stopifnot(hasName(model, f))
  }
}

#' @importFrom stats runif
generate_random_points.hawkes <- function(model) {
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

get_lower_bound.hawkes <- function(model) {
  lower_bound <- c(K = 1e-100, beta = 1e-100, c = 1e-100, theta = 1e-100, N = 1)
  lower_bound[names(lower_bound) %in% get_param_names(model)]
}

get_upper_bound.hawkes <- function(model) {
  upper_bound <- c(K = 10000, beta = 1.016 - 1e-100, c = 300, theta = 300, N = 1e7)
  upper_bound[names(upper_bound) %in% get_param_names(model)]
}

get_ampl_execution_options.hawkes <- function(model) {
  ''
}

#' @export
print.hawkes <- function(x, ...) {
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

#' @export
predict_final_popularity.hawkes <- function(model) {
  check_required_hawkes_fields(model, c('par', 'model_type', 'data', 'observation_time'))
  branching_factor <- get_branching_factor(model)
  if (branching_factor >= 1) {
    warning('Branching Factor greater than 1, not possible to predict the size(super critical)')
    return(rep(Inf, length(model$data)))
  } else {
    # calculating the expected size of first level of descendants
    a1s <- get_a1(model)
    # calculating the final value
    final_popularity <- vapply(model$data, nrow, FUN.VALUE = NA_integer_) + a1s / (1 - branching_factor)
    return(final_popularity)
  }
}

#' @export
get_viral_score.hawkes <- function(model, mu) {
  check_required_hawkes_fields(model, c('par', 'model_type'))
  if (mu == 0) return(0)
  branching_factor <- get_branching_factor(model)
  # branching factor cannot be larger than 1
  if (branching_factor >= 1) return(Inf)
  mu / (1 - branching_factor)
}

# ampl related functions --------------------------------------------------

# check the cascades before feeding it into ampl for fitting hawkes
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

  if (length(data) == 1 && (is.null(observation_time))) {
    observation_time <- max(data[[1]]$time)
  }

  if (all(is.null(observation_time)) ||
      any(observation_time <= 0) ||
      (length(observation_time) > 1 && length(observation_time) < length(data))) {
    stop('Please double check the observation time!')
  }
  # ampl doesn't recognize infinity so set to a large number
  # didn't set to machine max as it will cause error in AMPL
  observation_time[is.infinite(observation_time)] <- 1e20
  if (length(observation_time) == 1) {
    observation_time <- rep(observation_time, length(data))
  }

  data <- lapply(seq_along(data), function(i) {
    hist <- data[[i]]
    if (observation_time[i] < hist$time[nrow(hist)]) {
      warning('The provided observation time is smaller than the last observed event! Attempt to slice the data.')
      hist <- hist[hist$time < observation_time[i], ]
    }
    new_row <- data.frame(time = observation_time[i], magnitude = 0)
    hist <- rbind(hist[, c('time', 'magnitude')], new_row)
    offset_same_time(hist)
  })

  # sanity check
  if (any(sapply(data, function(d) is.unsorted(d$time)))) stop('Something went wrong! The time is not strictly increasing.')

  data
}

get_ampl_data_output.hawkes <- function(model) {
  # preprocess the data here
  data <- preprocess_data(model$data, model$observation_time)

  # prepare output datas here
  lengthes <- sapply(data, nrow)
  indexed_data <- do.call(rbind.data.frame, lapply(seq_along(data), function(i) {
    hist <- data[[i]]
    data.frame(index1 = i, index2 = seq(nrow(hist)), magnitude = hist$magnitude, time = hist$time)
  }))
  start_zeros <- sapply(data, function(hist) {
    which(hist$time > 0)[1] - 1
  })

  c(
    ampl_output_from_r('HL', length(data), 'atomic'),
    ampl_output_from_r('L', lengthes, 'vector'),
    ampl_output_from_r('ML', max(lengthes), 'atomic'),
    ampl_output_from_r(names = c('magnitude', 'time'), var = indexed_data, 'data.frame'),
    ampl_output_from_r('J0', start_zeros, 'vector')
  )
}

get_ampl_model_output.hawkes <- function(model) {
  ## if HawkesN, lower bound for population is at least as many as I have seen
  max_N <- max(sapply(model$data, nrow))
  if ("N" %in% names(model$init_par) && "N" %in% names(model$lower_bound) && model$lower_bound[["N"]] < max_N)
    model$lower_bound[["N"]] <- max_N

  ## correct initial parameters out of bounds -- should not happen, but ...
  model$init_par[is.finite(model$init_par) & (model$init_par > model$upper_bound)] <- model$upper_bound[is.finite(model$init_par) & (model$init_par > model$upper_bound)]
  model$init_par[is.finite(model$init_par) & (model$init_par < model$lower_bound)] <- model$lower_bound[is.finite(model$init_par) & (model$init_par < model$lower_bound)]

  ## first describe the data
  output.PARAM <- paste(
    '# define data',
    'param HL > 0;',
    'param ML > 0;',
    'param L {1..HL} >= 0;',
    'param magnitude {1..HL,1..ML} >= 0;',
    'param time {1..HL,1..ML} >= 0;',
    'param J0 {1..HL} >= 0;',
    '',
    sep = '\n'
  )

  ## next describe our initial parameters, if we have them
  ## assume that the variables we use are the ones in the init_par
  output.VAR <- "#define parameters"
  for (var in names(model$init_par)) {
    if (is.finite(model$init_par[var])) {
      output.VAR <- paste(output.VAR, paste('var ', var, ' := ', model$init_par[var],'; ', sep = ""), sep = '\n')
    }  else {
      output.VAR <- paste(output.VAR, paste('var ', var, '; ', sep = ""), sep = '\n')
    }
  }

  ## next, do the part of the objective function and kernel dependent constraints
  output.OJ <- paste(
    '',
    # allow for model-specific data output
    '# define objective function to maximize',
    'maximize Likelihood:',
    get_ampl_likelihood(model),
    sep = '\n'
  )
  output.contraints <- paste(
    '# define bounds and constraints',
    get_ampl_constraints(model),
    sep = '\n'
  )

  ## finally, construct boxes (bounds) on parameters
  output.BOX <- ""
  for (var in names(model$init_par)) {
    ## construct bound statement
    crt <- sprintf('subject to %s_limit:', var)
    if (is.finite(model$lower_bound[var])) crt <- paste(crt, model$lower_bound[var], "<=", sep = " ")
    crt <- paste(crt, var, sep = " ")
    if (is.finite(model$upper_bound[var])) crt <- paste(crt, "<=", model$upper_bound[var], sep = " ")
    crt <- paste(crt, ";", sep = "")
    if (is.finite(model$lower_bound[var]) && is.finite(model$upper_bound[var]))
      output.BOX <- paste(output.BOX, crt, sep = '\n')
  }

  ## finally peace it all together
  c(output.PARAM, output.VAR, output.OJ, output.contraints, output.BOX)
}
