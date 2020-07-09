#' Function for sampling from the powerlaw distribution of user influence
#' it is the equivalent of the Richter-Gutenberg distribution in the Helmstetter model
#' the powerlaw distribution was determined from the twitter data, from the #retweets
#' alpha = 2.016, xmin = 1. Draw n values
#' @param n The number of samples to be generated
#' @param alpha Powerlaw distribution parameters
#' @param mmin Powerlaw distribution parameters
#' @import poweRlaw
#' @return A single number, a random user magnitude
#' @export
#' @examples
#' generate_user_magnitude(n = 1)
generate_user_magnitude <- function(n, alpha = 2.016, mmin = 1) {
  if ( is.null(.globals$user_infl) ) {
    .globals$user_infl <- conpl$new()
  }
  .globals$user_infl$setXmin(mmin)
  .globals$user_infl$setPars(alpha)

  return(dist_rand(.globals$user_infl, n))
}

##################################### SIMULATION OF THE HAWKES PROCESS ##############################################

# Generate Hawkes porcess events without the background rate
generate_series_no_background_rate <- function(par, model_type, Tmax = Inf, maxEvents = NULL, M = NULL, history_init = NULL, tol = 1e-5) {
  # determine if this is a marked model or not
  marked <- T
  if (substr(model_type, 1, 1) != 'm') {
    marked <- F
    M <- 1
  }

  if (is.null(M)) {
    M <- generate_user_magnitude(1)
  }

  # initial event, M and time 0
  history <- as.data.frame(matrix(c(M, 0), nrow=1, byrow = TRUE))
  colnames(history) <- c("magnitude", "time")

  if (!is.null(history_init)) {
    history <- history_init[, c("magnitude", "time")]
  }
  t <- utils::tail(history$time, n = 1)
  model <- new_hawkes(model_type = model_type, par = par)
  while(t <= Tmax){
    # generate the time of the next event, based on the previous events
    t.lambda.max <- t
    model$data <- list(history)
    intensityMax <- get_model_intensity_at(model, t = t.lambda.max)
    ## if intensityMax is too small then cut
    if (intensityMax < tol) break

    ## first, sample one event time from the maximum intensity
    r <- runif(1)
    t <- t - log(r) / intensityMax
    if (t > Tmax) break

    ## then perform rejection sampling
    s <- runif(1)
    model$data <- list(history)
    thr <- get_model_intensity_at(model, t = t) / intensityMax

    if (s <= thr) {
      ## if here, it means we accept the event
      # generate the influence of the next event, by sampling the powerlaw distribution of the #retweets
      mag <- ifelse(marked, generate_user_magnitude(n = 1), 1)

      # add the next event to the history
      event <- matrix( c(mag, t), nrow=1, byrow = T)
      colnames(event) <- c("magnitude", "time")

      ## now, find out where to insert the new event
      pos <- findInterval(x = event[1, "time"], vec = unlist(history$time) )
      ## if there is a tail in history (i.e. there are later events than the one sampled just now)
      if (pos < nrow(history) ) {
        history <- rbind(history[1:pos,], event, history[(pos+1):nrow(history),])
      } else {
        history <- rbind(history[1:pos,], event)
      }
    }

    if (!is.null(maxEvents) && nrow(history) >= maxEvents)
      break
  }

  return(history)
}

get_new_magnitude <- function(model_type = NULL) {
  if (!is.null(model_type) && substr(model_type, 1, 1) == 'm') {
    return(generate_user_magnitude(1))
  } else {
    return(1)
  }
}

generate_immigrant_event_series.CONST <- function(par, model_type, Tmax) {
  stopifnot('lambda' %in% names(par) || !is.infinite(Tmax))
  # first event is at time 0
  immigrants <- data.frame(magnitude = get_new_magnitude(model_type$hawkes_decay_type), time = 0)

  repeat {
    new_event <- data.frame(magnitude = get_new_magnitude(model_type$hawkes_decay_type),
                            time = utils::tail(immigrants$time, n = 1) + stats::rexp(1, rate = par[['lambda']]))
    if (new_event$time > Tmax) break()
    immigrants <- rbind(immigrants, new_event)
  }
  immigrants
}

generate_immigrant_event_series <- function(par, model_type, Tmax) {
  switch (model_type$hawkes_immigrant_type,
    CONST = generate_immigrant_event_series.CONST(par, model_type, Tmax)
  )
}

#' Main function to generate a Hawkes process sequence. It allows intermediary
#' saves and continuing a stopped simulation. Creates a CSV file with two
#' columns, each row is an event: (magnitude, time)
#' @param model A model class object with par and model_type presented. par and
#' model_type are not requird once this is given
#' @param par A named vector of model parameters, K, alpha, beta, mmin, c, theta - parameters of the Hawkes kernel
#' @param model_type Model type
#' @param sim_no The number of simulated cascades
#' @param cores The number of cores (processes) used for simulation
#' @param Tmax Maximum time of simulation.
#' @param maxEvents Maximum number of events to be simulated.
#' @param M Magnitude of the initial event
#' @param tol Simulation stops when intensity smaller than tol.
#' @param return_as_object wether return the cascades within a model object
#' @param init_history If given, the simulation will start after the last
#' event in the given init_history
#' @return A list of data.frames where each data.frame is a simulated event
#' cascade with the given model by default. Or a model object with the data.frames
#' if return_as_object is True
#' @export
#' @examples
#' generate_series(model_type = 'EXP',
#'                              par = c(K = 0.9, theta = 1),
#'                              sim_no = 10, Tmax = Inf)
generate_series <- function(model, par, model_type, sim_no = 1, cores = 1, Tmax = Inf,
                            maxEvents = NULL, M = NULL, tol = 1e-5, return_as_object = F,
                            init_history = NULL) {
  # stopifnot(is.null(history_init) || is.data.frame(history_init))
  if (!missing(model) && (!missing(par) || !missing(model_type))) {
    stop('Please either provide a model or (par, model_type) instead of both.')
  } else if (!missing(model)) {
    check_required_hawkes_fields(model, c('par', 'model_type'))
    par <- model$par
    model_type <- model$model_type
  }
  model_type <- interpret_model_type(model_type)

  data <- mclapply(seq(sim_no), function(sim_iter) {
    immigrant_events <- data.frame(magnitude = ifelse(is.null(M), get_new_magnitude(model_type$hawkes_decay_type), M),
                                   time = 0)
    if (!is.null(model_type$hawkes_immigrant_type)) {
      if (!is.null(init_history)) stop('Cannot simulate Hawkes processes with both background rate and given history')
      immigrant_events <- generate_immigrant_event_series(par, model_type, Tmax)
    }
    if (!is.null(model_type$hawkes_decay_type) && is.null(init_history)) {
      cascades <- lapply(seq(nrow(immigrant_events)), function(i) {
        generate_series_no_background_rate(par[get_param_names(new_hawkes(model_type = model_type$hawkes_decay_type))],
                                                        model_type$hawkes_decay_type, Tmax = Tmax,
                                                        maxEvents = maxEvents, M = M,
                                                        history_init = immigrant_events[i, ], tol = tol)
      })
      cascade <- do.call(rbind, cascades)
      if (nrow(immigrant_events) > 1) cascade[order(cascade$time), ]
    } else if (!is.null(model_type$hawkes_decay_type) && !is.null(init_history)) {
      cascade <- generate_series_no_background_rate(par[get_param_names(new_hawkes(model_type = model_type$hawkes_decay_type))],
                                                    model_type$hawkes_decay_type, Tmax = Tmax,
                                                    maxEvents = maxEvents, M = M,
                                                    history_init = init_history, tol = tol)
    } else {
      cascade <- immigrant_events
    }
    cascade
  }, mc.cores = cores)

  if (!return_as_object) {
    data
  } else {
    new_hawkes(model_type = model_type, par = par, data = data, observation_time = Tmax)
  }
}
