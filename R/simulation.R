# Function for sampling from the powerlaw distribution of user influence
# it is the equivalent of the Richter-Gutenberg distribution in the Helmstetter model
# the powerlaw distribution was determined from the twitter data, from the #retweets
# alpha = 2.016, xmin = 1. Draw n values
#' @import poweRlaw
#' @export
generate_user_influence <- function(n, alpha = 2.016, mmin = 1) {

  if ( is.null(.globals$user_infl) ) {
    .globals$user_infl <- conpl$new()
  }
  .globals$user_infl$setXmin(mmin)
  .globals$user_infl$setPars(alpha)

  return(dist_rand(.globals$user_infl, n))
}

##################################### SIMULATION OF THE HAWKES PROCESS ##############################################

# calculates the influence of a given event at the givent time. The event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
kernelFct <- function(event, t, par = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T, model_type='PL') {
  switch(model_type,
         PL = .kernelFct.PL(event, t, par = par, alpha = alpha, mmin = mmin, inclusive = inclusive),
         EXP = .kernelFct.EXP(event, t, par = par, alpha = alpha, mmin = mmin, inclusive = inclusive),
         EXPN = .kernelFct.EXPN(event, t, par = par, alpha = alpha, mmin = mmin, inclusive = inclusive),
         PLN = .kernelFct.PLN(event, t, par = par, alpha = alpha, mmin = mmin, inclusive = inclusive),
         stop(sprintf("Unimplemented kernel option '%s'.", model_type)))
}

# calculates the influence of a given event at the givent time the event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
.kernelFct.PL <- function(event, t, par = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2), alpha = 2.016, mmin = 1, inclusive = T) {
  par <- unlist(par)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- par["K"]

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ par["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + par["c"])^(1+par["theta"])

  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  if (!inclusive) {
    val[t == ti] = 0
    val[mi == mmin] = 0
  }

  (val)
}

# calculates the influence of a given event at the givent time the event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
.kernelFct.EXP <- function(event, t, par = c(K = 0.024, beta = 0.5, theta = 0.2), alpha = 2.016, mmin = 1, inclusive = T) {
  par <- unlist(par)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- par["K"]

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ par["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- par["theta"] * (exp(-par["theta"] * (t - ti)))

  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  if (!inclusive) {
    val[t == ti] = 0
    val[mi == mmin] = 0
  }

  (val)
}

# calculates the influence of a given event at the givent time the event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
.kernelFct.EXPN <- function(event, t, par = c(K = 0.024, beta = 0.5, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T) {
  par <- unlist(par)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  ## compute correponding Nt at the current time t
  if (inclusive) {
    Nt <- min(sum(ti <= t), par["N"])
  } else {
    Nt <- min(sum(ti < t), par["N"])
  }

  # f(p_j) part - virality of a diffusion. Constant for a given diffusion. Furthermore, discount for available events.
  fun_f <- par["K"] * (1 - Nt / par["N"])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ par["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- par["theta"] * (exp(-par["theta"] * (t - ti)))

  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  if (!inclusive) {
    val[t == ti] = 0
    val[mi == mmin] = 0
  }

  return(val)
}

# calculates the influence of a given event at the givent time the event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
.kernelFct.PLN <- function(event, t, par = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T) {
  par <- unlist(par)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  ## compute correponding Nt at the current time t
  if (inclusive) {
    Nt <- min(sum(ti <= t), par["N"])
  } else {
    Nt <- min(sum(ti < t), par["N"])
  }

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- par["K"] * (1 - Nt / par["N"])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ par["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + par["c"])^(1+par["theta"])

  val = fun_f * fun_ro * fun_psi
  val[t<ti] = 0
  val[mi<mmin] = 0
  if (!inclusive) {
    val[t == ti] = 0
    val[mi == mmin] = 0
  }

  (val)
}

## The CIF function is necessary to simulate a non-stationary (NON-HOMOGENUOUS) poisson process
# conditional intensity function - the CIF
# in here need to calculate the conditional intensity, by using the kernels
CIF = function(x, history, ...) {
  subst <- history[history$time <= x,]
  return(sum(kernelFct(event = subst, t = x, ...)))
}

#' Main function to generate a Hawkes process sequence. It allows intermediary
#' saves and continuing a stopped simulation. Creates a CSV file with two
#' columns, each row is an event: (magnitude, time)
#' @param model a model class object with par and model_type presented. par and
#' model_type are not requird once this is given
#' @param par a named vector of model parameters, K, alpha, beta, mmin, c, theta - parameters of the Hawkes kernel
#' @param model_type model type
#' @param sim_no the number of simulated cascades
#' @param cores the number of cores (processes) used for simulation
#' @param M - magnitude of the initial event (in case no initial history
#'   provided)
#' @param history_init An initial history can be provided (R structure obtained
#'   from a previous call of this function). Its purpose is to allow custom
#'   initializations and to continue simulation of stopped processes.
#' @param Tmax - maximum time of simulation.
#' @export
generate_hawkes_event_series <- function(model, par, model_type, sim_no = 1, cores = 1, alpha = 2.016, mmin = 1, M = 10000, Tmax = 10, history_init = NULL, maxEvents = NA) {
  if (!missing(model) && (!missing(par) || !missing(model_type))) {
    stop('Please either give model or (par, model_type) instead of both.')
  } else if (!missing(model)) {
    check_required_hawkes_model_fields(model, c('par', 'model_type'))
    par <- model$par
    model_type <- model$model_type
  } else {
    model <- new_hawkes_model(model_type = model_type, par = par)
  }

  # determine if this is a marked model or not
  if (substr(model_type, 1, 1) == 'm') {
    stopifnot('beta' %in% names(par))
    model_type <- substr(model_type, 2, nchar(model_type))
  } else {
    par[['beta']] <- 0
  }

  data <- mclapply(seq(sim_no), function(sim_iter) {
    # initial event, magnitude M and time 0
    history <- as.data.frame(matrix(c(M, 0), nrow=1, byrow = TRUE))
    colnames(history) <- c("magnitude", "time")
    t <- history[1,]$time

    if (!is.null(history_init)) {
      history <- history_init
      t <- history[nrow(history), "time"]
    }

    while(t <= Tmax){
      # generate the time of the next event, based on the previous events
      t.lambda.max <- t
      intensityMax <- CIF(x = t.lambda.max, history = history, par = par, alpha = alpha, mmin = mmin, model_type = model_type)
      ## if intensityMax is too small then cut
      if (intensityMax < 1e-5) break

      ## first, sample one event time from the maximum intensity
      r <- runif(1)
      t <- t - log(r) / intensityMax
      if (t > Tmax) break

      ## then perform rejection sampling
      s <- runif(1)

      thr <- CIF(x = t, history = history, par = par, alpha = alpha, mmin = mmin, model_type = model_type) / intensityMax

      if (s <= thr) {
        ## if here, it means we accept the event
        # generate the influence of the next event, by sampling the powerlaw distribution of the #retweets
        mag <- generate_user_influence(n = 1, alpha = alpha, mmin = mmin)

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

      if (!is.na(maxEvents) && nrow(history) >= maxEvents)
        break;
    }

    return(history)
  }, mc.cores = cores)

  model$data <- data

  return(model)
}
