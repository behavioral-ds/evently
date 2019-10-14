# Function for sampling from the powerlaw distribution of user influence
# it is the equivalent of the Richter-Gutenberg distribution in the Helmstetter model
# the powerlaw distribution was determined from the twitter data, from the #retweets
# alpha = 2.016, xmin = 1. Draw n values
#' @import poweRlaw
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
kernelFct <- function(event, t, params = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T, model_type='PL') {
  switch(model_type,
         PL = .kernelFct.PL(event, t, params = params, alpha = alpha, mmin = mmin, inclusive = inclusive),
         EXP = .kernelFct.EXP(event, t, params = params, alpha = alpha, mmin = mmin, inclusive = inclusive),
         EXPN = .kernelFct.EXPN(event, t, params = params, alpha = alpha, mmin = mmin, inclusive = inclusive),
         PLN = .kernelFct.PLN(event, t, params = params, alpha = alpha, mmin = mmin, inclusive = inclusive),
         stop(sprintf("Unimplemented kernel option '%s'.", model_type)))
}

# calculates the influence of a given event at the givent time the event is a 2
# elements list (mi, ti). Inclusive parameter means that if an event is present
# at time t, then it contributes to the conditional intensity. If inclusive ==
# F, then events at time t are removed.
.kernelFct.PL <- function(event, t, params = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2), alpha = 2.016, mmin = 1, inclusive = T) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- params["K"]

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ params["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + params["c"])^(1+params["theta"])

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
.kernelFct.EXP <- function(event, t, params = c(K = 0.024, beta = 0.5, theta = 0.2), alpha = 2.016, mmin = 1, inclusive = T) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- params["K"]

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ params["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- params["theta"] * (exp(-params["theta"] * (t - ti)))

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
.kernelFct.EXPN <- function(event, t, params = c(K = 0.024, beta = 0.5, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  ## compute correponding Nt at the current time t
  if (inclusive) {
    Nt <- min(sum(ti <= t), params["N"])
  } else {
    Nt <- min(sum(ti < t), params["N"])
  }

  # f(p_j) part - virality of a diffusion. Constant for a given diffusion. Furthermore, discount for available events.
  fun_f <- params["K"] * (1 - Nt / params["N"])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ params["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- params["theta"] * (exp(-params["theta"] * (t - ti)))

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
.kernelFct.PLN <- function(event, t, params = c(K = 0.024, beta = 0.5, c = 0.001, theta = 0.2, N = 1000), alpha = 2.016, mmin = 1, inclusive = T) {
  params <- unlist(params)
  # the event has 2 components: (magnitude_i, time_i)
  mat_event = matrix(unlist(event), ncol = 2, byrow = F)
  mi = mat_event[,1]
  ti = mat_event[,2]

  ## compute correponding Nt at the current time t
  if (inclusive) {
    Nt <- min(sum(ti <= t), params["N"])
  } else {
    Nt <- min(sum(ti < t), params["N"])
  }

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- params["K"] * (1 - Nt / params["N"])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi / mmin) ^ params["beta"]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + params["c"])^(1+params["theta"])

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
#' @param params a named vector of model parameters, K, alpha, beta, mmin, c, theta - parameters of the Hawkes kernel
#' @param M - magnitude of the initial event (in case no initial history
#'   provided)
#' @param history_init An initial history can be provided (R structure obtained
#'   from a previous call of this function). Its purpose is to allow custom
#'   initializations and to continue simulation of stopped processes.
#' @param Tmax - maximum time of simulation.
#' @param filename - file to which save the CSV file with the simulation.
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @export
generate_Hawkes_event_series <- function(params, model_type, alpha = 2.016, mmin = 1, M=10000, Tmax = 10, filename = NULL, history_init = NULL, maxEvents = NA, immigrants = NULL) {
  params <- unlist(params)
  # determine if this is a marked model or not
  if (substr(model_type, 1, 1) == 'm') {
    stopifnot('beta' %in% names(params))
    model_type <- substr(model_type, 2, nchar(model_type))
  } else {
    params[['beta']] <- 0
  }
  saveInterval <- 1000
  history <- NULL

  ## if no immigrant was given, then use the M at time 0
  if (is.null(immigrants)){
    immigrants <- as.data.frame(matrix(c(M, 0), nrow=1, byrow = TRUE))
    colnames(immigrants) <- c("magnitude", "time")
  }

  # initial event, magnitude M and time 0
  if ( !is.null(history_init)) {
    history <- history_init
    t <- history[nrow(history), "time"]
  }

  # maybe we have to continue simulation
  if (is.null(filename)) {
    if (is.null(history)) {
      history <- immigrants
      t <- history[1,]$time
      colnames(history) <- c("magnitude", "time")
    }
  } else {
    if (file.exists(filename)) {
      # means we are continuing the simulation
      # check first if there is some initialized history
      if ( !is.null(history)) {
        warning(sprintf("You gave me both a simulation to continue (yes, file %s exists) and an initial history. Ignoring initial history!", filename))
      }
      history <- read.table(file = filename, header = T)
      colnames(history) <- c("magnitude", "time")
      t <- history$time[nrow(history)]
      cat(sprintf("--> Loaded %d events from file, simulation time %.3f.\n", nrow(history), t))
    } else {
      # means we mearly want to save results to file
      # check first if there is some initialized history
      if (is.null(history)) {
        history <- immigrants
        t <- history$time[1]
        colnames(history) <- c("magnitude", "time")
      }
      cat(sprintf("--> Will save progress to history file %s, every %d events!\n", filename, saveInterval))
    }
  }

  while(t <= Tmax){
    # generate the time of the next event, based on the previous events
    t.lambda.max <- t
    intensityMax <- CIF(x = t.lambda.max, history = history, params = params, alpha = alpha, mmin = mmin, model_type = model_type)

    ## first, sample one event time from the maximum intensity
    r <- runif(1)
    t <- t - log(r) / intensityMax

    ## then perform rejection sampling
    s <- runif(1)
    thr <- CIF(x = t, history = history, params = params, alpha = alpha, mmin = mmin, model_type = model_type) / intensityMax
    if ( s <= thr) {
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

    # save progress
    if ((nrow(history) %% saveInterval == 0) && (!is.null(filename))) {
      write.table(x = history, file = filename, sep = "\t", row.names = F, col.names = T)
    }

    if (!is.na(maxEvents) && nrow(history) >= maxEvents)
      break;
  }

  # if a filename was provided, write the history to the file
  if (!is.null(filename)) {
    write.table(x = history, file = filename, sep = "\t", row.names = F, col.names = T)
  }

  # cat(sprintf("\n--> Simulation done!\n"))
  return(history)
}
