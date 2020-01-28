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

  if (is.null(observation_time)) stop('Please specify an observation time!')
  if (observation_time <= 0) stop('Observation time must be greater than 0!')
  if (is.infinite(observation_time)) {
    # ampl doesn't recognize infinity so set to a large number
    # didn't set to machine max as it will cause error in AMPL
    observation_time <- 1e20
  }

  data <- lapply(data, function(hist) {
    if (observation_time < hist$time[nrow(hist)]) {
      stop('The provided observation time is smaller than the last observed event.')
    }
    new_row <- data.frame(time = observation_time, magnitude = 0)
    hist <- rbind(hist[, c('time', 'magnitude')], new_row)
    offset_same_time(hist)
  })

  # sanity check
  if (any(sapply(data, function(d) is.unsorted(d$time)))) stop('something went wrong..')

  data
}

output_dat <- function(model, file) {
  # preprocess the data here
  data <- preprocess_data(model$data, model$observation_time)

  file.create(file)
  lengthes <- sapply(data, nrow)

  output_string <- paste0("param HL:= ", length(data), ";\n param L := ",
                          paste(seq_along(lengthes), lengthes, collapse = ' '),
                          ";\nparam ML:= ", max(lengthes),
                          ";\nparam: magnitude time :=")

  tmp_string_array <- rep(NA, times = sum(lengthes) * 5)
  k <- 1
  for (i in seq_along(data)) {
    for (j in seq(lengthes[i])) {
      tmp_string_array[k] <- '\n'
      tmp_string_array[k + 1] <- as.character(i)
      tmp_string_array[k + 2] <- as.character(j)

      tmp_string_array[k + 3] <- as.character(data[[i]]$magnitude[j])
      tmp_string_array[k + 4] <- as.character(data[[i]]$time[j])

      k <- k + 5
    }
  }
  tmp_string_array <- paste0(tmp_string_array, collapse = ' ')
  output_string <- paste0(output_string, tmp_string_array, ";\n param J0 := ")

  for (i in seq_along(data)) {
    output_string <- paste(output_string, i, min(which(data[[i]]$time > 0)) - 1, sep = " ")
  }
  output_string <- paste(output_string, ";", sep = "")

  # allow for model-specific data output
  output_string <- paste0(output_string, get_ampl_data_output(model))

  write(output_string, file = file)
}

output_mod <- function(model, file) {
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
    get_ampl_model_output(model),
    '# define objective function to maximize',
    'maximize Likelihood:',
    get_ampl_likelihood(model),
    '',
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
  output <- paste(output.PARAM, output.VAR, output.OJ, output.BOX, sep = '\n')
  write(output, file = file)
}
