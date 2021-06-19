# this script implements methods for HawkesN with PL kernel (unmarked, PLN, and marked, mPLN)

get_param_names.hawkes_PLN <- function(model) {
  c('K', 'c', 'theta', 'N')
}

get_ampl_likelihood.hawkes_PLN <- function(model) {
  paste(
    'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(N-i+1) - log(N)',
    '+ log(K) + log(sum{j in ind[cn,i]..i-1} ((time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100 ))',
    ')',
    '- K * sum {i in 1..L[cn]-1} (sum{j in i..L[cn]-1}( (N - j) * ( (time[cn,j] - time[cn,i] + c)^(-1*theta) - (time[cn,j+1] - time[cn,i] + c)^(-1*theta) ))) / (theta * N));'
  )
}

get_ampl_constraints.hawkes_PLN <- function(model) {
  ''
}


get_param_names.hawkes_mPLN <- function(model) {
  c('K', 'beta', 'c', 'theta', 'N')
}

get_ampl_likelihood.hawkes_mPLN <- function(model) {
  paste(
    'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(N-i+1) - log(N)',
    '+ log(K) + log(sum{j in ind[cn,i]..i-1} (magnitude[cn,j]^beta * (time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100 ))',
    ')',
    '- K * sum {i in 1..L[cn]-1} ( magnitude[cn,i]^beta * sum{j in i..L[cn]-1}( (N - j) * ( (time[cn,j] - time[cn,i] + c)^(-1*theta) - (time[cn,j+1] - time[cn,i] + c)^(-1*theta) ))) / (theta * N));'
  )
}

get_ampl_constraints.hawkes_mPLN <- function(model) {
  ''
}

#' @export
get_branching_factor.hawkes_PLN <- function(model) {
  model$par[['K']]* (1 / model$par[['theta']]) * (1 / model$par[['c']])^model$par[['theta']]
}

#' @export
get_branching_factor.hawkes_mPLN <- function(model) {
  # assuming alpha = 2.016
  (model$par[['K']] * 1.016 / (1.016-model$par[['beta']]) ) * (1 / model$par[['theta']]) * (1 / model$par[['c']])^model$par[['theta']]
}

#' @export
get_a1.hawkes_PLN <- function(model) {
  stop('This method does not exist for PLN')
}

#' @export
get_a1.hawkes_mPLN <- function(model) {
  stop('This method does not exist for mPLN')
}

#' @export
predict_final_popularity.hawkes_PLN <- function(model, data = NULL, observation_time = NULL) {
  stop('This method does not exist for PLN')
}

#' @export
predict_final_popularity.hawkes_mPLN <- function(model, data = NULL, observation_time = NULL) {
  stop('This method does not exist for mPLN')
}

#' @export
get_model_intensity_at.hawkes_PLN <- function(model, t, cascade_index = 1) {
  get_model_intensity_at.hawkes_mPLN(list(model_type = 'mPLN', par = c(model$par, beta = 0), data = model$data),
                         t = t, cascade_index = cascade_index)
}

#' @export
get_model_intensity_at.hawkes_mPLN <- function(model, t, cascade_index = 1) {
  event <- model$data[[cascade_index]][model$data[[cascade_index]]$time <= t, ]
  par <- model$par
  mi <- event$magnitude
  ti <- event$time

  ## compute correponding Nt at the current time t
  Nt <- min(sum(ti <= t), par[["N"]])

  # f(p_j) part - virality of a video. Constant for a given video
  fun_f <- par[["K"]] * (1 - Nt / par[["N"]])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi) ^ par[["beta"]]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- 1 / (t - ti + par[["c"]])^(1+par[["theta"]])

  val <- fun_f * fun_ro * fun_psi

  sum(val)
}
