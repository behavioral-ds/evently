# this script implements methods for HawkesN with EXP kernel (unmarked, EXPN, and marked, mEXPN)

get_param_names.hawkes_EXPN <- function(model) {
  c('K', 'theta', 'N')
}

get_ampl_likelihood.hawkes_EXPN <- function(model) {
  paste(
    'sum {cn in 1..HL} ((  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(N-i+1))',
    '- (L[cn] - J0[cn] - 1) * log(N) + (L[cn] - J0[cn] - 1) * log(K * theta) + sum{i in J0[cn]+1..L[cn]-1} (log(sum{j in 1..i-1} (exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100 )))',
    '- K * sum {i in 1..L[cn]-1} (sum{j in i..L[cn]-1}((N - j) / N * (exp(-1*theta*(time[cn,j] - time[cn,i])) - exp(-1*theta * (time[cn,j+1] - time[cn,i])))))));'
  )
}

get_ampl_constraints.hawkes_EXPN <- function(model) {
  '' # no constraint for HawkesN
}

get_param_names.hawkes_mEXPN <- function(model) {
  c('K', 'beta', 'theta', 'N')
}

get_ampl_likelihood.hawkes_mEXPN <- function(model) {
  paste(
    'sum {cn in 1..HL} ((  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(N-i+1))',
    '- (L[cn] - J0[cn] - 1) * log(N) + (L[cn] - J0[cn] - 1) * log(K * theta) + sum{i in J0[cn]+1..L[cn]-1} (log(sum{j in 1..i-1} (exp(beta * log(magnitude[cn,j] + 1e-100)) * exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100 )))',
    '- K * sum {i in 1..L[cn]-1} (exp(beta * log(magnitude[cn,i]  + 1e-100)) * sum{j in i..L[cn]-1}((N - j) / N * (exp(-1*theta*(time[cn,j] - time[cn,i])) - exp(-1*theta * (time[cn,j+1] - time[cn,i])))))));'
  )
}

get_ampl_constraints.hawkes_mEXPN <- function(model) {
  '' # no constraint for HawkesN
}

#' @export
get_branching_factor.hawkes_EXPN <- function(model) {
  model$par[['K']]
}

#' @export
get_branching_factor.hawkes_mEXPN <- function(model) {
  # assuming alpha = 2.016
  (model$par[['K']] * 1.016) / (1.016 - model$par[['beta']])
}

#' @export
get_a1.hawkes_EXPN <- function(model) {
  stop('This method does not exist for EXPN')
}

#' @export
get_a1.hawkes_mEXPN <- function(model) {
  stop('This method does not exist for mEXPN')
}

#' @export
predict_final_popularity.hawkes_EXPN <- function(model, data = NULL, observation_time = NULL) {
  stop('This method does not exist for EXPN')
}

#' @export
predict_final_popularity.hawkes_mEXPN <- function(model, data = NULL, observation_time = NULL) {
  stop('This method does not exist for mEXPN')
}

#' @export
get_model_intensity_at.hawkes_EXPN <- function(model, t, cascade_index = 1) {
  get_model_intensity_at.hawkes_mEXPN(list(model_type = 'mEXPN', par = c(model$par, beta = 0), data = model$data),
                         t = t, cascade_index = cascade_index)
}

#' @export
get_model_intensity_at.hawkes_mEXPN <- function(model, t, cascade_index = 1) {
  event <- model$data[[cascade_index]]
  event <- event[event$time <= t, ]
  par <- model$par
  mi <- event$magnitude
  ti <- event$time

  ## compute correponding Nt at the current time t
  Nt <- min(sum(ti <= t), par[["N"]])


  # f(p_j) part - virality of a diffusion. Constant for a given diffusion. Furthermore, discount for available events.
  fun_f <- par[["K"]] * (1 - Nt / par[["N"]])

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- mi ^ par[["beta"]]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- par[["theta"]] * (exp(-par[["theta"]] * (t - ti)))

  val <- fun_f * fun_ro * fun_psi

  sum(val)
}
