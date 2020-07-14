# this script implements methods for Hawkes with EXP kernel (unmarked, EXP and marked, mEXP)

get_param_names.hawkes_EXP <- function(model) {
  c('K', 'theta')
}

get_ampl_likelihood.hawkes_EXP <- function(model) {
  paste('sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(K * theta * sum {j in 1..i-1} (exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
    '))',
    '- K * sum {i in 1..L[cn]-1} ((1 - exp(-1 * theta * (time[cn,L[cn]] - time[cn,i])))));')
}

get_ampl_constraints.hawkes_EXP <- function(model) {
  'subject to branching_factor: K <= 1;'
}

get_param_names.hawkes_mEXP <- function(model) {
  c('K', 'beta', 'theta')
}

get_ampl_likelihood.hawkes_mEXP <- function(model) {
  paste('sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (',
        'log(K * theta * sum {j in 1..i-1} (magnitude[cn,j]^beta * exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
        '))',
        '- K * sum {i in 1..L[cn]-1} (magnitude[cn,i]^beta * (1 - exp(-1 * theta * (time[cn,L[cn]] - time[cn,i])))));')
}

get_ampl_constraints.hawkes_mEXP <- function(model) {
  'subject to branching_factor: K * 1.016 + beta <= 1.016;'
}

#' @export
get_branching_factor.hawkes_EXP <- function(model) {
  model$par[['K']]
}

#' @export
get_branching_factor.hawkes_mEXP <- function(model) {
  # assuming alpha = 2.016
  (model$par[['K']] * 1.016) / (1.016 - model$par[['beta']])
}

#' @export
get_a1.hawkes_EXP <- function(model) {
  processed_data <- preprocess_data(model$data, model$observation_time)
  vapply(processed_data, function(history) {
    sum(1 / (exp((history$time[nrow(history)] - history$time[-nrow(history)]) * model$par[['theta']]))) * model$par[['K']]
  }, FUN.VALUE = NA_real_)
}

#' @export
get_a1.hawkes_mEXP <- function(model) {
  processed_data <- preprocess_data(model$data, model$observation_time)
  vapply(processed_data, function(history) {
    sum((history$magnitude[-nrow(history)]) ^ model$par[['beta']] / (exp((history$time[nrow(history)] - history$time[-nrow(history)]) * model$par[['theta']]))) * model$par[['K']]
  }, FUN.VALUE = NA_real_)
}

#' @export
predict_final_popularity.hawkes_EXP <- function(model, data = NULL, observation_time = NULL) {
  NextMethod()
}

#' @export
predict_final_popularity.hawkes_mEXP <- function(model, data = NULL, observation_time = NULL) {
  NextMethod()
}

#' @export
get_model_intensity_at.hawkes_EXP <- function(model, t, cascade_index = 1) {
  get_model_intensity_at.hawkes_mEXP(list(model_type = 'mEXP', par = c(model$par, beta = 0), data = model$data),
                         t = t, cascade_index = cascade_index)
}

#' @export
get_model_intensity_at.hawkes_mEXP <- function(model, t, cascade_index = 1) {
  event <- model$data[[cascade_index]]
  event <- event[event$time <= t, ]
  par <- model$par
  mi <- event$magnitude
  ti <- event$time

  fun_f <- par[["K"]]

  # ro(m_i) part - the influence of the user of the event
  fun_ro <- (mi) ^ par[["beta"]]

  # psi(t, ti) part - the decaying / relaxation kernel
  fun_psi <- par[["theta"]] * (exp(-par[["theta"]] * (t - ti)))

  val <- fun_f * fun_ro * fun_psi

  sum(val)
}
