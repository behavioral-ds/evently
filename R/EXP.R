# this script implements methods for Hawkes with EXP kernel (unmarked, EXP and marked, mEXP)

get_param_names.hawkes_EXP <- function(model) {
  c('K', 'theta')
}

get_ampl_likelihood.hawkes_EXP <- function(model) {
  paste('sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (log(K) + log(theta)',
    '+ log(sum {j in 1..i-1} (exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
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
  paste('sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (log(K) + log(theta)',
        '+ log(sum {j in 1..i-1} (magnitude[cn,j]^beta * exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
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
predict_final_popularity.hawkes_EXP <- function(model) {
  NextMethod()
}

#' @export
predict_final_popularity.hawkes_mEXP <- function(model) {
  NextMethod()
}
