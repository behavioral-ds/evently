# this script implements methods for Hawkes with PL kernel (unmarked, PL, and marked, mPL)

get_param_names.hawkes_PL <- function(model) {
  c('K', 'c', 'theta')
}

get_ampl_likelihood.hawkes_PL <- function(model) {
  paste(
    'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(K)',
    '+ log(sum {j in 1..i-1} ((time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100))',
    ')',
    '- K * sum {i in 1..L[cn]-1} (( (1 / c)^theta - ( 1 / (time[cn,L[cn]] - time[cn,i] + c))^theta )) / theta);'
  )
}

get_ampl_constraints.hawkes_PL <- function(model) {
  'subject to branching_factor: K * (1 / theta) * (1 / c)^theta <= 1;'
}


get_param_names.hawkes_mPL <- function(model) {
  c('K', 'beta', 'c', 'theta')
}

get_ampl_likelihood.hawkes_mPL <- function(model) {
  paste(
    'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(K)',
    '+ log(sum {j in 1..i-1} (magnitude[cn,j]^beta * (time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100))',
    ')',
    '- K * sum {i in 1..L[cn]-1} (magnitude[cn,i]^beta * ( (1 / c)^theta - ( 1 / (time[cn,L[cn]] - time[cn,i] + c))^theta )) / theta);'
  )
}

get_ampl_constraints.hawkes_mPL <- function(model) {
  'subject to branching_factor: ( K * 1.016 / (1.016-beta) ) * (1 / theta) * (1 / c)^theta <= 1;'
}

#' @export
get_branching_factor.hawkes_PL <- function(model) {
  model$par[['K']]* (1 / model$par[['theta']]) * (1 / model$par[['c']])^model$par[['theta']]
}

#' @export
get_branching_factor.hawkes_mPL <- function(model) {
  # assuming alpha = 2.016
  (model$par[['K']] * 1.016 / (1.016-model$par[['beta']]) ) * (1 / model$par[['theta']]) * (1 / model$par[['c']])^model$par[['theta']]
}

#' @export
get_a1.hawkes_PL <- function(model) {
  processed_data <- preprocess_data(model$data, model$observation_time)
  vapply(processed_data, function(history) {
    sum(1 / (model$par[['theta']] * ((history$time[nrow(history)] + model$par[['c']] - history$time[-nrow(history)]) ^ model$par[['theta']]))) * model$par[['K']]
  }, FUN.VALUE = NA_real_)
}

#' @export
get_a1.hawkes_mPL <- function(model) {
  processed_data <- preprocess_data(model$data, model$observation_time)
  vapply(processed_data, function(history) {
    sum((history$magnitude[-nrow(history)]) ^ model$par[['beta']] / (model$par[['theta']] * ((history$time[nrow(history)] + model$par[['c']] - history$time[-nrow(history)]) ^ model$par[['theta']]))) * model$par[['K']]
  }, FUN.VALUE = NA_real_)
}

#' @export
predict_final_popularity.hawkes_PL <- function(model) {
  NextMethod()
}

#' @export
predict_final_popularity.hawkes_mPL <- function(model) {
  NextMethod()
}
