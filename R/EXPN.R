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

