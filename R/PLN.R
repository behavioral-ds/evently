# this script implements methods for HawkesN with PL kernel (unmarked, PLN, and marked, mPLN)

get_param_names.hawkes_PLN <- function(model) {
  c('K', 'c', 'theta', 'N')
}

get_ampl_likelihood.hawkes_PLN <- function(model) {
  paste(
    'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(N-i+1) - log(N)',
    '+ log(K) + log(sum{j in 1..i-1} ((time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100 ))',
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
    '+ log(K) + log(sum{j in 1..i-1} (magnitude[cn,j]^beta * (time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100 ))',
    ')',
    '- K * sum {i in 1..L[cn]-1} ( magnitude[cn,i]^beta * sum{j in i..L[cn]-1}( (N - j) * ( (time[cn,j] - time[cn,i] + c)^(-1*theta) - (time[cn,j+1] - time[cn,i] + c)^(-1*theta) ))) / (theta * N));'
  )
}

get_ampl_constraints.hawkes_mPLN <- function(model) {
  ''
}
