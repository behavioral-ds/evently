# this script implements methods for Hawkes with EXP kernel

get_param_names.hawkes_EXP <- function(model) {
  c('K', 'beta', 'theta')
}

get_likelihood.hawkes_EXP <- function(model) {
  paste('sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (log(K) + log(theta)',
    '+ log(sum {j in 1..i-1} (magnitude[cn,j]^beta * exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
    '))',
    '- K * sum {i in 1..L[cn]-1} (magnitude[cn,i]^beta * (1 - exp(-1 * theta * (time[cn,L[cn]] - time[cn,i])))));')
}

get_constrains.hawkes_EXP <- function(model) {
  'subject to branching_factor: K * 1.016 + beta <= 1.016;'
}
