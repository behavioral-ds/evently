# this script implements methods for Hawkes with multiple parts: a specific kernel function + a background rate type

get_param_names.hawkes_MULTI <- function(model) {
  decay_type <- model$model_type$hawkes_decay_type
  immigrant_type <- model$model_type$hawkes_immigrant_type

  c(get_param_names(new_hawkes(model_type = decay_type)),
    get_param_names(new_hawkes(model_type = immigrant_type)))
}

get_ampl_likelihood.hawkes_MULTI <- function(model) {
  combined_type <-
    paste(
      model$model_type$hawkes_decay_type,
      model$model_type$hawkes_immigrant_type,
      sep = '_'
    )
  switch (
    combined_type,
    'EXP_CONST' = paste(
      'sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (log(lambda + K * theta * sum {j in 1..i-1} (exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)))',
      '- lambda * (time[cn,L[cn]] - time[cn, 1]) - K * sum {i in 1..L[cn]-1} ((1 - exp(-1 * theta * (time[cn,L[cn]] - time[cn,i])))));'
    ),
    'PL_CONST' =   paste(
      'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(lambda + K*sum {j in 1..i-1} ((time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100)))',
      '- lambda * (time[cn,L[cn]] - time[cn, 1]) - K * sum {i in 1..L[cn]-1} (( (1 / c)^theta - ( 1 / (time[cn,L[cn]] - time[cn,i] + c))^theta )) / theta);'
    ),
    'mEXP_CONST' =   paste(
      'sum {cn in 1..HL} ( (1 - 0^(L[cn]-J0[cn]-1)) * (sum {i in J0[cn]+1..L[cn]-1} (log(lambda + K * theta * sum {j in 1..i-1} (magnitude[cn,j]^beta * exp(-1 * theta * (time[cn,i] - time[cn,j]))) + 1e-100)',
      '))',
      '- lambda * (time[cn,L[cn]] - time[cn, 1]) - K * sum {i in 1..L[cn]-1} (magnitude[cn,i]^beta * (1 - exp(-1 * theta * (time[cn,L[cn]] - time[cn,i])))));'
    ),
    'mPL_CONST' =   paste(
      'sum {cn in 1..HL} (  (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (log(lambda + K * sum {j in 1..i-1} (magnitude[cn,j]^beta * (time[cn,i] - time[cn,j] + c) ^ (-1-theta)) + 1e-100))',
      ')',
      '- lambda * (time[cn,L[cn]] - time[cn, 1]) - K * sum {i in 1..L[cn]-1} (magnitude[cn,i]^beta * ( (1 / c)^theta - ( 1 / (time[cn,L[cn]] - time[cn,i] + c))^theta )) / theta);'
    ),
    'EXPN_CONST' = stop('HawkesN models don\'t work with background intensities.'),
    'PLN_CONST' = stop('HawkesN models don\'t work with background intensities.'),
    'mEXPN_CONST' = stop('HawkesN models don\'t work with background intensities.'),
    'mPLN_CONST' = stop('HawkesN models don\'t work with background intensities.'),
    stop('Unknown model type!')
  )
}

get_ampl_constraints.hawkes_MULTI <- function(model) {
  decay_type <- model$model_type$hawkes_decay_type
  immigrant_type <- model$model_type$hawkes_immigrant_type

  c(get_ampl_constraints(new_hawkes(model_type = decay_type)),
    get_ampl_constraints(new_hawkes(model_type = immigrant_type)))
}

generate_random_points.hawkes_MULTI <- function(model) {
  decay_type <- model$model_type$hawkes_decay_type
  immigrant_type <- model$model_type$hawkes_immigrant_type

  cbind(generate_random_points(new_hawkes(model_type = decay_type)),
        generate_random_points(new_hawkes(model_type = immigrant_type)))
}

get_lower_bound.hawkes_MULTI <- function(model) {
  decay_type <- model$model_type$hawkes_decay_type
  immigrant_type <- model$model_type$hawkes_immigrant_type

  c(get_lower_bound(new_hawkes(model_type = decay_type)),
    get_lower_bound(new_hawkes(model_type = immigrant_type)))
}

get_upper_bound.hawkes_MULTI <- function(model) {
  decay_type <- model$model_type$hawkes_decay_type
  immigrant_type <- model$model_type$hawkes_immigrant_type

  c(get_upper_bound(new_hawkes(model_type = decay_type)),
    get_upper_bound(new_hawkes(model_type = immigrant_type)))
}
