context('Fit on cascades')

test_that('only list of dataframe is accepted', {
  cascade <- data.frame(time = seq(10), magnitude = seq(10))
  expect_error(fit_series(cascade))
})

test_that('fitting works', {
  set.seed(888)
  params <- c(K = 1.3, theta = 1, N = 100)
  sims <- generate_Hawkes_event_series(params = params, model_type = 'EXPN')
  fitted <- fit_series(data = list(sims), model_type = 'EXPN', observation_time = Inf)
  expect_equal(fitted$convergence, 0)
  expect_equal(nrow(fitted$data[[1]]), nrow(sims))
})

test_that('simulating and fitting are working for EXPN', {
  set.seed(888)
  cut_time <- 7
  params <- c(K = 1.3, theta = 1, N = 100)
  sims <- lapply(seq(500), function(i) generate_Hawkes_event_series(params = params, model_type = 'EXPN', Tmax = cut_time))
  sims <- lapply(sims, function(.x) .x[.x$time < cut_time, ])
  fitted <- lapply(seq(10), function(i) {fit_series(sims[(50*(i-1) +1):(50*i)], cores = 10, model_type = 'EXPN', lower_bound = c(K = 1e-100, theta = 1e-100, N = 100), upper_bound = c(K = 1e+6, theta = 1e+3, N = 100), observation_time = cut_time)})
  expect_true(fitted[[1]]$init_par[[1]] != fitted[[1]]$par[[1]])
  fitted_params <- do.call(rbind.data.frame, lapply(fitted, function(.x) as.list(.x[['par']])))
  expect_true(all(abs(apply(fitted_params, 2, mean) - params) < 1e-1))
})

test_that('infinity observation time works for EXPN', {
  set.seed(888)
  params <- c(K = 5, theta = 1, N = 50)
  sims <- lapply(seq(50), function(i) generate_Hawkes_event_series(params = params, model_type = 'EXPN', Tmax = Inf))
  fitted <- fit_series(sims, cores = 10, model_type = 'EXPN', observation_time = Inf)
  expect_true(fitted[[1]]$init_par[[1]] != fitted[[1]]$par[[1]])
  fitted_params <- do.call(rbind.data.frame, lapply(fitted, function(.x) as.list(.x[['par']])))
  expect_true(all(abs(apply(fitted_params, 2, mean) - params) < 5e-1))
})

test_that('compute holdout log-likelihood works', {
  cut_time <- 7
  params <- c(K = 1.3, theta = 1)
  model <- new_hawkes_model(model_type = 'EXP', data = list(data.frame(time = 0, magnitude = 1)), observation_time = Inf, par = params)
  nll <- get_hawkes_neg_likelihood_value(model)
  expect_less_than(abs(nll - params[['K']]), 1e-10)
})
