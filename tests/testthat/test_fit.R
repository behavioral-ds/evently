context('Fit on cascades')
cascade <- data.frame(time = seq(10), magnitude = seq(10))

test_that('only list of dataframe is accepted', {
  expect_error(fit_series(cascade))
})

test_that('fitting works', {
  params <- c(K = 1.3, theta = 1, N = 100)
  sims <- generate_Hawkes_event_series(params = params, model_type = 'EXPN', M = 100)
  fitted <- fit_series(data = list(sims), model_type = 'EXP', observation_time = 1e10)
})

test_that('simulating and fitting are working for EXP') {
  set.seed(888)
  cut_time <- 7
  params <- c(K = 1.3, beta = 0.3, theta = 1, N = 100)
  sims <- lapply(seq(500), function(i) generate_Hawkes_event_series(params = params, model_type = 'mEXPN', M = 100, Tmax = cut_time))
  sims <- lapply(sims, function(.x) .x[.x$time < cut_time, ])
  fitted <- lapply(seq(10), function(i) {fit_series(sims[(50*(i-1) +1):(50*i)], cores = 10, model_type = 'mEXPN', lower_bound = c(K = 1e-100, beta = 1e-100, theta = 1e-100, N = 100), upper_bound = c(K = 1e+6, beta = 1.016 - 1e-100, theta = 1e+3, N = 100), observation_time = cut_time)})

  fitted_params <- do.call(rbind.data.frame, lapply(fitted, function(.x) as.list(.x[['par']])))
  expect_true(all(abs(apply(fitted_params, 2, mean) - params) < 1e-2))
}
