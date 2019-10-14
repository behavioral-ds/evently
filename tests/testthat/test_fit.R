context('Fit on cascades')
cascade <- data.frame(time = seq(10), magnitude = seq(10))

test_that('only list of dataframe is accepted', {
  expect_error(fit_series(cascade))
})

test_that('simulating and fitting are working for EXP') {
  cut_time <- 7
  params <- c(K = 1.3, beta = 0.3, theta = 1, N = 100)
  sims <- lapply(seq(200), function(i) generate_Hawkes_event_series(params = params, model_type = 'mEXPN', M = 100, Tmax = cut_time))
  sims <- lapply(sims, function(.x) .x[.x$time < cut_time, ])
  fitted <- lapply(seq(10), function(i) {fit_series(sims[(20*(i-1) +1):(20*i)], cores = 10, model_type = 'mEXPN', lower_bound = c(K = 1e-100, beta = 1e-100, theta = 1e-100, N = 100), upper_bound = c(K = 1e+6, beta = 1.016 - 1e-100, theta = 1e+3, N = 100), observation_time = cut_time)})

  fitted_params <- do.call(rbind.data.frame, lapply(fitted, function(.x) as.list(.x[['par']])))
  summary(fitted_params)
}
