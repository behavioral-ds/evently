context('Fit on cascades')

test_that('fit on a dataframe is allowed', {
  cascade <- data.frame(time = seq(10), magnitude = seq(10))
  expect_equal(fit_series(cascade, model_type = 'EXP')$convergence, 0)
})

test_that('data without magnitudes is allowed', {
  cascade <- data.frame(time = seq(10))
  expect_equal(fit_series(cascade, model_type = 'EXP')$convergence, 0)
})

test_that('Debug is possible', {
  cascade <- data.frame(time = seq(0, 10), magnitude = seq(11))
  expect_error(fit_series(list(cascade), model_type = 'EXP', observation_time = Inf, debug = T), regexp = 'Debugging is on!')
})

test_that('Errors for only single event cascades', {
  cascade <- data.frame(time = 0, magnitude = 1)
  expect_error(fit_series(list(cascade), model_type = 'EXP'), regexp = 'Please double check the observation time!')
})

test_that('Data sliced when observation time is smaller than last event time', {
  cascade <- data.frame(time = c(0, 5), magnitude = c(1, 1))
  expect_warning(fitted <- fit_series(list(cascade), model_type = 'EXP', observation_time = 1))
  expect_warning(get_hawkes_neg_likelihood_value(fitted))
})

test_that('fitting works', {
  set.seed(888)
  par <- c(K = 1.3, theta = 1, N = 100)
  par_data_frame <- as.data.frame(t(par))
  sims <- generate_series(par = par, model_type = 'EXPN')
  sims <- lapply(sims, function(.x) cbind(.x, data.frame(test = rep(1, nrow(.x)))))
  fitted <- fit_series(data = sims, model_type = 'EXPN', observation_time = Inf, init_pars = par_data_frame)
  expect_equal(fitted$convergence, 0)
  expect_equal(nrow(fitted$data[[1]]), nrow(sims[[1]]))
  expect_true(all(fitted$init_par == par))

  # test parallel fitting
  par_data_frame <- rbind(par_data_frame, par_data_frame)
  fitted <- fit_series(data = sims, model_type = 'EXPN', observation_time = Inf, cores = 2, init_pars = par_data_frame)
  expect_equal(fitted$convergence, 0)
  expect_equal(nrow(fitted$data[[1]]), nrow(sims[[1]]))
  expect_true(all(fitted$init_par == par))
  fitted <- fit_series(data = sims, model_type = 'EXPN', observation_time = Inf, cores = 2, init_pars = par_data_frame, parallel_type = 'FORK')
  expect_equal(fitted$convergence, 0)
  expect_equal(nrow(fitted$data[[1]]), nrow(sims[[1]]))
  expect_true(all(fitted$init_par == par))

  # a testing data.frame
  test_sim <- list(data.frame(time = c(0, 1, 1, 1, 2, 2, 2), magnitude = rep(1, 7)))
  fitted <- fit_series(test_sim, model_type = 'EXP', observation_time = 2)
  expect_equal(fitted$convergence, 0)
})

test_that('fitting works for PL', {
  set.seed(888)
  par <- c(K = 1.3, c = 1, theta = 1)
  sims <- generate_series(par = par, model_type = 'PL')
  fitted <- fit_series(data = sims, model_type = 'PL', observation_time = Inf)
  expect_equal(fitted$convergence, 0)
  expect_equal(nrow(fitted$data[[1]]), nrow(sims[[1]]))
})

test_that('simulating and fitting are working for EXPN', {
  set.seed(888)
  cut_time <- 7
  par <- c(K = 1.3, theta = 1, N = 100)
  sims <- generate_series(par = par, model_type = 'EXPN', Tmax = cut_time, sim_no = 500)
  sims <- lapply(sims, function(.x) .x[.x$time < cut_time, ])
  fitted <- lapply(seq(10), function(i) {fit_series(sims[(50*(i-1) +1):(50*i)], cores = 1, model_type = 'EXPN', lower_bound = c(K = 1e-100, theta = 1e-100, N = 100), upper_bound = c(K = 1e+6, theta = 1e+3, N = 100), observation_time = cut_time)})
  expect_true(fitted[[1]]$init_par[[1]] != fitted[[1]]$par[[1]])
  fitted_par <- do.call(rbind.data.frame, lapply(fitted, function(.x) as.list(.x[['par']])))
  expect_true(all(abs(apply(fitted_par, 2, mean) - par) < 1e-1))
})

test_that('infinity observation time works for EXPN', {
  set.seed(888)
  par <- c(K = 5, theta = 1, N = 50)
  sims <- generate_series(par = par, model_type = 'EXPN', Tmax = Inf, sim_no = 50)
  fitted <- fit_series(sims, cores = 1, model_type = 'EXPN', observation_time = Inf)
  expect_true(fitted$init_par[[1]] != fitted$par[[1]])
  expect_true(all(abs(fitted$par- par) < 5e-1))
})

test_that('random observation times works for EXPN', {
  set.seed(889)
  par <- c(K = 5, theta = 0.2, N = 50)
  observation_times <- runif(50, 5, 20)
  sims <- lapply(seq(50), function(i) generate_series(par = par, model_type = 'EXPN', Tmax = observation_times[i], sim_no = 1)[[1]])
  fitted <- fit_series(sims, cores = 1, model_type = 'EXPN', observation_time = observation_times)
  expect_true(fitted$init_par[[1]] != fitted$par[[1]])
  expect_true(all(abs(fitted$par- par) < 5e-1))
})


test_that('compute holdout log-likelihood works', {
  cut_time <- 7
  par <- c(K = 1.3, theta = 1)
  model <- new_hawkes(model_type = 'EXP', data = list(data.frame(time = 0, magnitude = 1)), observation_time = Inf, par = par)
  nll <- get_hawkes_neg_likelihood_value(model)
  expect_lt(abs(nll - par[['K']]), 1e-10)
  nll <- get_hawkes_neg_likelihood_value(par = model$par, data = model$data, model_type = model$model_type, observation_time = model$observation_time)
  expect_lt(abs(nll - par[['K']]), 1e-10)
  nll <- get_hawkes_neg_likelihood_value(model, observation_time = 1)
  expect_lt(abs(nll - par[['K']] * (1 - exp(-1 * par['theta']))), 1e-5)
  nll <- get_hawkes_neg_likelihood_value(par = model$par, data = model$data[[1]], model_type = model$model_type, observation_time = model$observation_time)
  expect_lt(abs(nll - par[['K']]), 1e-10)
})

test_that('fit with constant background rate works', {
  cascade <- data.frame(time = seq(0, 10), magnitude = seq(11))
  fitted <- fit_series(cascade, model_type = c('EXP', 'CONST'))
  expect_s3_class(fitted, 'hawkes')
})

test_that('fit with limit_event should work', {
  cascade <- data.frame(time = seq(0, 10), magnitude = seq(11))
  expect_error(fit_series(cascade, model_type = 'EXP', limit_event = list()))
  expect_error(fit_series(cascade, model_type = 'EXP', limit_event = list(type = 'wrongtype')))
})


test_that('fit DMM should work', {
  datas <- list(
    data.frame(time = seq(0, 10), magnitude = seq(11)),
    data.frame(time = 0, magnitude = 1),
    data.frame(time = seq(0, 10), magnitude = seq(11)),
    data.frame(time = seq(0, 10), magnitude = seq(11))
  )
  res <- fit_series(datas, model_type = 'DMM')
  expect_false(any(is.null(res$par)))
})
