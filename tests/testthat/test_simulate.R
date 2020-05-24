context('Simulation')

test_that('max time should be smaller Tmax when set', {
  max_time <- 1
  sim <- generate_hawkes_event_series(par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = max_time)[[1]]
  expect_lte(sim$time[nrow(sim)], max_time)
})

test_that('simulation should accept model class object', {
  model <- new_hawkes(par = c(K = 2, theta = 1), model_type = 'EXP')
  expect_error(generate_hawkes_event_series(model = model, par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = 1))
  data <- generate_hawkes_event_series(model = model, Tmax = 1)
  expect_true(is.list(data))
  expect_s3_class(data[[1]], 'data.frame')
  expect_equal(data[[1]]$magnitude[1], 1)
})

test_that('simulation works', {
  # EXP
  set.seed(888)
  par <- c(K = 0.9, theta = 1)
  sims <- generate_hawkes_event_series(par = par, model_type = 'EXP', Tmax = Inf, sim_no = 1000)
  sizes <- sapply(sims, nrow)
  expect_lt(abs(mean(sizes) - 1/(1-get_branching_factor(new_hawkes(par = par, model_type = 'EXP')))), 0.2)
})

test_that('compute intensity values correctly', {
  cascade <- data.frame(time = c(0, 1), magnitude = c(1, 2))
  model1 <- new_hawkes(model_type = 'EXP', par = c(K = 0.9, theta = 1), data = cascade)
  v1 <- get_model_intensity_at(model1, t = 3)
  model2 <- new_hawkes(model_type = 'mEXP', par = c(K = 0.9, beta = 0.5, theta = 1), data = cascade)
  v2 <- get_model_intensity_at(model2, t = 3)
  model3 <- new_hawkes(model_type = 'PL', par = c(K = 0.9, theta = 1, c = 2), data = cascade)
  v3 <- get_model_intensity_at(model3, t = 3)
  model4 <- new_hawkes(model_type = 'mPL', par = c(K = 0.9, beta = 0.5, theta = 1, c = 2), data = cascade)
  v4 <- get_model_intensity_at(model4, t = 3)
  model5 <- new_hawkes(model_type = 'EXPN', par = c(K = 0.9, theta = 1, N = 10), data = cascade)
  v5 <- get_model_intensity_at(model5, t = 3)
  model6 <- new_hawkes(model_type = 'mEXPN', par = c(K = 0.9, beta = 0.5, theta = 1, N = 10), data = cascade)
  v6 <- get_model_intensity_at(model6, t = 3)
  model7 <- new_hawkes(model_type = 'PLN', par = c(K = 0.9, theta = 1, c = 2, N = 10), data = cascade)
  v7 <- get_model_intensity_at(model7, t = 3)
  model8 <- new_hawkes(model_type = 'mPLN', par = c(K = 0.9, beta = 0.5, theta = 1, c = 2, N = 10), data = cascade)
  v8 <- get_model_intensity_at(model8, t = 3)

  expect_lt(abs(v1 - 0.1666101), 1e-4)
  expect_lt(abs(v2 - 0.2170621), 1e-4)
  expect_lt(abs(v3 - 0.09225), 1e-4)
  expect_lt(abs(v4 - 0.1155495), 1e-4)

  expect_lt(abs(v5 - 0.1332881), 1e-4)
  expect_lt(abs(v6 - 0.1736496), 1e-4)
  expect_lt(abs(v7 - 0.0738), 1e-4)
  expect_lt(abs(v8 - 0.09243961), 1e-4)
})

test_that('simulation with immgrants works', {
  # test CONST
  set.seed(888)
  par <- c(K = 0, theta = 1, lambda = 0.5)
  sims <- generate_hawkes_event_series(par = par, model_type = c('CONST', 'EXP'), Tmax = 10, sim_no = 1000)
  sizes <- sapply(sims, nrow)
  expect_lt(abs(mean(sizes) - 10*par[['lambda']] - 1), 0.2) # HPP of intensity 5

  set.seed(888)
  par <- c(K = 0, theta = 1, lambda = 0.5)
  sims <- generate_hawkes_event_series(par = par, model_type = c('CONST'), Tmax = 10, sim_no = 1000)
  sizes <- sapply(sims, nrow)
  expect_lt(abs(mean(sizes) - 10*par[['lambda']] - 1), 0.2) # HPP of intensity 5
})
