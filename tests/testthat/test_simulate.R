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
