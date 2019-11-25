context('Simulation')

test_that('max time should be smaller Tmax when set', {
  max_time <- 1
  sim <- generate_hawkes_event_series(par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = max_time)[[1]]
  expect_lte(sim$time[nrow(sim)], max_time)
})

test_that('simulation should accept model class object', {
  model <- new_hawkes_model(par = c(K = 2, theta = 1), model_type = 'EXP')
  expect_error(generate_hawkes_event_series(model = model, par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = 1))
  data <- generate_hawkes_event_series(model = model, Tmax = 1)
  expect_true(is.list(data))
  expect_s3_class(data[[1]], 'data.frame')
  expect_equal(data[[1]]$magnitude[1], 1)
})
