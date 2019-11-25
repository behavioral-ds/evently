context('Simulation')

test_that('max time should be smaller Tmax when set', {
  max_time <- 1
  sim <- generate_hawkes_event_series(par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = max_time)$data[[1]]
  expect_lte(sim$time[nrow(sim)], max_time)
})

test_that('simulation should accept model class object', {
  model <- new_hawkes_model(par = c(K = 2, theta = 1), model_type = 'EXP')
  expect_error(generate_hawkes_event_series(model = model, par = c(K = 2, theta = 1), model_type = 'EXP', Tmax = 1))
  model_output <- generate_hawkes_event_series(model = model, Tmax = 1)
  expect_s3_class(model_output, 'hawkes_model')
  expect_true(is.list(model_output$data))
  expect_s3_class(model_output$data[[1]], 'data.frame')
})
