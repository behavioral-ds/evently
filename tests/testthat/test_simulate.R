context('Simulation')

test_that('max time should be smaller Tmax when set', {
  max_time <- 1
  sim <- generate_hawkes_event_series(params = c(K = 2, theta = 1), model_type = 'EXP', Tmax = max_time)
  expect_lte(sim$time[nrow(sim)], max_time)
})
