context('Extract diffusion features from cascades')

test_that('fit on a list of cascades individually', {
  data <- list(data.frame(time = c(0, 1), magnitude = c(1, 1)),
               data.frame(time = c(0, 2), magnitude = c(1, 1)))
  group_fits <- group_fit_series(data = data, model_type = 'PL')

  expect_less_than(abs(group_fits[[1]]$value - fit_series(data = data[1], model_type = 'PL')$value), 1e-3)
})
