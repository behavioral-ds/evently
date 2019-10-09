context('Fit on cascades')
cascade <- data.frame(time = seq(10), magnitude = seq(10))

test_that('only list of dataframe is accepted', {
  expect_error(fit_series(cascade))
})
