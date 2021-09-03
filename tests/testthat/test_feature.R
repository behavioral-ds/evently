context('Extract diffusion features from cascades')

test_that('fit on a list of cascades individually', {
  data <- list(data.frame(time = c(0, 1), magnitude = c(1, 1)),
               data.frame(time = c(0, 2), magnitude = c(1, 1)))
  for (i in seq(5)) {
    data <- c(data, data)
  }
  data <- c(data, list(data.frame(time = c(0, 2), magnitude = c(1, 1))))
  names(data) <- c(rep('1', 32), rep('2', 32), '3')
  group_fits <- group_fit_series(data = data)

  expect_s3_class(group_fits, 'hawkes.group.fits')
})
