context('Test plotting functions')

test_that('plot_kernel_function works', {
  model <- new_hawkes(par = c(K = 0.9, theta = 1), model_type = 'EXP')
  expect_s3_class(plot_kernel_function(fitted_models = list(model) ), 'ggplot')
})
