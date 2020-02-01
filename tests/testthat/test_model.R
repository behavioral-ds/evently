context('Model class')

test_that('model object shoud be created', {
  model <- new_hawkes(par = c(K = 1, theta = 0.3), model_type = 'EXP')
  expect_s3_class(model, c('hawkes_EXP'))
  expect_s3_class(model, c('hawkes'))
})

test_that('model object shoud not be created due to wrong parameters', {
  expect_warning(new_hawkes(par = c(1, 0.3), model_type = 'EXP'))
  expect_error(new_hawkes(par = c(K = 1, theta = 0.3, N = 100), model_type = 'EXP'))
})

test_that('branching factor is correctly computed', {
  model1 <- new_hawkes(par = c(K = 0.1, theta = 0.3), model_type = 'EXP')
  model2 <- new_hawkes(par = c(K = 0.1, beta = 1, theta = 0.3), model_type = 'mEXP')
  model3 <- new_hawkes(par = c(K = 0.1, c = 0.2, theta = 0.3), model_type = 'PL')
  model4 <- new_hawkes(par = c(K = 0.1, beta = 1, c = 0.2, theta = 0.3), model_type = 'mPL')

  model5 <- new_hawkes(par = c(K = 0.1, theta = 0.3, N = 100), model_type = 'EXPN')
  model6 <- new_hawkes(par = c(K = 0.1, beta = 1, theta = 0.3, N = 100), model_type = 'mEXPN')
  model7 <- new_hawkes(par = c(K = 0.1, c = 0.2, theta = 0.3, N = 100), model_type = 'PLN')
  model8 <- new_hawkes(par = c(K = 0.1, beta = 1, c = 0.2, theta = 0.3, N = 100), model_type = 'mPLN')

  expect_equal(get_branching_factor(model1), 0.1, tolerance = 1e-6)
  expect_equal(get_branching_factor(model2), 6.35, tolerance = 1e-6)
  expect_equal(get_branching_factor(model3), 0.5402189, tolerance = 1e-6)
  expect_equal(get_branching_factor(model4), 34.3039, tolerance = 1e-6)

  expect_equal(get_branching_factor(model5), 0.1, tolerance = 1e-6)
  expect_equal(get_branching_factor(model6), 6.35, tolerance = 1e-6)
  expect_equal(get_branching_factor(model7), 0.5402189, tolerance = 1e-6)
  expect_equal(get_branching_factor(model8), 34.3039, tolerance = 1e-6)
})
