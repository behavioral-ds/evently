context('Extract cascades from raw tweet json objects')

test_that('Cascades can be extracted from raw tweets', {
  datas <- parse_raw_tweets_to_cascades(system.file('extdata', 'tweets_anonymized.jsonl', package = 'evently'))
  for (data in datas) {
    expect_s3_class(data, 'data.frame')
  }
  expect_s3_class(fit_series(datas, model_type = 'mPL', cores = 1, observation_time = Inf), 'hawkes')
})
