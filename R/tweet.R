# this script hosts methods that are related to extracting data from raw tweet json objects
# crawled from Twitter API

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE

#' This function extracts cascades from a given jsonl file where each line is a tweet
#' json object. Please refer to the Twitter developer documentation:
#' https://developer.twitter.com/en/docs/tweets/data-dictionary/overview/tweet-object
#' @param path File path to the tweets jsonl file
#' @param batch Number of tweets to be read for processing at each iteration, choose
#' the best number for your memory load. Defaults to at most 10000 tweets each iteration.
#' @param cores Number of cores to be used for processing each batch in parallel.
#' @param output_path If provided, the index.csv and data.csv files which define the cascaddes
#' will be generated. In index.csv, each row is a cascade where events can be obtained from data.csv
#' by corresponding indics (start_ind to end_ind). Defaults to NULL.
#' @param keep_user Twitter user ids will be kept.
#' @param keep_absolute_time Keep the absolute tweeting times.
#' @param keep_text Keep the tweet text.
#' @param progress The progress will be reported if set to True (default)
#' @param return_as_list If true then a list of cascades (data.frames) will be returned.
#' @param save_temp If temporary files should be generated while processing. Processing can be resumed on failures.
#' @param api_version Version of Twitter API used for collecting the tweets.
#' @return If return_as_list is TRUE then a list of data.frames where each data.frame is a retweet cascade.
#' Otherwise there will be no return.
#' @import parallel
#' @export
parse_raw_tweets_to_cascades <- function(path, batch = 100000, cores = 1, output_path = NULL,
                                         keep_user = F, keep_absolute_time = F, keep_text = F, progress = T,
                                         return_as_list = T, save_temp = F, api_version=1) {
  check_required_packages(c('jsonlite', 'data.table', 'bit64'))
  library(data.table)
  # a helper function
  zero_if_null <- function(count) {
    ifelse(is.null(count), 0, count)
  }

  if (api_version == 1) {
    parse_tweet <- function(tweet, keep_text = F) {
      tryCatch({
        json_tweet <- jsonlite::fromJSON(tweet)
        id <- json_tweet$id_str
        magnitude <- zero_if_null(json_tweet$user$followers_count)
        user_id <- json_tweet$user$id_str
        screen_name <- json_tweet$user$screen_name
        retweet_id <- NA
        if (keep_text) text <- json_tweet$text
        if (!is.null(json_tweet[['retweeted_status']])) {
          # if this tweet is a retweet, get original tweet's information
          retweet_id <- json_tweet$retweeted_status$id_str
          if (keep_text) text <- NA
        }
        res <- list(id = id, magnitude = magnitude, user_id = user_id,
             screen_name = screen_name, retweet_id = retweet_id)
        if (keep_text) res[['text']] <- text
        res
      },
      error = function(e) {
        warning(sprintf('Error processing json: %s', e))
      })
    }
  } else if (api_version == 2) {
    parse_tweet <- function(tweet, keep_text = F) {
      tryCatch({
        json_tweet <- jsonlite::fromJSON(tweet)
        if (is.null(json_tweet$includes) || is.null(json_tweet$includes$users)) {
          stop('The author information is required!')
        }
        id <- json_tweet$data$id
        magnitude <- zero_if_null(json_tweet$includes$users$public_metrics$followers_count)
        user_id <- json_tweet$data$author_id
        screen_name <- json_tweet$user$screen_name
        retweet_id <- NA
        if (keep_text) text <- json_tweet$data$text
        if (!is.null(json_tweet$data$referenced_tweets) && json_tweet$data$referenced_tweets$type == 'retweeted') {
          # if this tweet is a retweet, get original tweet's information
          retweet_id <- json_tweet$data$referenced_tweets$id
          if (keep_text) text <- NA
        }
        res <- list(id = id, magnitude = magnitude, user_id = user_id,
                    screen_name = screen_name, retweet_id = retweet_id)
        if (keep_text) res[['text']] <- text
        res
      },
      error = function(e) {
        warning(sprintf('Error processing json: %s', e))
      })
    }
  } else {
    stop('Unknown API version!')
  }

  con <- file(path, "r")
  i <- 1
  total_tweets <- 0
  processed_tweets_batch <- list()
  temp_prefix <- 'processed_tweets_tmp_'
  repeat {
    tweets <- readLines(con, n = batch)
    if (progress) cat(sprintf('Total tweets processed so far: %s; tweets to process at this iteration: %s',
                      total_tweets, length(tweets)))
    total_tweets <- total_tweets + length(tweets)
    if (length(tweets) == 0) break()

    if (save_temp && !file.exists(file.path(output_path, sprintf('%s%s.csv', temp_prefix, i)))) {
      stopifnot(!is.null(output_path))

      processed_tweets_batch_list <- data.table::rbindlist(mclapply(tweets, parse_tweet, keep_text = keep_text, mc.cores = cores))
      # save this intermediate results in case the function fails
      data.table::fwrite(processed_tweets_batch_list, file = file.path(output_path, sprintf('%s%s.csv', temp_prefix, i)))
      rm(processed_tweets_batch_list) # to clear up memory
    } else if (!save_temp) {
      processed_tweets_batch[[i]] <- data.table::rbindlist(mclapply(tweets, parse_tweet, keep_text = keep_text, mc.cores = cores))
    }
    cat('\r')
    rm(tweets) # to clear up memory
    i <- i + 1
  }
  cat('\n')
  close(con)

  if (save_temp) processed_tweets_batch <- lapply(list.files(path = output_path, pattern = temp_prefix),
                                                  function(f) data.table::fread(file.path(output_path, f)))
  processed_tweets <- data.table::as.data.table(data.table::rbindlist(processed_tweets_batch))
  processed_tweets[is.na(retweet_id), retweet_id := id]
  processed_tweets <- processed_tweets[(retweet_id %in% id) & !is.na(id)] #id could be NA due to processing errors
  processed_tweets[, absolute_time := melt_snowflake(id)$timestamp_ms]
  setorder(processed_tweets, retweet_id, absolute_time)
  processed_tweets[, time := absolute_time - absolute_time[1], retweet_id]
  processed_tweets[, time := as.double(time)/1000]
  processed_tweets[, index := 1:nrow(processed_tweets)]
  if (keep_text) {
    index <- processed_tweets[, .(start_ind = index[1],
                                  end_ind = index[length(index)],
                                  tweet_time = absolute_time[1]/1000,
                                  text = text[1]), retweet_id][, c('start_ind', 'end_ind', 'tweet_time', 'text')]
  } else {
    index <- processed_tweets[, .(start_ind = index[1],
                                  end_ind = index[length(index)],
                                  tweet_time = absolute_time[1]/1000), retweet_id][, c('start_ind', 'end_ind', 'tweet_time')]
  }
  processed_tweets[, absolute_time := absolute_time/1000]
  kept_columns <- c('time', 'magnitude')
  if (keep_user) kept_columns <- c(kept_columns, 'user_id', 'screen_name')
  if (keep_absolute_time) kept_columns <- c(kept_columns, 'absolute_time')

  data <- processed_tweets[, kept_columns, with = F]

  if (!is.null(output_path)) {
    fwrite(index, file = file.path(output_path, 'index.csv'))
    fwrite(data, file = file.path(output_path, 'data.csv'))

    tmp_files <- list.files(output_path, pattern = 'processed_tweets_tmp')
    file.remove(tmp_files)
  }

  if (return_as_list) {
    cascade_sizes <- index$end_ind - index$start_ind + 1
    # return as a list of datas
    datas <- split(data, rep(1:nrow(index), cascade_sizes))
    return(unname(datas))
  }
}
