# this script hosts methods that are related to extracting data from raw tweet json objects
# crawled from Twitter API

#' This function extracts cascades from a given jsonl file where each line is a tweet
#' json object. Please refer to the Twitter developer documentation:
#' https://developer.twitter.com/en/docs/tweets/data-dictionary/overview/tweet-object
#' @param path File path to the tweets jsonl file
#' @param keep_user Twitter user ids will be kept
#' @param keep_absolute_time Keep the absolute tweeting times
#' @param progress A progress bar will present if set to True (default)
#' @return A list of data.frames where each data.frame is a retweet cascade
#' @export
parse_raw_tweets_to_cascades <- function(path, keep_user = F, keep_absolute_time = F, progress = T) {
  check_required_packages('jsonlite')
  con <- file(path, "r")
  tweets <- readLines(con, n = -1)
  close(con)

  # use environments as dictionaries for fast lookup
  cascades_time <- new.env() # retweet times relative to the original tweet
  cascades_magnitude <- new.env() # the number of followers
  cascades_user <- new.env() # user id
  cascades_tweet_time <- new.env() # absolute time of the tweet

  # a helper function
  zero_if_null <- function(count) {
    ifelse(is.null(count), 0, count)
  }

  # compute the time difference between two time strings from tweet objects
  convert_time_epoch <- function(t) {
    as.numeric(strptime(t, "%a %b %d %T %z %Y", tz = 'GMT'))
  }

  if (progress) pb <- utils::txtProgressBar(min = 0, max = length(tweets), style = 3)
  for (k in seq_along(tweets)) {
    tweet <- tweets[[k]]
    if (progress) utils::setTxtProgressBar(pb, k)
    tryCatch({
      json_tweet <- jsonlite::fromJSON(tweet)
      current_id <- json_tweet$id_str
      current_time <- convert_time_epoch(json_tweet$created_at)

      if (!is.null(json_tweet[['retweeted_status']])) {
        # if this tweet is a retweet, get original tweet's information
        original_time <- convert_time_epoch(json_tweet$retweeted_status$created_at)
        original_id <- json_tweet$retweeted_status$id_str

        if (is.null(cascades_time[[original_id]])) {
          # if the original tweet hasn't been added yet, we add it here

          cascades_time[[original_id]] <- 0 # as this is relative time, the time for initial tweet is 0
          cascades_magnitude[[original_id]] <- zero_if_null(json_tweet$retweeted_status$user$followers_count)
          cascades_user[[original_id]] <- json_tweet$retweeted_status$user$id_str
          cascades_tweet_time[[original_id]] <- original_time
        }
        cascades_time[[original_id]][length(cascades_time[[original_id]]) + 1] <- current_time - original_time
        cascades_magnitude[[original_id]][length(cascades_magnitude[[original_id]]) + 1] <- zero_if_null(json_tweet$user$followers_count)
        cascades_user[[original_id]][length(cascades_user[[original_id]]) + 1] <- json_tweet$user$id_str
      } else {
        # if this tweet is not a retweet
        if (is.null(cascades_time[[current_id]])) {
          cascades_time[[current_id]] <- 0

          cascades_magnitude[[current_id]] <- zero_if_null(json_tweet$user$followers_count)
          cascades_user[[current_id]] <- json_tweet$user$id_str
          cascades_tweet_time[[current_id]] <- current_time
        } else {
          # this tweet might be added already via its retweets, but we use this magnitude instead as it's the status when this tweet happened
          cascades_magnitude[[current_id]][1] <- zero_if_null(json_tweet$user$followers_count)
        }
      }
    },
    error = function(e) {
      warning(sprintf('Error processing json: %s', e))
    })
  }
  if (progress) close(pb)

  # define vectors for final outputs
  res_time <- c()
  res_mag <- c()
  res_user <- c()
  res_start <- c()
  res_end <- c()
  res_tweet_time <- c()
  i <- 0

  # finally, let's extract all cascades
  for (cascade in names(cascades_time)) {
    tmp_time <- cascades_time[[cascade]]
    tmp_mg <- cascades_magnitude[[cascade]]
    tmp_user <- cascades_user[[cascade]]

    sorted_indices <- order(tmp_time)
    res_time <- c(res_time, tmp_time[sorted_indices])
    res_mag <- c(res_mag, tmp_mg[sorted_indices])
    res_user <- c(res_user, tmp_user[sorted_indices])

    res_start <- c(res_start, i + 1)
    res_end <- c(res_end, i + length(tmp_time))
    res_tweet_time <- c(res_tweet_time, cascades_tweet_time[[cascade]])
    i <- i + length(tmp_time)
  }

  # formatted as two dataframes in case we want to output as the usual csv formats
  index <- data.frame(start_ind = res_start, end_ind = res_end, tweet_time = res_tweet_time)
  data <- data.frame(magnitude = res_mag, time = res_time)
  cascade_sizes <- index$end_ind - index$start_ind + 1
  if (keep_user) data <- cbind(data, data.frame(user = res_user, stringsAsFactors = F))
  if (keep_absolute_time) {
    data <- cbind(data,
                  data.frame(absolute_time = rep(res_tweet_time, cascade_sizes) + data$time,
                             stringsAsFactors = F))
  }
  # return as a list of datas
  datas <- split(data, rep(1:nrow(index), cascade_sizes))
  return(unname(datas))
}
