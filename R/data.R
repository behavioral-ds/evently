#' `#auspol` hash tagged retweet diffusions
#'
#' A dataset containing retweet event diffusions relating to tweets hash tagged with `auspol`
#'
#' @format A list of 3333 data frames with three columns:
#' \describe{
#'   \item{time}{Relative retweet times w.r.t the initial tweet}
#'   \item{magnitude}{The number of followers a corresponding Twitter user has}
#'   \item{user}{An anonymized Twitter user id who created the tweet/retweet}
#' }
"auspol"
