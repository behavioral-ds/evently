# This script hosts functions for extracting diffusion features from groups of cascades

#' Given a list of cascades, this function fits each cascade individually by calling [fit_series].
#' If the given cascades are in a named list, the names will be regarded as groups and the result will be reformatted as a list of
#' group fits.
#' @param data A (named) list of data.frame(s) where each data.frame is an event cascade with
#' event tims and event magnitudes (optional). The list names (if present) will be used for grouping cascades with same names.
#' @param model_type A string representing the model type, e.g. EXP for Hawkes processes with an exponential kernel function
#' @param observation_times A list or an vector of observation times for the cascades in data
#' @param cores The number of cores used for parallel fitting, defaults to 1 (non-parallel)
#' @param ... Check the available arguments of {fit_series}
#' @return A list of model obejcts where each object fits on an invidual cascade in data
#' @export
group_fit_series <- function(data, model_type, observation_times = NULL, cores = 1, ...) {
  stopifnot(is.list(data))
  stopifnot(length(data) == length(observation_times) || is.null(observation_times) || length(observation_times) == 1)
  if (length(observation_times) == 1) {
    observation_times <- rep(observation_times, length(data))
  }

  fits <- mclapply(seq_along(data), function(i) fit_series(data = data[[i]], model_type = model_type,
                                                           observation_time = observation_times[[i]],
                                                           cores = 1, ...), mc.cores = cores)
  # data is named then start grouping cascades
  if (!is.null(names(data))) {
    fits <- split(fits, f = names(data))
    for (i in seq_along(fits)) {
      class(fits[[i]]) <- 'hawkes.group.fits'
    }
  } else {
    class(fits) <- 'hawkes.group.fits'
  }

  fits
}

construct_temporal_features <- function(cascades) {
  sizes <- sapply(cascades, nrow)
  times <- sapply(cascades, function(.x) .x$time[nrow(.x)])
  intervals <- unlist(lapply(cascades, function(.x) diff(.x$time)))
  if (length(intervals) == 0) intervals <- -1

  magnitudes <- unlist(lapply(cascades, function(.x) .x$magnitude))
  rename <- function(x, name) {
    names(x) <- paste(name, names(x))
    x
  }
  c(rename(summary(sizes), 'size'), rename(summary(times), 'final event time'),
    rename(summary(intervals), 'event time interval'), rename(summary(magnitudes), 'user magnitude'),
    c(`number of cascades` = length(cascades)))
}

#' Given a list of group-fits produced by 'group_fit_series', this function generates features
#' for each group-fit by summarizing the fitted parameters.
#' @param list_fits A list of group fits returned by {group_fit_series}
#' @param data A named list of cascades
#' @return A data frame of features for each group. If features are all -1, it means all the
#' fits of the group are NAs
#' @export
generate_features <- function(list_fits, data = NULL) {
  # determine if list_fits is a list of hawkes.group.fits
  stopifnot(is.list(list_fits) && all(sapply(list_fits, function(fits) 'hawkes.group.fits' %in% class(fits))))
  conver_to_feature <- function(values, param) {
    summarized <- as.list(summary(values))
    names(summarized) <- paste(param, names(summarized))
    summarized
  }

  params <- get_param_names(list_fits[[1]][[1]])

  # # v1 simple summary
  # res <- lapply(list_fits, function(fits) {
  #   do.call(c, lapply(params, function(param) {
  #     param_values <- sapply(fits, function(single_fit) single_fit$par[[param]])
  #     param_values <- param_values[!is.na(param_values)]
  #     if (length(param_values) == 0) param_values <- -1 # assign -1 to all NAs
  #     conver_to_feature(param_values, param)
  #   }))
  # })
  # v2 discretization
  params_quantiles <- lapply(params, function(param) {
    all_param <- unlist(lapply(list_fits, function(fits) {
      sapply(fits, function(single_fit) single_fit$par[[param]])
    }))
    unique(quantile(all_param, seq(0, 1, by = 0.05), na.rm = T, names = F))
  })
  names(params_quantiles) <- params
  res <- lapply(list_fits, function(fits) {
    do.call(c, lapply(params, function(param) {
      param_values <- sapply(fits, function(single_fit) single_fit$par[[param]])
      param_values <- as.numeric(param_values[!is.na(param_values)])
      quantile_counts <- as.list(table(cut(param_values, breaks = params_quantiles[[param]], include.lowest = T)))
      names(quantile_counts) <- paste(param, seq(length(quantile_counts)))
      quantile_counts
    }))
  })
  res_df <- do.call(rbind.data.frame, res)
  res_df <- cbind(data.frame(id = rownames(res_df)), res_df)


  if (!is.null(data)) {
    data_features <- do.call(rbind.data.frame,
                             lapply(split(unname(data), names(data)),
                                 function(.x) as.list(construct_temporal_features(.x))))
    res_df <- cbind(res_df, data_features)
  }
  rownames(res_df) <- NULL
  res_df
}

# Compute order 1 wassersterin distance between the empirical distributions of v1 and v2
wasserstein1d <- function(v1, v2) {
  ecdf1 <- stats::ecdf(v1)
  ecdf2 <- stats::ecdf(v2)
  v <- sort(c(v1, v2))
  diff_v <- diff(v)
  sum(sapply(utils::head(v, n = -1), function(p) abs(ecdf1(p) - ecdf2(p))) * diff_v)
}

compute_fits_distance <- function(fits1, fits2) {
  stopifnot(class(fits1) == 'hawkes.group.fits' && class(fits2) == 'hawkes.group.fits')
  params <- names(fits1[[1]]$par)

  # compute distance for each variable and sum as the final distance
  sum(sapply(params, function(param) {
    param1 <- sapply(fits1, function(f) f$par[[param]])
    param2 <- sapply(fits2, function(f) f$par[[param]])
    wasserstein1d(param1, param2)
  }))
}

#' Given a list of grouped fits, compute a distance matrix
#' @param group_fits A list of grouped fits returned by {group_fit_series}
#' @return A dist matrix of pairwise distances between each group-fit
#' @export
fits_dist_matrix <- function(group_fits) {
  group_no <- length(group_fits)
  d <- apply(utils::combn(group_no, 2), 2, function(idx) {
    compute_fits_distance(group_fits[[idx[1]]], group_fits[[idx[2]]])
  })
  attr(d, 'Size') <- group_no
  if (!is.null(names(group_fits))) {
    attr(d, 'Labels') <- names(group_fits)
  }
  attr(d, 'Diag') <- FALSE
  attr(d, 'Upper') <- FALSE
  class(d) <- 'dist'
  d
}

#' @export
print.hawkes.group.fits <- function(x, ...) {
  cat(sprintf('Total fits: %s\n', length(x)))
  cat(sprintf('Model type: %s\n', x[[1]]$model_type))
}
