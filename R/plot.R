# this script implements methods for plotting figures given fitted models

#' Plot a Hawkes process and its intensity function
#' @param model A model object where data, model_type and par are required
#' @param cascade This cascade data.frame will be used if provided
#' @param cascade_index Determine which cascade in the list of cascades to plot
#' default to the first cascade
#' @return A ggplot object
#' @export
#' @examples
#' par <- c(K = 0.95, theta = 1)
#' data <- generate_series(model_type = 'EXP',
#'                                      par = par,
#'                                      sim_no = 1, Tmax = Inf)
#' plot_event_series(new_hawkes(model_type = 'EXP',
#'                              par = par,
#'                              data = data))
plot_event_series <- function(model, cascade = NULL, cascade_index = 1) {
  check_required_packages('ggplot2')
  if (!is.null(cascade)) {
    stopifnot(is.data.frame(cascade) || (is.list(cascade) && is.data.frame(cascade[[1]])))
    model$data <- if (is.data.frame(cascade)) list(cascade) else cascade
  }
  check_required_hawkes_fields(model, c('model_type', 'data', 'par'))

  cascade <- model$data[[cascade_index]]
  g <- ggplot2::ggplot() +
    ggplot2::stat_function(data = data.frame(x = c(0, max(cascade$time) * 1.2)), ggplot2::aes(.data$x), n = 1000,
                           fun = Vectorize(function(t) get_model_intensity_at(model, t = t)), color = 'red')
  built_g <- ggplot2::ggplot_build(g)
  max_intensity <- max(built_g$data[[1]]$y)
  g +
    ggplot2::geom_point(data = cascade,
                        ggplot2::aes(x = .data$time,
                                     y = .data$magnitude / max(.data$magnitude) * max_intensity /2 + max_intensity),
                        size = 4) +
    ggplot2::geom_segment(data = cascade, linetype = 2,
                          ggplot2::aes(x = .data$time,
                                       y = .data$magnitude / max(.data$magnitude) * max_intensity/2 + max_intensity,
                                       xend = .data$time, yend = 0)) +
    ggplot2::xlab('time') + ggplot2::ylab('intensity') +
    ggplot2::theme_bw()
}

#' Plot the kernel functions of Hawkes processes
#' @param fitted_models A list of fitted model objects to plot the kernel functions
#' @return A ggplot object
#' @export
plot_kernel_function <- function(fitted_models) {
  check_required_packages('ggplot2')
  library(ggplot2)
  cut_off_intensity <- 1e-3
  lapply(fitted_models, function(model) check_required_hawkes_fields(model, c('model_type', 'par')))
  data <- list(data.frame(magnitude = 1, time = 0)) # dummy data to get the kernel function from a Hawkes intensity function

  kernel_functions <- lapply(fitted_models, function(model) {
    model$data <- list(model$data[[1]][1,])
    Vectorize(function(t) get_model_intensity_at(model, t = t))
  })
  cut_off_time <- max(sapply(kernel_functions, function(f) {
    uniroot(f = function(x) f(x) - cut_off_intensity, interval = c(0, 1e6))$root
  }))
  g <- ggplot2::ggplot(data = data.frame(x = c(0, cut_off_time)), aes(x))
  for (f_i in seq_along(kernel_functions)) {
    g <- g + stat_function(fun = kernel_functions[[f_i]], aes(color = !!as.character(f_i)))
  }
  g + xlab('relative time') + ylab('kernel function value') + labs(color='model') + theme_bw()
}
