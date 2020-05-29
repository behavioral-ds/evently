# this script implements methods for plotting figures given fitted models

#' Plot a Hawkes process and its intensity function
#' @param model A model object where data, model_type and par are required
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
plot_event_series <- function(model, cascade_index = 1) {
  check_required_packages('ggplot2')
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
    ggplot2::geom_segment(data = cascade,
                          ggplot2::aes(x = .data$time,
                                       y = .data$magnitude / max(.data$magnitude) * max_intensity/2 + max_intensity,
                                       xend = .data$time, yend = 0)) +
    ggplot2::xlab('time') + ggplot2::ylab('intensity')
}
