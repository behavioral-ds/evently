# this script implements methods for plotting figures given fitted models

#' Plot a Hawkes process and its intensity function
#' @param model a model object where data, model_type and par are required
#' @param cascade_index determine which cascade in the list of cascades to plot
#' default to the first cascade
plot_event_series <- function(model, cascade_index = 1) {
  check_required_packages('ggplot2')
  check_required_hawkes_fields(model, c('model_type', 'data', 'par'))

  cascade <- model$data[[cascade_index]]
  g <- ggplot2::ggplot() +
    ggplot2::stat_function(data = data.frame(x = c(0, max(cascade$time) * 1.2)), aes(x), n = 1000,
                           fun = Vectorize(function(t) get_model_intensity_at(model, t = t)), color = 'red')
  built_g <- ggplot2::ggplot_build(g)
  max_intensity <- max(built_g$data[[1]]$y)
  g +
    ggplot2::geom_point(data = cascade, aes(x = time, y = magnitude / max(magnitude) * max_intensity /2 + max_intensity), size = 4) +
    ggplot2::geom_segment(data = cascade, aes(x = time, y = magnitude / max(magnitude) * max_intensity/2 + max_intensity, xend = time, yend = 0)) +
    xlab('time') + ylab('intensity')
}
