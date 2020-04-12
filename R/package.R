#' Fitting Hawkes processes with AMPL
#'
#'
#' @docType package
#' @name evently
NULL

.globals <- new.env(parent = emptyenv())
.globals$execution <- if (Sys.which('ampl') != '') 'ampl' else NULL
.globals$tmp <- '/tmp'
