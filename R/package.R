#' Fitting Hawkes processes with AMPL
#'
#'
#' @docType package
#' @name evently
NULL

.globals <- new.env(parent = emptyenv())
.globals$execution <- sprintf('export PATH=$PATH:%s; ampl', Sys.getenv('AMPL_PATH'))
.globals$tmp <- '/tmp'
