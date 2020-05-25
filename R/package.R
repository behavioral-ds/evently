#' Fitting Hawkes processes with AMPL
#'
#' @description To learn more about evently, start with the vignettes:
#' `browseVignettes(package = "evently")`
#' @docType package
#' @name evently
NULL

.globals <- new.env(parent = emptyenv())
.globals$execution <- sprintf('export PATH=$PATH:%s; ampl', Sys.getenv('AMPL_PATH'))
.globals$tmp <- '/tmp'
utils::globalVariables('.data')
