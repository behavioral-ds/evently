# *************************************************
#                     Setup
# *************************************************

#' @importFrom utils menu
.onAttach <- function(libname, pkgname) {
  if (Sys.getenv('AMPL_PATH') != '' || Sys.which('ampl') != '') {
    .globals$execution <- sprintf('export PATH=$PATH:%s; ampl', Sys.getenv('AMPL_PATH'))
  } else {
    packageStartupMessage("\n********************************************************")
    packageStartupMessage("  This package requires AMPL and ipopt")
    packageStartupMessage("  But they are not in your PATH environment.")
    packageStartupMessage("  Please specify its binary folder path in ~/.Renviron")
    packageStartupMessage("  Please also make sure the ipopt binary is in the ")
    packageStartupMessage("  same folder.")
    packageStartupMessage("********************************************************\n")
    packageStartupMessage('It seems AMPL is not found in your PATH environment, do you want to install it now?')
    if (interactive()) {
      installation_choice <- menu(c('yes', 'no'))
      if (installation_choice == 1) {
        ampl_path <- readline(prompt = paste0('Enter the path where AMPL should be place [', Sys.getenv('HOME'), ']: '))
        setup_ampl(ampl_path)
      } else {
        warning('AMPL is missing! You might not be able to fit models.')
      }
    } else {
      warning('AMPL is missing! You might not be able to fit models.')
    }
  }
}
