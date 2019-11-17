# *************************************************
#                     Setup
# *************************************************

.onAttach <- function(libname, pkgname) {
  if (Sys.which('ampl') != '') {
    .globals$execution <- 'ampl'
  } else {
    packageStartupMessage("\n********************************************************")
    packageStartupMessage("  This package requires AMPL and ipopt")
    packageStartupMessage("  But they are not in your PATH environment.")
    packageStartupMessage("  Please specify its binary folder path by executing:")
    packageStartupMessage("  setup_ampl(ampl_path);")
    packageStartupMessage("  Please also make sure the ipopt binary is in the same")
    packageStartupMessage("  folder and ma57 linear solver is installed.")
    packageStartupMessage("********************************************************\n")
  }
}
