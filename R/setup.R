# *************************************************
#                     Setup
# *************************************************

generate_ampl_execution <- function(path) {
  paste0('export PATH=$PATH:', path, '; ampl')
}

.onAttach <- function(libname, pkgname) {
  if (any(Sys.which(c('ampl', 'ipopt')) == '')) {
    packageStartupMessage("\n********************************************************")
    packageStartupMessage("  This package works with AMPL and ipopt")
    packageStartupMessage("  But they are not in your PATH.")
    packageStartupMessage("  Please specify their binary folder paths with:")
    packageStartupMessage("  setup_ampl(ampl_path, ipopt_bin_path, ipopt_lib_path); set_tmp_folder(path);")
    packageStartupMessage("********************************************************\n")
  } else {
    .globals$execution <- generate_ampl_execution(Sys.getenv('PATH'))
  }
}
