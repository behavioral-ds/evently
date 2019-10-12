# *************************************************
#                     Setup
# *************************************************

generate_ampl_execution <- function(path) {
  paste0('export PATH=$PATH:', path, '; ampl')
}

.onAttach <- function(libname, pkgname) {
  # a list of potential places to look for ampl and ipopt
  # currently only assumes a folder named ampl and a folder named coinipopt
  # in the parent directory
  potential_ampl <- c(paste(dirname(getwd()), 'ampl', sep = '/'))
  potential_ipopt <- c(paste(dirname(getwd()), 'coinipopt', sep = '/'),
                       paste(dirname(getwd()), 'CoinIpopt', sep = '/'))
  potential_ampl <- Filter(file.exists, potential_ampl)
  potential_ipopt <- Filter(file.exists, potential_ipopt)

  if (length(potential_ipopt) == 1 && length(potential_ampl) == 1) {
    setup_ampl(ampl_path = potential_ampl, ipopt_bin_path = paste0(potential_ipopt, '/build/bin'),
               ipopt_lib_path = paste0(potential_ipopt, '/build/lib'))
  } else {
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
}
