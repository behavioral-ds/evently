# some utility functions

check_required_packages <- function(pkg_name) {
  pkg_name_to_install <- Filter(function(p) !requireNamespace(p, quietly = TRUE), pkg_name)
  if (length(pkg_name_to_install) > 0) {
    cat(sprintf("package(s) %s is required, do you want to install them now?", paste(pkg_name_to_install, collapse = ',')))
    if (interactive()) {
      installation_choice <- menu(c('yes', 'no'))
      if (installation_choice == 1) {
        install.packages(pkg_name_to_install)
      } else {
        stop(sprintf("package(s) %s is required. \nPlease run `install.packages(c(%s))` to use this functionality.",
                     paste(pkg_name_to_install, collapse = ','), paste(sprintf('"%s"', pkg_name_to_install), collapse = ', ')),
             call. = FALSE)
      }
    } else {
      stop(sprintf("package(s) %s is required. \nPlease run `install.packages(c(%s))` to use this functionality.",
                   paste(pkg_name_to_install, collapse = ','), paste(sprintf('"%s"', pkg_name_to_install), collapse = ', ')),
           call. = FALSE)
    }
  }
  for (pkg in pkg_name) suppressMessages(library(pkg, character.only = TRUE))
}
