# some utility functions

check_required_packages <- function(pkg_name) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(sprintf("package %s is required for parse_raw_tweets_to_cascades(). \nPlease install.packages(\"%s\") to use this functionality.", pkg_name, pkg_name),
         call. = FALSE)
  }
}
