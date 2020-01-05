library(testthat)
library(evently)
if (Sys.which('ampl') == '') setup_ampl(Sys.getenv('HOME'))

test_check("evently")
