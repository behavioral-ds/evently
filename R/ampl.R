#' Set up the AMPL environment by downloading an AMPL demo version and the compiled
#' ipopt binary. Only supports UNIX compatible OSs.
#' @param ampl_path the path where the AMPL folder will be placed
#' @importFrom utils download.file
#' @importFrom utils untar
#' @importFrom utils unzip
#' @export
setup_ampl <- function(ampl_path) {
  machine_ident <- paste0(Sys.info()[['sysname']], Sys.info()[['machine']])
  if (missing(ampl_path) || ampl_path == '') ampl_path <- Sys.getenv('HOME')
  ampl_final_path <- paste0(ampl_path, '/ampl')

  # download AMPL
  ampl_suffix <- switch (machine_ident,
    'Linuxx86_64' = 'ampl.linux64',
    'Darwinx86_64' = 'ampl.macosx64',
    stop('Cannot infer the operating system your are using. Autodownload failed.')
  )
  ampl_download_url <- paste0('https://ampl.com/demo/', ampl_suffix, '.tgz')
  download.file(ampl_download_url, destfile = '/tmp/ampl.tgz')
  untar('/tmp/ampl.tgz', exdir = ampl_path)
  file.rename(from = paste0(ampl_path, '/', ampl_suffix),
              to = ampl_final_path)

  # download ipopt
  ipopt_download_url <- switch (machine_ident,
    'Linuxx86_64' = 'https://ampl.com/dl/open/ipopt/ipopt-linux64.zip',
    'Darwinx86_64' = 'https://ampl.com/dl/open/ipopt/ipopt-osx.zip',
    stop('Cannot infer the operating system your are using. Autodownload failed.')
  )
  download.file(ipopt_download_url, destfile = '/tmp/ipopt.zip')
  unzip(zipfile = '/tmp/ipopt.zip', exdir = ampl_final_path)
  Sys.chmod(paste0(ampl_final_path, '/ipopt'), '777', use_umask = FALSE)

  # add ampl to path environment
  write(paste0('PATH=', Sys.getenv('PATH'), ':', ampl_final_path), file = paste0(Sys.getenv('HOME'), '/.Renviron'), append = TRUE)

  .globals$execution <- paste0('export PATH=$PATH:', ampl_final_path, '; ampl')
}

#' Set up the folder for placing temporary files, defaults to /tmp
#' @param path a string of path to the folder where you would like to place temporary files
#' @export
set_tmp_folder <- function(path) {
  .globals$tmp <- path
}

# preparation code --------------------------------------------------------
prepare_tmp_files <- function() {
  pid <- Sys.getpid()
  dat <- sprintf('%s/tmp-%s.dat', .globals$tmp, pid)
  mod <- sprintf('%s/tmp-%s.mod', .globals$tmp, pid)
  run <- sprintf('%s/tmp-%s.run', .globals$tmp, pid)
  res <- sprintf('%s/res-%s.txt', .globals$tmp, pid)
  if (file.exists(dat)) file.remove(dat)
  if (file.exists(mod)) file.remove(mod)
  if (file.exists(run)) file.remove(run)
  if (file.exists(res)) file.remove(res)

  list(dat = dat, mod = mod, run = run, res = res, pid = pid)
}


# ampl main code ----------------------------------------------------------

# fit with ampl
ampl_run <- function(model = model, solver = 'ipopt', goal = 'fit') {
  tmp_files <- prepare_tmp_files()
  output_dat(model = model, file = tmp_files$dat)
  output_mod(model = model, file = tmp_files$mod)

  # what's the expected result
  res <- switch (goal,
    fit = .run(tmp_files = tmp_files, model = model, solver = solver),
    nll = .run_get_likelihood(tmp_files = tmp_files, model = model)
  )

  file.remove(tmp_files$dat, tmp_files$mod, tmp_files$run, tmp_files$res)

  return(res)
}

# Runs AMPL, with the given model (mod) and data file (dat). "solver" choices
# the solver (default "minos"). Returns the fitted parameters, the value of the
# log.likelihood function and the exit status.
#' @importFrom utils read.csv
.run <- function(tmp_files, model, solver = "minos") {
  var_names <- paste(names(model$init_par), collapse = ", ")
  solver_text <- sprintf("option solver %s;", solver)
  if (solver == "ipopt")
    solver_text <- paste(solver_text, "options ipopt_options 'print_level=1 max_iter=1000';", sep = '\n')
  # print("halt_on_ampl_error=yes causing errors; not sure why but removed for now")
  # linear_solver=ma57
  if (solver == "knitro")
    solver_text <- paste(solver_text, "options knitro_options 'ms_enable=1 honorbnds=1';", sep = '\n')
  content <- paste(
    sprintf("option TMPDIR \"%s\";", .globals$tmp),
    solver_text,
    paste('model ', tmp_files$mod, ';', sep = ''),
    paste('data ', tmp_files$dat, ';', sep = ''),
    "solve;",
    sprintf("display %s, Likelihood, solve_exitcode > %s/res-%s.txt;", var_names, .globals$tmp, tmp_files$pid),
    "exit;",
    sep = "\n"
  )

  write(content, tmp_files$run)
  ret <- NULL

  ## run AMPL with the configs we created
  tryCatch(expr = {
    system(paste(.globals$execution, tmp_files$run), ignore.stdout = T, ignore.stderr = F)

    tmp <- read.csv(tmp_files$res, sep = '=', header = FALSE)
    ret <- tmp[,2]
    names(ret) <- trimws(levels(tmp[,1])[tmp[,1]])
  }, error = function(e) {
    warning(sprintf("AMPL execution failed! error %s: ", e$message))
  })
  if (!is.null(ret)) {
    ## construct the return list
    model$value <- -ret[["Likelihood"]]

    ret <- ret[!names(ret) %in% "Likelihood"]

    model$convergence <- ret[["solve_exitcode"]]

    ret <- ret[!names(ret) %in% "solve_exitcode"]

    model$par <- ret
  }
  return(model)
}

.run_get_likelihood <- function(tmp_files, model) {
  content <- paste(
    paste('model ', tmp_files$mod, ';', sep = ''),
    paste('data ', tmp_files$dat, ';', sep = ''),
    # use display the whole expression instead of display Likelihood, as display
    # Likelihood seems load the whole processing model into memory which slows
    # down this extremely
    sprintf("display %s > %s/res-%s.txt;", gsub(';', '', get_ampl_likelihood(model)), .globals$tmp, tmp_files$pid),
    "exit;",
    sep = "\n"
  )
  run_tmp <- sprintf('%s/tmp-%s.run', .globals$tmp, tmp_files$pid)
  res_tmp <- sprintf('%s/res-%s.txt', .globals$tmp, tmp_files$pid)
  write(content, run_tmp)

  system(paste(.globals$execution, run_tmp))

  tmp <- readChar(res_tmp, file.info(res_tmp)$size)
  ret <- as.numeric(unlist(strsplit(tmp, '='))[2])
  names(ret) <- c('Likelihood')

  # return neg log likelihood
  return(-ret[["Likelihood"]])
}
