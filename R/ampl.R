#' Set up the AMPL environment by downloading an AMPL demo version and the compiled
#' ipopt binary. Only supports UNIX compatible OSs.
#' @param ampl_path The path where the AMPL folder will be placed
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
  download_file(ampl_download_url, destfile = '/tmp/ampl.tgz')
  untar('/tmp/ampl.tgz', exdir = ampl_final_path)
  if (length(list.files(ampl_final_path)) == 1) {
    # 2021/01/29: the file structure of ampl.tgz changed, fix it
    pre_dir <- file.path(ampl_final_path, list.files(ampl_final_path))
    file.rename(file.path(pre_dir, list.files(pre_dir)), file.path(ampl_final_path, list.files(pre_dir)))
  }
  # file.rename(from = paste0(ampl_path, '/', ampl_suffix),
  #            to = ampl_final_path)

  # download ipopt
  ipopt_download_url <- switch (machine_ident,
    'Linuxx86_64' = 'https://ampl.com/dl/open/ipopt/ipopt-linux64.zip',
    'Darwinx86_64' = 'https://ampl.com/dl/open/ipopt/ipopt-osx.zip',
    stop('Failed to guess the operating system you are using. Autodownload failed.')
  )
  download_file(ipopt_download_url, destfile = '/tmp/ipopt.zip')
  unzip(zipfile = '/tmp/ipopt.zip', exdir = ampl_final_path)
  Sys.chmod(paste0(ampl_final_path, '/ipopt'), '777', use_umask = FALSE)

  # add ampl to path environment
  write(paste0('AMPL_PATH=', ampl_final_path), file = paste0(Sys.getenv('HOME'), '/.Renviron'), append = TRUE)
  Sys.setenv(AMPL_PATH = ampl_final_path)
  .globals$execution <- paste0('export PATH=$PATH:', ampl_final_path, '; ampl')
}

# download file with error handling
#' @importFrom utils download.file
download_file <- function(url, destfile) {
  status <- 1
  tryCatch({
    status <- download.file(url, destfile = destfile)
  }, error = function(e) {
    print('Failed to download, retry with wget...')
  })
  if (status != 0) {
    # failed to download, retry with wget
    download.file(url, destfile = destfile, method = 'wget')
  }
}

#' Set up the folder for placing temporary files, defaults to /tmp
#' @param path A string of path to the folder where you would like to place temporary files
#' @export
set_tmp_folder <- function(path) {
  .globals$tmp <- path
}

# preparation code --------------------------------------------------------
#' Prepare the temporary auxilixry files for AMPL
#' @param type One of "mod" (AMPL model file), "dat" (AMPL data file),
#' "run" (AMPL run file) and "res" (file hosts AMPL runned output)
prepare_tmp_file <- function(type) {
  pid <- Sys.getpid()
  f <- sprintf('%s/tmp-%s.%s', .globals$tmp, pid, type)
  if (file.exists(f)) file.remove(f)
  f
}


# ampl main code ----------------------------------------------------------

# fit with ampl
ampl_run <- function(model = model, solver = 'ipopt', dat_file, mod_file, goal = 'fit',
                     debug = F, ampl_execution = .globals$execution) {
  if (missing(dat_file)) {
    dat_file <- output_dat(model)
    # if the file is created within this function then remove once finishes
    if (!debug) on.exit(file.remove(dat_file), add = T)
  }
  if (missing(mod_file)) {
    mod_file <- output_mod(model)
    # if the file is created within this function then remove once finishes
    if (!debug) on.exit(file.remove(mod_file), add = T)
  }
  # double check AMPL path before running
  if (length(ampl_execution) == 0) stop('The path to AMPL execution is empty!')

  tmp_run <- prepare_tmp_file(type = 'run')
  tmp_res <- prepare_tmp_file(type = 'res')
  tmp_files <- list(dat = dat_file, mod = mod_file, run = tmp_run, res = tmp_res)
  # what's the expected result
  res <- switch (goal,
    fit = .run(model = model, solver = solver, tmp_files = tmp_files, debug, ampl_execution),
    nll = .run_get_likelihood(model = model, tmp_files = tmp_files, debug, ampl_execution)
  )

  return(res)
}

# Runs AMPL, with the given model (mod) and data file (dat). "solver" choices
# the solver (default "minos"). Returns the fitted parameters, the value of the
# log.likelihood function and the exit status.
#' @importFrom utils read.csv
.run <- function(model, solver = "minos", tmp_files, debug, ampl_execution, ...) {
  arguments <- list(...) # allow future extensions

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
    get_ampl_execution_options(model),
    paste('model ', tmp_files$mod, ';', sep = ''),
    paste('data ', tmp_files$dat, ';', sep = ''),
    "solve;",
    sprintf("display %s, Likelihood, solve_exitcode > %s;", var_names, tmp_files$res),
    "exit;",
    sep = "\n"
  )

  write(content, tmp_files$run)
  ret <- NULL

  # allow to debug the ampl execution
  if (debug) {
    stop(sprintf('Debugging is on! Please find the following AMPL files for inspection:\nModel: %s\nData: %s\nExecution: %s\n',
                 tmp_files$mod, tmp_files$dat, tmp_files$run), call. = F)
  }
  ## run AMPL with the configs we created
  tryCatch(expr = {
    system(paste(ampl_execution, tmp_files$run), ignore.stdout = T, ignore.stderr = F)

    tmp <- read.csv(tmp_files$res, sep = '=', header = FALSE, stringsAsFactors = FALSE)
    ret <- tmp[,2]
    names(ret) <- trimws(tmp[,1])
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

  file.remove(tmp_files$run, tmp_files$res)
  return(model)
}

.run_get_likelihood <- function(model, tmp_files, debug, ampl_execution, ...) {
  arguments <- list(...) # allow future extensions

  content <- paste(
    get_ampl_execution_options(model),
    paste('model ', tmp_files$mod, ';', sep = ''),
    paste('data ', tmp_files$dat, ';', sep = ''),
    # use display the whole expression instead of display Likelihood, as display
    # Likelihood seems load the whole processing model into memory which slows
    # down this extremely
    sprintf("display %s > %s;", gsub(';', '', get_ampl_likelihood(model)), tmp_files$res),
    "exit;",
    sep = "\n"
  )

  write(content, tmp_files$run)

  # allow to debug the ampl execution
  if (debug) {
    stop(sprintf('Debugging is on! Please find the following AMPL files for inspection:\nModel: %s\nData: %s\nExecution: %s\n',
                 tmp_files$mod, tmp_files$dat, tmp_files$run), call. = F)
  }

  system(paste(ampl_execution, tmp_files$run))

  tmp <- readChar(tmp_files$res, file.info(tmp_files$res)$size)
  ret <- as.numeric(unlist(strsplit(tmp, '='))[2])
  names(ret) <- c('Likelihood')

  file.remove(tmp_files$res, tmp_files$run)
  # return neg log likelihood
  return(-ret[["Likelihood"]])
}
