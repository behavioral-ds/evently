#' @export
setup_ampl <- function(ampl_path, ipopt_bin_path, ipopt_lib_path) {
  .globals$execution <- generate_ampl_execution(paste0(ampl_path, ':', ipopt_bin_path, ':', ipopt_lib_path))
}

#' @export
set_tmp_folder <- function(path) {
  .globals$tmp <- path
}

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

# fit with ampl
ampl_run <- function(model = model, solver = 'ipopt', ...) {
  tmp_files <- prepare_tmp_files()
  output_dat(data = model$data, file = tmp_files$dat)
  output_mod(model = model, file = tmp_files$mod)

  res <- .run(tmp_files = tmp_files, model = model, solver = solver)
  file.remove(tmp_files$dat, tmp_files$mod, tmp_files$run, tmp_files$res)

  return(res)
}

ampl_get_neg_likelihood_value <- function(model, ...) {
  tmp_files <- prepare_tmp_files()
  output_dat(data = model$data, file = tmp_files$dat)
  output_mod(model = model, file = tmp_files$mod)

  neg_likelihood <- .run_get_likelihood(tmp_files = tmp_files, model = model)
  file.remove(tmp_files$dat, tmp_files$mod, tmp_files$run, tmp_files$res)

  return(neg_likelihood)
}

# Runs AMPL, with the given model (mod) and data file (dat). "solver" choices
# the solver (default "minos"). Returns the fitted parameters, the value of the
# log.likelihood function and the exit status.
#' @importFrom utils read.csv
.run <- function(tmp_files, model, solver = "minos") {
  var_names <- paste(names(model$init_par), collapse = ", ")
  solver_text <- sprintf("option solver %s;", solver)
  if (solver == "ipopt")
    solver_text <- paste(solver_text, "options ipopt_options 'linear_solver=ma57 print_level=1 max_iter=1000';", sep = '\n')
  # print("halt_on_ampl_error=yes causing errors; not sure why but removed for now")
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
    system(paste(.globals$execution, tmp_files$run), ignore.stdout = T, ignore.stderr = T)

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
  return(-ret["Likelihood"])
}
