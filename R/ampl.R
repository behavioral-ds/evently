
#' @export
setup_ampl <- function(ampl_path, ipopt_path) {
  .globals$execution <- generate_ampl_execution(paste0(ampl_path, ';', ipopt_path))
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

#' Runs AMPL, with the given model (mod) and data file (dat). "solver" choices
#' the solver (default "minos"). Returns the fitted parameters, the value of the
#' log.likelihood function and the exit status.
.run <- function(tmp_files, model, solver = "minos") {
  ## need to populate other solvers here: "minor", "baron", "ipopt"
  model <- paste('model ', tmp_files$mod, ';', sep = '')
  data <- paste('data ', tmp_files$dat, ';', sep = '')

  var_names <- paste(names(model$init_par), collapse = ", ")
  solver_text <- sprintf("option solver %s;", solver)
  if (solver == "ipopt")
    solver_text <- paste(solver_text, "options ipopt_options 'linear_solver=ma57 print_level=1 max_iter=1000';", sep = '\n')
  print("halt_on_ampl_error=yes causing errors; not sure why but removed for now")
  if (solver == "knitro")
    solver_text <- paste(solver_text, "options knitro_options 'ms_enable=1 honorbnds=1';", sep = '\n')
  content <- paste(
    sprintf("option TMPDIR \"%s\";", .globals$tmp),
    solver_text,
    model,
    data,
    "solve;",
    sprintf("display %s, Likelihood, solve_exitcode > %s/res-%s.txt;", var_names, .globals$tmp, tmp_files$pid),
    "exit;",
    sep = "\n"
  )

  write(content, tmp_files$run)
  ret <- NULL
  ## run AMPL with the configs we created
  tryCatch(expr = {
    system(paste(.globals$execution, tmp_files$run))

    tmp <- read.csv(tmp_files$res, sep = '=', header = FALSE)
    ret <- tmp[,2]
    names(ret) <- trimws(levels(tmp[,1])[tmp[,1]])
  }, error = function(e) {
    warning(sprintf("AMPL execution failed! error %s: ", e$message))
  })
  if (!is.null(ret)) {
    ## construct the return list
    print("return neg.log.likelihood")
    model$value <- -ret[["Likelihood"]]

    ret <- ret[!names(ret) %in% "Likelihood"]

    model$convergence <- ret[["solve_exitcode"]]

    ret <- ret[!names(ret) %in% "solve_exitcode"]

    model$par <- ret
  }
  return(model)
}

output_mod <- function(model, file) {
  ## if HawkesN, lower bound for population is at least as many as I have seen
  if ("N" %in% names(model$init_par) && "N" %in% names(model$lower_bound) && model$lower_bound[["N"]] < nrow(model$data))
    model$lower_bound[["N"]] <- nrow(model$data) - 1

  ## correct initial parameters out of bounds -- should not happen, but ...
  model$init_par[is.finite(model$init_par) & (model$init_par > model$upper_bound)] <- model$upper_bound[is.finite(model$init_par) & (model$init_par > model$upper_bound)]
  model$init_par[is.finite(model$init_par) & (model$init_par < model$lower_bound)] <- model$lower_bound[is.finite(model$init_par) & (model$init_par < model$lower_bound)]

  ## first describe the data
  output.PARAM <- get_mod_param_def(model$data)

  ## next describe our initial parameters, if we have them
  ## assume that the variables we use are the ones in the init_params
  output.VAR <- "#define parameters"
  for (var in names(model$init_par)) {
    if (is.finite(model$init_par[var])) {
      output.VAR <- paste(output.VAR, paste('var ', var, ' := ', model$init_par[var],'; ', sep = ""), sep = '\n')
    }  else {
      output.VAR <- paste(output.VAR, paste('var ', var, '; ', sep = ""), sep = '\n')
    }
  }

  ## next, do the part of the objective function and kernel dependent constraints
  output.OJ <- output_mod_likelihood(model)

  ## finally, construct boxes (bounds) on parameters
  output.BOX <- ""
  for (var in names(model$init_par)) {
    ## construct bound statement
    crt <- sprintf('subject to %s_limit:', var)
    if (is.finite(model$lower_bound[var])) crt <- paste(crt, model$lower_bound[var], "<=", sep = " ")
    crt <- paste(crt, var, sep = " ")
    if (is.finite(model$upper_bound[var])) crt <- paste(crt, "<=", model$upper_bound[var], sep = " ")
    crt <- paste(crt, ";", sep = "")
    if (is.finite(model$lower_bound[var]) && is.finite(model$upper_bound[var]))
      output.BOX <- paste(output.BOX, crt, sep = '\n')
  }

  ## finally peace it all together
  output <- paste(output.PARAM, output.VAR, output.OJ, output.BOX, sep = '\n')
  write(output, file = file)
}

