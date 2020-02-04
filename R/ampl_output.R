ampl_output_from_r <- function(names, var, type) {
  name <- paste(names, collapse = ' ')
  output_string_head <- paste0(ifelse(length(names) == 1, "param ", "param: "),
                               name, ":= ")
  # based on the class of given var
  if (type == 'atomic') {
    output_string_tail <- var
  } else if (type == 'vector') {
    output_string_tail <- paste(seq_along(var), var, collapse = ' ')
  } else if (type == 'data.frame') {
    tmp_string_array <- rep(NA, times = nrow(var) * (ncol(var) + 1))
    k <- 1
    for (i in seq(nrow(var))) {
      tmp_string_array[k] <- '\n'
      for (j in seq(ncol(var))) {
        tmp_string_array[k+j] <- var[i, j]
      }
      k <- k + ncol(var) + 1
    }
    output_string_tail <- paste0(tmp_string_array, collapse = ' ')
  } else {
    stop('Unkown type to output as ampl data file!')
  }

  # paste together
  paste0(output_string_head, output_string_tail, ';\n')
}

output_dat <- function(model, file) {
  if (missing(file)) {
    file <- prepare_tmp_file(type = 'dat')
  }
  # allow model-specific data output
  declarations <- get_ampl_data_output(model)
  output_string <- paste0(declarations, collapse = '')
  write(output_string, file = file)
  return(file)
}

output_mod <- function(model, file) {
  if (missing(file)) {
    file <- prepare_tmp_file(type = 'mod')
  }
  # allow model-specific model output
  declarations <- get_ampl_model_output(model)
  output <- paste0(declarations, collapse = '')
  write(output, file = file)
  return(file)
}
