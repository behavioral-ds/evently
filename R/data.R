
output_dat.data.frame <- function(data, file) {
  file.create(file)
  output_string <- paste("param L := ", nrow(data), ";\nparam: magnitude time :=\n", sep = "")
  for (i in 1:nrow(data)) {
    output_string <- paste(output_string, "     ", i,"    ", data$magnitude[i], "  ", data$time[i], sep = "")
    if (i != nrow(data)) {
      output_string <- paste(output_string, "\n", sep = "")
    }
  }
  output_string <- paste(output_string, ";", sep = "")
  write(output_string, file = file)
}

output_dat.list <- function(data, file) {
  file.create(file)
  lengthes <- sapply(data, nrow)

  output_string <- paste0("param HL:= ", length(data), ";\n param L := ",
                          paste(seq_along(lengthes), lengthes, collapse = ' '),
                          ";\nparam ML:= ", max(lengthes),
                          ";\nparam: magnitude time :=")

  tmp_string_array <- rep(NA, times = sum(lengthes) * 5)
  k <- 1
  for (i in seq_along(data)) {
    for (j in seq(lengthes[i])) {
      tmp_string_array[k] <- '\n'
      tmp_string_array[k + 1] <- as.character(i)
      tmp_string_array[k + 2] <- as.character(j)

      tmp_string_array[k + 3] <- as.character(data[[i]]$magnitude[j])
      tmp_string_array[k + 4] <- as.character(data[[i]]$time[j])

      k <- k + 5
    }
  }
  tmp_string_array <- paste0(tmp_string_array, collapse = ' ')
  output_string <- paste0(output_string, tmp_string_array, ";\n param J0 := ")

  for (i in seq_along(data)) {
    output_string <- paste(output_string, i, min(which(data[[i]]$time > 0)) - 1, sep = " ")
  }
  output_string <- paste(output_string, ";", sep = "")

  write(output_string, file = file)
}

get_mod_param_def.data.frame <- function(data) {
  paste(
    '# define data',
    'param L > 1;',
    'param magnitude {1..L} >= 0;',
    'param time {1..L} >= 0;',
    'param J0 := min {i in 1..L : time[i] > 0} i - 1;',
    '',
    sep = '\n'
  )
}

get_mod_param_def.list <- function(data) {
  paste(
    '# define data',
    'param HL > 0;',
    'param ML > 0;',
    'param L {1..HL} >= 0;',
    'param magnitude {1..HL,1..ML} >= 0;',
    'param time {1..HL,1..ML} >= 0;',
    'param J0 {1..HL} >= 0;',
    '',
    sep = '\n'
  )
}

output_dat.default <- function(data, ...) {
  stop('unknown training data format!')
}

# function dispatcher
output_dat <- function(data, ...) {
  UseMethod('output_dat')
}

get_mod_param_def <- function(data) {
  UseMethod('get_mod_param_def')
}
