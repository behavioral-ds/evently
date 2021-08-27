# implementation of CIKM'20: https://github.com/qykong/dual-mixture-hawkes-processes


get_param_names.hawkes_DMM <- function(model) {
  c('par') # placeholder
}

############# utility functions ##############
random_init_probabilities <- function(k) {
  values <- runif(k, 1, 100)
  vsum <- sum(values)
  values/vsum
}

log_sum_exp <- function(x) {
  m <- max(x)
  if (is.infinite(m) && m < 0) return(-Inf)
  m + log(sum(exp(x - m)))
}

########## mixturePLkernel helper functions ######################
# defines the log-likelihood function
get_ampl_likelihood.hawkes_mixturePLkernel <- function(model) {
  paste(
    'sum {cn in 1..HL} ( sum {cli in 1..CL} (  PKH[cli, cn] * (1 - 0^(L[cn]-J0[cn]-1)) *  (sum {i in J0[cn]+1..L[cn]-1} (',
    'log(sum {j in 1..i-1} (theta[cli] * c[cli]^theta[cli] * (time[cn,i] - time[cn,j] + c[cli]) ^ (-1-theta[cli])) + 1e-100)',
    '  ))',
    '));'
  )
}

# customized AMPL data output
get_ampl_data_output.hawkes_mixturePLkernel <- function(model) {
  output <- NextMethod()
  c(output,
    ampl_output_from_r('PKH', model$p_k_h_index, 'data.frame'))
}

# customized AMPL model output
get_ampl_model_output.hawkes_mixturePLkernel <- function(model) {
  paste('param HL > 0; param CL := ', model$cluster_number, '; ',
        'param ML > 0;',
        'param L {1..HL} >= 0;',
        'param magnitude {1..HL,1..ML} >= 0;',
        'param time {1..HL,1..ML} >= 0;',
        'param ind {1..HL,1..ML} >= 0;',
        'param J0 {1..HL} >= 0;',
        'param PKH {1..CL,1..HL} >=0;',
        'var c {1..CL} >= 0;',
        'var theta {1..CL} >= 0;',
        if (is.finite(model$init_par[[1]])) paste(glue::glue('let {get_param_names(model)} := {model$init_par[get_param_names(model)]};'), collapse = "\n") else '',
        'maximize Likelihood:',
        get_ampl_likelihood(model),
        # parameter constraints
        paste(glue::glue('subject to c_limit_{seq(model$cluster_number)}: c[{seq(model$cluster_number)}] <= 200;'), collapse = '\n'),
        paste(glue::glue('subject to theta_limit_{seq(model$cluster_number)}: theta[{seq(model$cluster_number)}] <= 5;'), collapse = '\n'),
        sep = '\n'
  )
}

# defines the parameters (`model$cluster_number` sets of theta and c)
get_param_names.hawkes_mixturePLkernel <- function(model) {
  c(sprintf('theta[%s]', seq(model$cluster_number)),
    sprintf('c[%s]', seq(model$cluster_number)))
}

# provide random initialisations for c and theta
generate_random_points.hawkes_mixturePLkernel <- function(model) {
  params <- get_param_names(model)
  inits <- lapply(seq_along(params), function(x) runif(10, min = .Machine$double.eps, max = 3))
  inits <- as.data.frame(do.call(cbind, inits))
  names(inits) <- params

  inits[9, ] <- NA
  inits[10, ] <- Inf
  inits
}

######### kernel mixture model fitting ##################
KMMEM <- function(data, k, max_iter = 10, ipopt_max_iter = 1000) {
  ps <- random_init_probabilities(k)

  params <- generate_random_points(new_hawkes(model_type='mixturePLkernel', model_vars = list(cluster_number = k)))[1,]
  ns <- names(params)

  kernel_log_likelihood <- function(cluster_no, param, hist) {
    c <- param[[sprintf('c[%s]', cluster_no)]]
    theta <- param[[sprintf('theta[%s]', cluster_no)]]

    if (nrow(hist) == 1) stop()
    sum(sapply(seq(2, nrow(hist)), function(.x) log(sum(theta * c^theta * (hist$time[.x] - hist$time[seq(.x-1)] + c) ^ (-1-theta)) +.Machine$double.xmin )))
  }

  iter <- 1
  total_likelihood <- -.Machine$double.xmax
  print('start em on kernel functions....')

  repeat {
    qs <- lapply(seq_along(data), function(l) sapply(seq(k), function(.x) (log(ps[.x]) + kernel_log_likelihood(.x, params, data[[l]]))))
    new_total_likelihood <- sum(sapply(seq_along(data), function(l) log_sum_exp(qs[[l]]) ))
    print(sprintf('iteration %s: %s', iter, round(new_total_likelihood, digits = 2)))
    if (abs(new_total_likelihood - total_likelihood) < 1e-2 || iter >= max_iter) {
      break
    } else {
      iter <- iter + 1
      total_likelihood <- new_total_likelihood
    }
    p_k_l <- lapply(seq(k), function(.x) sapply(seq_along(data), function(l) if (is.infinite(qs[[l]][.x])) 0 else 1/(sum(exp(qs[[l]][-.x] - qs[[l]][.x])) + 1) ))
    p_k_h_index <- lapply(seq(k), function(.x) {
      lapply(seq_along(data), function(l) list(k = .x, h = l, p = p_k_l[[.x]][[l]]))
    })
    p_k_h_index <- do.call(function(...) rbind.data.frame(..., make.row.names = FALSE),
                           p_k_h_index)

    ps <- sapply(seq(k), function(.x) sum(p_k_l[[.x]])/length(data))

    params <- as.data.frame(as.list(params))
    colnames(params) <- ns

    res <- fit_series(data, model_type = 'mixturePLkernel', init_pars = params,
                      observation_time = Inf, model_vars = list(cluster_number = k,
                                                                p_k_h_index = p_k_h_index),
                      cores = 1, ipopt_max_iter = ipopt_max_iter)
    params <- res$par
  }
  params <- lapply(seq(k), function(.x) c(theta = params[[sprintf('theta[%s]', .x)]],
                           c = params[[sprintf('c[%s]', .x)]]))
  return(list(probability = ps, params = params, p_k_h_index = p_k_h_index, total_likelihood = total_likelihood))
}

KMMEM_repeat <- function(..., times = 10, cores = 1) {
  models <- mclapply(seq(times), function(i) {
    tryCatch({
      KMMEM(...)
    }, error = function(e) {
      print(e)
      return(list(total_likelihood = -.Machine$double.xmax))
    })
  }, mc.cores = cores)
  models[[which.max(sapply(models, function(x) x[['total_likelihood']]))]]
}


######### Borel mixture model fitting ##################

BMMEM <- function(sizes, k = 1, max_iter = 50) {
  counts <- as.data.frame(table(sizes), stringsAsFactors = F)
  names(counts) <- c('size', 'count')
  counts$size <- as.integer(counts$size)

  ps <- random_init_probabilities(k)
  n_stars <- runif(k, 0, 1)

  borel <- function(x, lam) {
    stats::dpois(x-1, x*lam)/x
  }
  iter <- 1
  total_likelihood <- -.Machine$double.xmax
  repeat {
    qs <- lapply(seq(k), function(.x) sapply(seq_along(counts$size), function(l) ps[.x] * borel(counts$size[l], n_stars[.x])))
    new_total_likelihood <- sum(sapply(seq_along(counts$size), function(l) log(sum(sapply(seq(k), function(.x) qs[[.x]][l])))) * counts$count)
    if (new_total_likelihood - total_likelihood < 1e-2 || iter >= max_iter) {
      break
    } else {
      iter <- iter + 1
      total_likelihood <- new_total_likelihood
    }

    s_q <- sapply(seq_along(counts$size), function(l) sum(sapply(seq(k), function(.x) qs[[.x]][l])))
    p_k_l <- lapply(seq(k), function(.x) sapply(seq_along(counts$size), function(l) qs[[.x]][l]/s_q[l] ))
    n_stars <- sapply(seq(k), function(.x) sum((counts$size-1)*p_k_l[[.x]] * counts$count)/sum((counts$size)*p_k_l[[.x]] * counts$count) )

    ps <- sapply(seq(k), function(.x) sum(p_k_l[[.x]] * counts$count)/sum(counts$count) )
  }

  return(list(n_star = n_stars,
              p = ps,
              total_likelihood = total_likelihood))
}

BMMEM_repeat <- function(..., times = 10, cores = 1) {
  models <- mclapply(seq(times), function(i) {BMMEM(...) }, mc.cores = cores)
  models[[which.max(sapply(models, function(x) x[['total_likelihood']]))]]
}

######## fit the two mixture models together #########
fit_series_by_model.hawkes_DMM <- function(model, cores, init_pars,
                                           parallel_type, .init_no,
                                           ipopt_max_iter = 1000, ...) {
  hists <- model$data
  clusters <- model$cluster_no
  times <- model$times
  if (is.null(times)) times <- 10

  sizes <- sapply(hists, nrow)
  cat('start BMM\n')
  if (!is.null(clusters)) {
    BMM_clusters <- clusters[1]
    n_star_p <- BMMEM_repeat(sizes, k = BMM_clusters, cores = cores, times = times)
    if (length(clusters) == 2) kernel_clusters <- clusters[2] else kernel_clusters <- clusters
  } else {
    # test k from 2 to 10 and choose the best by AIC
    trials <- lapply(seq(2, 10), function(k) {
      n_star_p <- BMMEM_repeat(sizes, k = k, cores = cores, times = times)
      list(n_star_p = n_star_p,
           aic_n_star = -n_star_p$total_likelihood * 2 + k * 2,
           k = k)
    })
    best_ind <- which.min(sapply(trials, function(x) x[['aic_n_star']]))
    BMM_clusters <- trials[[best_ind]]$k
    kernel_clusters <- trials[[best_ind]]$k
    n_star_p <- trials[[best_ind]]$n_star_p
    cat(sprintf('best number of clusters is %s\n', BMM_clusters))
  }

  cat('BMM done; start clustering kernel functions\n')

  # remove single event cascades as they won't be computed in KMM anyway
  keeped_hists <- hists[sapply(hists, function(h) nrow(h) >= 2)]
  if (!is.null(model$max_event_length) && model$max_event_length > 0) {
    cat(sprintf('Capping number of events in KMM to %s', model$max_event_length))
    keeped_hists <- lapply(keeped_hists, function(hist) hist[seq(min(model$max_event_length, nrow(hist))), ])
  }

  # if no cascades left then return here
  if (length(keeped_hists) == 0) {
    model$par <- list(n_star = n_star_p$n_star,
                      p_n_star = n_star_p$p,
                      params = NA,
                      p_params = NA,
                      kernel_clusters = kernel_clusters,
                      BMM_clusters = BMM_clusters)
    return(model)
  }
  cat('start KMM\n')
  kernel_clusters <- min(length(keeped_hists), kernel_clusters)

  res <- KMMEM_repeat(keeped_hists, k = kernel_clusters, times = times,
                      cores = cores, ipopt_max_iter=ipopt_max_iter)
  cat('done KMM\n')
  params <- res
  model$par <- list(n_star = n_star_p$n_star,
                    n_star_probability = n_star_p$p,
                    kernel_params = res$params,
                    kernel_params_probability = res$probability,
                    kernel_clusters = kernel_clusters,
                    BMM_clusters = BMM_clusters)
  model$value <- res$value
  model$upper_bound <- NULL
  model$lower_bound <- NULL
  model$init_par <- NULL
  model
}
