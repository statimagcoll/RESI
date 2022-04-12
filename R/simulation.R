#' Functions for simulations
#'
#' Function genearting longitudinal data
#' based on linear mixed-effects model.
#'
#' Function for the correlation structure for the errors within a subject
#' @param rho.e the correlation between two errors
#' @param time the values of time variables for a subject
ARMtx <- function(time, rho.e){
  mtx <- matrix(NA, length(time), length(time))
  for (j in 1:length(time)){
    for (k in 1:length(time)){
      mtx[j, k] <- rho.e^(abs(time[j] - time[k]))
    }
  }
  return(mtx)
}

#' Function for simulating continuous outcomes
#' @param N total sample size (i.e., num of subjects)
#' @param S = the values of RESI for the fixed effects
#' @param pi = proportion or probability of being assigned to trt group
#' @param ni_range = range of num of measurements
#' @param rho.G = correlation coef between random intercepts and slopes
#' @param sigma0: SD of random intercepts
#' @param sigma.t: SD of random slopes
#' @param sigma.e: SD of errors
#' @param rho.e: correlation coef of the errors within a subject
#' @param fixed.design: whether the trt assignment is fixed by design (TRUE) or random. use togther with pi
#' @export
sim_data_cont = function(N, S, pi, ni_range, rho.G, sigma0, sigma.t, sigma.e, rho.e, fixed.design = TRUE){

# 1. GENERATING VARIABLE VALUES
  # number of measurements
  if (length(ni_range) == 1 | ni_range[1] == ni_range[2]) ni = rep(ni_range[1], times = N)
  if (length(ni_range) == 2 & ni_range[1] != ni_range[2]) ni = sample(ni_range[1]:ni_range[2], N, replace = TRUE)

  # subject ID
  id = rep(1:N, times = ni)

  # genearte the observation times (assuming everyone has the baseline measurement.)
  time_list = sapply(ni, function(x) c(0, sort(sample(1:(max(ni_range)-1), x-1, replace=FALSE))))
  time = c(unlist(time_list))
  # c(unlist(sapply(ni-1, seq, from = 0)))  # minus one since we start at time = 0

  # generate the trt group
  if (fixed.design) {
    trt = sample(rep(0:1, times = c(ceiling(N*pi), N - ceiling(N*pi)) ), replace = FALSE)
  } else {trt = rbinom(N, 1, pi)}
  trt = rep(trt, times = ni)

  # Cov of random effects
  G_mat = matrix(c(sigma0^2, rho.G*sigma.t*sigma0,
                   rho.G*sigma.t*sigma0, sigma.t^2), 2, 2)

  # generate random effects
  gamma = mvtnorm::rmvnorm(N, c(0, 0), sigma = G_mat)
  gamma_0 = rep(gamma[, 1], times = ni)
  gamma_t = rep(gamma[, 2], times = ni)

  # measurement errors
  if (length(ni_range) == 1 | ni_range[1] == ni_range[2]) { # if balanced design
    t = 0:(max(ni_range)-1)
    e =  sigma.e * mvtnorm::rmvnorm(N, rep(0, length(t)), sigma = ARMtx(t, rho.e))
  } else { # if unbalanced
    e =  sigma.e * c(unlist(sapply(time_list, function(x) mvtnorm::rmvnorm(1, rep(0, length(x)), sigma = ARMtx(x, rho.e))) ))
  }

# 2. CONVERTING RESI TO BETA
  ## The covariance matrix of Yi
  Sigma_y = matrix(NA, nrow = ni_range[2], ncol = ni_range[2])
  # Design matrix for random effects
  time = (0:(ni_range[2]-1))
  Z = cbind(1, time)
  # Cov(Y_i) = Z^T G Z + R_i
  Sigma_y = Z %*% G_mat %*% t(Z) + ARMtx(time = time, rho.e = rho.e)

  # The true cov of \hat{beta} given X
  ## 1. for variable 'trt'


  # the variance of \hat{\beta}_{trt} when trt = 1
  X_1 = cbind(1, time = time, trt = 1) # design matrix
  info_beta_1 = t(X_1) %*% Sigma_y %*% X_1 * N # the whole info matrix
  var_int_1 = solve(info_beta_1[1, 1]) # for the intercept
  var_time_1 = solve(info_beta_1[2, 2]) # for time variable
  var_trt_1 = solve(info_beta_1[3, 3]) # for trt

  # the variance of \hat{\beta}_{trt} when trt = 0
  X_0 = cbind(1, time = time, trt = 0) # design matrix
  info_beta_0 = t(X_0) %*% Sigma_y %*% X_0 * N # the whole info matrix
  var_int_0 = solve(info_beta_0[1, 1]) # for the intercept
  var_time_0 = solve(info_beta_0[2, 2]) # for the time variable
  # for \beta_trt, the var is just ZERO?
  var_trt_0 = 0

  # unconditional variances
  var_int = var_int_1 * pi + var_int_0 * (1-pi)
  var_trt = var_trt_1 * pi
  var_time = var_time_1 * pi + var_time_0 * (1-pi)

  # the beta's for a given RESI
  beta_int = sqrt(N * S^2 * var_int)
  beta_trt = sqrt(N * S^2 * var_trt)
  beta_time = sqrt(N * S^2 * var_time)
  beta = c(beta_int, beta_time, beta_trt) # YEAH!!!!

# 3. CALCULATING THE OUTCOME VALUES
  ## design matrix for fixed effects
  X = cbind(1, time, trt)

  y = X %*% beta + gamma_0 + gamma_t * time + e

# 4. RETURNING SIMULATED DATA
  data_sim = data.frame(id = id, num_obs = rep(ni, times = ni), time = time, trt = trt,
                        gamma_0 = gamma_0, gamma_t = gamma_t, error = e, y = y )

  return(list(data = data_sim, G = G_mat, N = N, ni = ni))
}


