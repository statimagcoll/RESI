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

#' Function for simulating longitudinal continuous outcomes
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
  # the trt group for each subject
  if (fixed.design) {
    trt = sample(rep(0:1, times = c(ceiling(N*pi), N - ceiling(N*pi)) ), replace = FALSE)
  } else {trt = rbinom(N, 1, pi)}
  # duplicate it for longitudinal data
  trt = rep(trt, times = ni)

  # Cov of random effects
  G_mat = matrix(c(sigma0^2, rho.G*sigma.t*sigma0,
                   rho.G*sigma.t*sigma0, sigma.t^2), 2, 2)

  # generate random effects
  gamma = mvtnorm::rmvnorm(N, c(0, 0), sigma = G_mat)
  gamma_0 = rep(gamma[, 1], times = ni)
  gamma_t = rep(gamma[, 2], times = ni)

  # measurement errors
  t = 0:(max(ni_range)-1)
  if (length(ni_range) == 1 | ni_range[1] == ni_range[2]) { # if balanced design
    e =  sigma.e * mvtnorm::rmvnorm(N, rep(0, length(t)), sigma = ARMtx(t, rho.e)) %>% t() %>% c()
  } else { # if unbalanced
    e =  sigma.e * c(unlist(sapply(time_list, function(x) mvtnorm::rmvnorm(1, rep(0, length(x)), sigma = ARMtx(x, rho.e))) ))
  }

# 2. CONVERTING RESI TO BETA
  ## The covariance matrix of Yi
  Sigma_y = matrix(NA, nrow = ni_range[2], ncol = ni_range[2])
  # Design matrix for random effects
  time_points = (0:(ni_range[2]-1))
  Z = cbind(1, time_points)
  # Cov(Y_i) = Z^T G Z + R_i
  Sigma_y = Z %*% G_mat %*% t(Z) + ARMtx(time = time_points, rho.e = rho.e)

  # The true cov of \hat{beta} given X
  # Using numeric methods to find the variance of each parameter estimator
  # (assuming the model is correctly specified)
  # The covariance matrix of \hat{beta}
  sum = 0
  rep = 1e5
  for (i in 1:rep){
    # X_i^T Sigma^{-1} T
    X = cbind(1, time = time_points, trt = rbinom(1, 1, pi))
    sum = sum + t(X) %*% solve(Sigma_y) %*% X
  }
  COV_beta = solve(sum) * rep # = N * (the asymptotic cov of beta hat)

  var_int = COV_beta[1, 1] / N
  var_time = COV_beta[2, 2] / N
  var_trt = COV_beta[3, 3] / N

  # true SD's
  sd_int = sqrt(var_int)
  sd_time = sqrt(var_time)
  sd_trt = sqrt(var_trt)
  sd = c(sd_int, sd_time, sd_trt)

  # the beta's for a given RESI
  beta_int = sqrt(N * S[1]^2 * var_int)
  beta_time = sqrt(N * S[2]^2 * var_time)
  beta_trt = sqrt(N * S[3]^2 * var_trt)
  beta = c(beta_int, beta_time, beta_trt) # YEAH!!!!

# 2.2 Deriving the true values of pm-RESI
  # The covariance of \hat{\beta} given X
  # assuming the model is correctly specified.
  # The covariance matrix of \sqrt{N}*(\hat{\beta} - \beta_0) under independence assumption

  # rep = 1e3
  # if (fixed.design) {
  #   trt_temp = sample(rep(0:1, times = c(ceiling(rep*pi), rep - ceiling(rep*pi)) ), replace = FALSE)
  # } else {trt_temp = rbinom(rep, 1, pi)}
  # # duplicate it for longitudinal data
  # trt_temp = rep(trt_temp, each = unique(ni))
  # X = cbind(1, time = time_points, trt = trt_temp)
  # sum = t(X) %*% solve(diag(rep(2, rep*unique(ni)))) %*% X

  sum = 0
  rep = 1e3
  for (i in 1:rep){
    # generate the trt group
    # the trt group for each subject
    if (fixed.design) {
      trt_temp = sample(rep(0:1, times = c(ceiling(N*pi), N - ceiling(N*pi)) ), replace = FALSE)
    } else {trt_temp = rbinom(N, 1, pi)}
    # duplicate it for longitudinal data
    trt_temp = rep(trt_temp, times = ni)
    X = cbind(1, time = time, trt = trt_temp)
    sum = sum + t(X) %*% X * (sigma.e^2 + sigma0^2)
  }

  COV_beta_ind = solve(sum) * rep  # = the variance of \hat{\beta|ind}

  tot_obs = sum(ni)
  var_int_ind = COV_beta_ind[1, 1]
  var_time_ind = COV_beta_ind[2, 2]
  var_trt_ind = COV_beta_ind[3, 3]
  # The true pm-RESI
  pm_resi_int = sqrt( S[1]^2 * var_int / var_int_ind)
  pm_resi_time = sqrt( S[2]^2 * var_time / var_time_ind)
  pm_resi_trt = sqrt( S[3]^2 * var_trt / var_trt_ind)
  pm_resi = c(pm_resi_int, pm_resi_time, pm_resi_trt)
  # The true ESS
  ESS_int = sum(ni) * var_int_ind / var_int
  ESS_time = sum(ni) * var_time_ind / var_time
  ESS_trt = sum(ni) * var_trt_ind / var_trt
  ESS = c(ESS_int, ESS_time, ESS_trt)

# 3. CALCULATING THE OUTCOME VALUES
  ## design matrix for fixed effects
  X = cbind(1, time, trt)

  y = X %*% beta + gamma_0 + gamma_t * time + e

# 4. RETURNING SIMULATED DATA
  data_sim = data.frame(id = id, num_obs = rep(ni, times = ni), time = time, trt = trt,
                        gamma_0 = gamma_0, gamma_t = gamma_t, error = e, y = y)

  return(list(data = data_sim, G = G_mat, N = N, ni = ni, true_beta = beta, true_sd = sd, cov_y = Sigma_y,
              pm_resi = pm_resi, ESS = ESS, info = "Function updated on 8/29/2022"))
}



#' Function for the simulations with logistic regression models
#' @param n total sample size (i.e., num of subjects)
#' @param S = the values of RESI for the grouping parameter x: logit(P) = alpha + x beta
#' @param p = the probability of the outcome being 1
#' @param pi = proportion or probability of being assigned to trt group
#' @param fixed.design: whether the trt assignment is fixed by design (TRUE) or random. use togther with pi
#' @export

sim_data_cs_binary <- function(n, S, r, alpha, m, p, pi, fixed.design, num.cores){
  # 1. DATA GENERATION
  # simulate x and y
  if (fixed.design) {
    sequence = rep(0:1, times = m*c(n - ceiling(n*pi), ceiling(n*pi)))
    x = sample(sequence, replace = FALSE) %>% matrix(nrow = n, ncol = m)
  } else {
    x = rbinom(n*m, 1, pi) %>% matrix(nrow = n, ncol = m)
  }




  #
}


