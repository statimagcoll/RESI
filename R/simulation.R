#' Functions for simulations
#'
#' Function genearting longitudinal data
#' based on linear mixed-effects model.
#'
#' Function for the correlation structure for the errors within a subject
#' @param rho.e the correlation between two errors
#' @param time the values of time variables for a subject
#' @export
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
#' @param true_sigma.e: the true covariance matrix of errors within each subject
#' @param work.sigma.e the working covariance matrix specified
#' @param COV_beta the asymptotic covariance for $\sqrt{N} \times (\hat{\beta} - beta0)$
#' @param COV_beta_ind the asymptotic covariance for $\sqrt{N} \times (\hat{\beta} - beta0)$ under independence assumption
#' @param rho.e: correlation coef of the errors within a subject
#' @param fixed.trt: whether the trt assignment is fixed by design (TRUE) or random. use togther with pi
#' @param fixed.time whether the time points are fixed (TRUE, defuault) or random. If random, each time point (except baseline) has 60% chance to be observed
#' @export
sim_data_cont = function(N, S, pi, ni_range, true_sigma.e, work_sigma.e, COV_beta, COV_beta_ind, fixed.trt = TRUE, fixed.time = TRUE){

  # 1. GENERATING VARIABLE VALUES
  # number of measurements
  if (length(ni_range) == 1 | ni_range[1] == ni_range[2]) ni = rep(ni_range[1], times = N)
  if (length(ni_range) == 2 & ni_range[1] != ni_range[2]) ni = sample(ni_range[1]:ni_range[2], N, replace = TRUE)

  # subject ID
  id = rep(1:N, times = ni)

  # genearte the observation times (assuming everyone has the baseline measurement.)
  if (fixed.time){
    time_list = sapply(ni, function(x) c(0, sort(sample(1:(max(ni_range)-1), x-1, replace=FALSE))))
    time = c(unlist(time_list))
  } else {
    time_list = sapply(ni, function(x) runif(x,min = 0, max = 4) %>% sort())
    time = c(unlist(time_list))
  }

  # generate the trt group
  # the trt group for each subject
  if (fixed.trt) {
    trt = sample(rep(0:1, times = c(ceiling(N*pi), N - ceiling(N*pi)) ), replace = FALSE)
  } else {trt = rbinom(N, 1, pi)}
  # duplicate it for longitudinal data
  trt = rep(trt, times = ni)

  # Cov of random effects
  # G_mat = matrix(c(sigma0^2, rho.G*sigma.t*sigma0,
  #                  rho.G*sigma.t*sigma0, sigma.t^2), 2, 2)
  #
  # # generate random effects
  # gamma = mvtnorm::rmvnorm(N, c(0, 0), sigma = G_mat)
  # gamma_0 = rep(gamma[, 1], times = ni)
  # gamma_t = rep(gamma[, 2], times = ni)

  # measurement errors
  t = 0:(max(ni_range)-1)
  # if (length(ni_range) == 1 | ni_range[1] == ni_range[2]) { # if balanced design
  #   e =  sigma.e * mvtnorm::rmvnorm(N, rep(0, length(t)), sigma = ARMtx(t, rho.e)) %>% t() %>% c()
  # } else { # if unbalanced
  #   e =  sigma.e * c(unlist(sapply(time_list, function(x) mvtnorm::rmvnorm(1, rep(0, length(x)), sigma = ARMtx(x, rho.e))) ))
  # }
  e = mvtnorm::rmvnorm(N, rep(0, length(t)), sigma = true_sigma.e) %>% t %>% c # convert it to vector (by rows)


  # 2. CONVERTING RESI to BETA

  # The covariance matrix of Yi
  Sigma_y = true_sigma.e
  # Sigma_y = matrix(NA, nrow = ni_range[2], ncol = ni_range[2])
  # Design matrix for random effects
  # Note: this is the variance for \sqrt(N)(\hat{\beta}_long - \beta_0)
  var_int = COV_beta[1, 1]
  var_time = COV_beta[2, 2]
  var_trt = COV_beta[3, 3]

  # true SD's
  sd_int = sqrt(var_int)
  sd_time = sqrt(var_time)
  sd_trt = sqrt(var_trt)
  sd = c(sd_int, sd_time, sd_trt)

  # the beta's for a given RESI
  beta_int = sqrt(S[1]^2 * var_int)
  beta_time = sqrt(S[2]^2 * var_time)
  beta_trt = sqrt(S[3]^2 * var_trt)
  beta = c(beta_int, beta_time, beta_trt) # YEAH!!!!

  # 2.2 Deriving the true values of pm-RESI
  # The covariance of \hat{\beta} given X
  # The covariance matrix of \hat{\beta} under independence assumption

  tot_obs = sum(ni)
  var_int_ind = COV_beta_ind[1, 1]
  var_time_ind = COV_beta_ind[2, 2]
  var_trt_ind = COV_beta_ind[3, 3]

  # The true ESS = tot_obs * w
  ESS_int = tot_obs * var_int_ind / var_int
  ESS_time = tot_obs * var_time_ind / var_time
  ESS_trt = tot_obs * var_trt_ind / var_trt
  ESS = c(ESS_int, ESS_time, ESS_trt)

  # The true pm-RESI
  pm_resi_int = sqrt( S[1]^2 * N / ESS_int)
  pm_resi_time = sqrt( S[2]^2 * N / ESS_time)
  pm_resi_trt = sqrt( S[3]^2 * N / ESS_trt)
  pm_resi = c(pm_resi_int, pm_resi_time, pm_resi_trt)


  # 3. CALCULATING THE OUTCOME VALUES
  ## design matrix for fixed effects
  X = cbind(1, time, trt)

  # y = X %*% beta + gamma_0 + gamma_t * time + e
  y = X %*% beta  + e


  # 4. RETURNING SIMULATED DATA
  # data_sim = data.frame(id = id, num_obs = rep(ni, times = ni), time = time, trt = trt,
  #                       gamma_0 = gamma_0, gamma_t = gamma_t, error = e, y = y)
  data_sim = data.frame(id = id, num_obs = rep(ni, times = ni), time = time, trt = trt,
                        error = e, y = y)

  # # NEW: 5. making observed time points random
  #   if (fixed.time == FALSE){
  #     data_sim$observed = rbinom(nrow(data_sim), 1, prob = 0.75)
  #     data_sim$observed[data_sim$time == 0] = 1
  #     data_sim = subset(data_sim, observed == 1)
  #   }

  # return(list(data = data_sim, G = G_mat, N = N, ni = ni, true_beta = beta, true_sd = sd, cov_y = Sigma_y,
  #             pm_resi = pm_resi, ESS = ESS, info = "Function updated on 8/29/2022 3:16pm"))
  return(list(data = data_sim, N = N, ni = ni, true_beta = beta, true_sd = sd, cov_y = Sigma_y,
              pm_resi = pm_resi, ESS = ESS, info = "Function updated on 11/7/2022 2:38pm"))
}


#' @export
AR1_str = function(rho.e, num_visit){
  time = 0:(num_visit - 1)
  mtx <- matrix(NA, length(time), length(time))
  for (j in 1:length(time)){
    for (k in 1:length(time)){
      mtx[j, k] <- rho.e^(abs(time[j] - time[k]))
    }
  }
  return(mtx)
}

#' @export
comp_str = function(diag = 1, off_diag = 0.5, num_visit){
  mtx = matrix(NA, ncol = num_visit, nrow = num_visit)
  diag(mtx) = diag
  mtx[upper.tri(mtx)|lower.tri(mtx)] = off_diag
  return(mtx)
}

#' @export
AR1_to_comp = function(AR1){
  m = nrow(AR1)
  FUN = function(rho) sum(AR1[upper.tri(AR1)]) - m*(m-1)/2 * rho
  rho = uniroot(FUN, c(0, 10))$root
  comp = comp_str(1, rho, num_visit = m)
  return(comp)
}

#' @export
comp_to_AR1 = function(comp){
  m = nrow(comp)
  comp_rho = unique(comp[lower.tri(comp)])
  eq = paste0((m-1):1, "*rho", "^", 1:(m-1))
  FUN = function(rho) m*(m-1)/2 * comp_rho -  eval(parse(text = paste0(eq, collapse = " + ")))
  ar1_rho = uniroot(FUN, c(0, 10))$root
  ar1 = AR1_str(rho.e = ar1_rho, num_visit = m)
  return(ar1)
}
