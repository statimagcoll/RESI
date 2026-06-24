# Asymptotic confidence intervals for RESI
# Implements the normal approximation and quadratic-form (Imhof) CIs described in
# Zhang et al. (2025) "On the Asymptotic Distribution of the Robust Effect Size Index"
# and the quadratic-form notes by Vandekar.
#
# Exported:
#   resi_pe_asymptotic()   -- point estimates + asymptotic CIs (both methods)
#
# Internal pipeline:
#   .resi_precompute()     -- model-level sandwich ingredients (once per model)
#   .resi_contrast()       -- per-contrast whitened vector + Sigma_R
#   .resi_ci_normal_signed()  -- signed normal CI  (coefficients, m1=1)
#   .resi_ci_normal_unsigned() -- unsigned truncated normal CI (anova)
#   .resi_ci_qf()          -- quadratic-form Imhof CI (anova)
#   .get_L_coef()          -- contrast matrix for a single coefficient
#   .get_L_anova()         -- contrast matrix for an anova term (Type 2)


# ============================================================
#  Model-level precomputation
# ============================================================

#' @noRd
.resi_precompute <- function(model, type = "HC3") {

  type <- match.arg(type, c("HC3", "const", "HC", "HC0", "HC1", "HC2",
                             "HC4", "HC4m", "HC5"))
  if (type == "HC") type <- "HC0"

  # ---- basic model objects ----
  X <- model.matrix(model)
  if (any(alias <- is.na(coef(model)))) X <- X[, !alias, drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)
  e <- residuals(model, "response")

  # ---- HC leverage weights (applied to mean-parameter rows of psi only) ----
  h <- hatvalues(model)
  one_m_h <- pmax(1 - h, .Machine$double.eps)
  p_int <- max(1L, as.integer(round(sum(h))))
  sqrtw <- switch(toupper(type),
    "HC0"  = rep(1, n),
    "HC1"  = rep(sqrt(n / max(n - p, 1L)), n),
    "HC2"  = 1 / sqrt(one_m_h),
    "HC3"  = 1 / one_m_h,
    "HC4"  = {
      nhp <- (n / p_int) * h
      one_m_h^(-pmin(4, nhp) / 2)
    },
    "HC4M" = {
      nhp <- (n / p_int) * h
      delta <- pmin(1, nhp) + pmin(1.5, nhp)
      one_m_h^(-delta / 2)
    },
    "HC5"  = {
      nhp <- (n / p_int) * h
      k <- 0.7
      deltaCap <- max(4, (n / p_int) * k * max(h))
      delta <- pmin(nhp, deltaCap)
      one_m_h^(-delta / 4)
    },
    rep(1, n)   # fallback (const handled via is_const below)
  )
  is_const <- (type == "const")

  # ============================================================
  #  lm branch: theta = (phi, beta),  m = p+1
  # ============================================================
  if (inherits(model, "lm") && !inherits(model, "glm")) {

    phi <- summary(model)$sigma^2
    m   <- p + 1L

    # --- per-observation psi and psi' ---
    psi_list <- lapply(seq_len(n), function(i) {
      xi <- X[i, , drop = FALSE]        # 1 x p
      ei <- e[i]
      psi_phi  <- (ei^2 - phi) / (2 * phi^2)
      psi_beta <- (ei / phi) * xi
      cbind(psi_phi, psi_beta)           # 1 x m
    })

    psiprime_list <- lapply(seq_len(n), function(i) {
      xi <- X[i, , drop = FALSE]
      ei <- e[i]
      d_phi_phi  <- (phi - 2 * ei^2) / (2 * phi^3)
      d_phi_beta <- -ei * xi / phi^2
      d_beta_phi <- t(d_phi_beta)
      d_beta_beta <- -crossprod(xi) / phi
      rbind(c(d_phi_phi, d_phi_beta),
            cbind(d_beta_phi, d_beta_beta))
    })

    # --- A and B matrices ---
    A_full <- -.resi_mean_list(psiprime_list)
    A_full <- .resi_sym(A_full)
    A_inv  <- .resi_safe_inv(A_full)

    # B: HC weights applied to mean-parameter (beta) rows only
    B_full <- if (is_const) {
      .resi_mean_outer(psi_list, rep(1, n), lm_phi_row = TRUE)
    } else {
      .resi_mean_outer(psi_list, sqrtw, lm_phi_row = TRUE)
    }
    B_full <- .resi_sym(B_full)

    cov_theta <- .resi_sym(A_inv %*% B_full %*% A_inv)

    # --- analytic dA/dtheta and dB/dtheta ---
    XtX_n    <- crossprod(X) / n
    XtX_n_hc <- if (is_const) XtX_n else crossprod(X * sqrtw) / n

    # dA/d_phi only; all dA/d_beta_j = 0 for homoscedastic lm (V=1)
    dA_phi <- rbind(
      cbind(-1 / phi^3,     matrix(0, 1, p)),
      cbind(matrix(0, p, 1), -XtX_n / phi^2)
    )
    dA_dtheta <- matrix(0, m * m, m)
    dA_dtheta[, 1] <- as.vector(dA_phi)

    # dB/d_phi only
    e4mean <- mean(e^4)
    e3mean <- mean(e^3)
    dB_phi <- rbind(
      cbind((phi^2 - 2 * e4mean) / phi^5,
            -(3 / (2 * phi^4)) * t(colMeans(X) * e3mean)),
      cbind(-(3 / (2 * phi^4)) * colMeans(X) * e3mean,
            -XtX_n_hc / phi^2)
    )
    dB_dtheta <- matrix(0, m * m, m)
    dB_dtheta[, 1] <- as.vector(dB_phi)

    theta_hat <- c(phi, coef(model))
    lm_model  <- TRUE

  # ============================================================
  #  GLM branch: theta = beta,  m = p
  # ============================================================
  } else if (inherits(model, "glm")) {

    m <- p
    lm_model <- FALSE
    phi <- NULL

    ef    <- sandwich::estfun(model)
    ef    <- ef[, colnames(X), drop = FALSE]
    w_vec <- weights(model, type = "working")
    mu_hat <- pmin(pmax(fitted(model), .Machine$double.eps),
                   1 - .Machine$double.eps)

    psi_list      <- lapply(seq_len(n), function(i)
      matrix(ef[i, ], nrow = 1))
    psiprime_list <- lapply(seq_len(n), function(i)
      -w_vec[i] * crossprod(X[i, , drop = FALSE]))

    A_full <- -.resi_mean_list(psiprime_list)
    A_full <- .resi_sym(A_full)
    A_inv  <- .resi_safe_inv(A_full)

    B_full <- if (is_const) {
      .resi_mean_outer(psi_list, rep(1, n), lm_phi_row = FALSE)
    } else {
      .resi_mean_outer(psi_list, sqrtw, lm_phi_row = FALSE)
    }
    B_full <- .resi_sym(B_full)

    cov_theta <- .resi_sym(A_inv %*% B_full %*% A_inv)

    # analytic dA/dtheta (= dB/dtheta under canonical link)
    dA_dtheta <- matrix(0, m * m, m)
    for (j in seq_len(m)) {
      vj <- w_vec * (1 - 2 * mu_hat) * X[, j]
      dA_dtheta[, j] <- as.vector(crossprod(X, vj * X) / n)
    }
    dB_dtheta <- dA_dtheta

    theta_hat <- coef(model)

  } else {
    stop(".resi_precompute: model class not supported (lm or glm only)")
  }

  list(
    model      = model,
    X          = X,
    n          = n,
    p          = p,
    m          = m,
    lm_model   = lm_model,
    phi        = phi,
    theta_hat  = theta_hat,
    A_inv      = A_inv,
    B_full     = B_full,
    cov_theta  = cov_theta,
    dA_dtheta  = dA_dtheta,
    dB_dtheta  = dB_dtheta,
    is_const   = is_const,
    type       = type
  )
}


# ============================================================
#  Per-contrast computation
# ============================================================

#' @noRd
#' @param precomp output of .resi_precompute()
#' @param L_model contrast matrix in *beta* (coefficient) space (m1 x p)
#' @param vcovmat_n Optional matrix n * vcovfunc(model) in beta space (p x p).
#'   When supplied, used for Sigma_beta so that Stilde targets the same population
#'   quantity as the vcovfunc-based point estimate. cov_theta (from precomp) is
#'   still used for Sigma_R via the delta method.
#' @return list with R_beta, Stilde, dR_dtheta (m1 x m), Sigma_R (m1 x m1),
#'         beta_hat, Sigma_beta, m1
.resi_contrast <- function(precomp, L_model, vcovmat_n = NULL) {

  m1     <- nrow(L_model)
  m      <- precomp$m
  n      <- precomp$n
  A_inv  <- precomp$A_inv
  Sig    <- precomp$cov_theta     # m x m  (used for Sigma_R)
  dA     <- precomp$dA_dtheta     # m^2 x m
  dB     <- precomp$dB_dtheta     # m^2 x m

  # Extend L to theta space: prepend column of zeros for phi if lm
  L <- if (precomp$lm_model) cbind(0, L_model) else L_model  # m1 x m

  # ---- point estimate ----
  beta_hat <- drop(L %*% precomp$theta_hat)           # m1

  # Sigma_beta: use vcovmat_n (vcovfunc-based, in beta space) when supplied so
  # that Stilde is consistent with the vcovfunc-based RESI point estimate.
  # Note: cov_theta[beta,beta] ≈ n * vcovHC, so for the default HC3 vcovfunc
  # this is numerically transparent; it matters when vcovfunc != type.
  Sigma_beta <- if (!is.null(vcovmat_n)) {
    .resi_sym(L_model %*% vcovmat_n %*% t(L_model))
  } else {
    .resi_sym(L %*% Sig %*% t(L))
  }

  # ---- EVD of Sigma_beta ----
  eig <- eigen(Sigma_beta, symmetric = TRUE)
  d   <- pmax(eig$values, .Machine$double.eps)         # eigenvalues (m1)
  V   <- eig$vectors                                   # m1 x m1

  sqrtd <- sqrt(d)
  Phat  <- V %*% diag(1 / sqrtd, m1) %*% t(V)         # Sigma_beta^{-1/2}

  R_beta <- Phat %*% beta_hat                          # m1-vector, whitened RESI
  Stilde <- sqrt(sum(R_beta^2))                        # ||R_beta||

  # ---- Lyapunov weight matrix ----
  # W[j,k] = 1 / (sqrt(d[j]) * sqrt(d[k]) * (sqrt(d[j]) + sqrt(d[k])))
  W <- outer(sqrtd, sqrtd, function(a, b) 1 / (a * b * (a + b)))  # m1 x m1

  VT_beta <- drop(t(V) %*% beta_hat)   # m1-vector (V-basis projection of beta_hat)

  # ---- dR_dtheta: m1 x m matrix ----
  direct <- Phat %*% L                   # (a) direct term: m1 x m

  Achain <- matrix(0, m1, m)
  Bchain <- matrix(0, m1, m)
  for (k in seq_len(m)) {
    # (b) A-chain
    dA_k     <- matrix(dA[, k], m, m)
    dSig_A   <- -(A_inv %*% dA_k %*% Sig + Sig %*% dA_k %*% A_inv)
    dSigB_A  <- L %*% dSig_A %*% t(L)
    M_A      <- t(V) %*% dSigB_A %*% V  # m1 x m1 in V-basis
    Achain[, k] <- -V %*% (W * M_A) %*% VT_beta

    # (c) B-chain (zero for parametric / const)
    if (!precomp$is_const) {
      dB_k     <- matrix(dB[, k], m, m)
      dSig_B   <- A_inv %*% dB_k %*% A_inv
      dSigB_B  <- L %*% dSig_B %*% t(L)
      M_B      <- t(V) %*% dSigB_B %*% V
      Bchain[, k] <- -V %*% (W * M_B) %*% VT_beta
    }
  }

  dR_dtheta <- direct + Achain + Bchain   # m1 x m

  # ---- Sigma_R = Cov(sqrt(n) R_hat) ----
  Sigma_R <- dR_dtheta %*% Sig %*% t(dR_dtheta)   # m1 x m1
  Sigma_R <- .resi_sym(Sigma_R)

  list(
    m1         = m1,
    n          = n,
    beta_hat   = beta_hat,
    Sigma_beta = Sigma_beta,
    R_beta     = R_beta,
    Stilde     = Stilde,
    dR_dtheta  = dR_dtheta,
    Sigma_R    = Sigma_R
  )
}


# ============================================================
#  CI functions
# ============================================================

#' Signed normal CI for a single coefficient (m1 = 1)
#' @noRd
.resi_ci_normal_signed <- function(contrast, alpha = 0.05) {
  Sigma_R <- contrast$Sigma_R  # scalar (1x1)
  n       <- contrast$n
  Stilde  <- contrast$Stilde   # |Z|/sqrt(n) sign preserved via R_beta sign
  R_beta  <- contrast$R_beta   # signed (scalar)

  # signed: S_pm = R_beta (= Z/sqrt(n) when m1=1, Sigma_beta=1)
  S_signed <- R_beta   # scalar
  se       <- sqrt(as.numeric(Sigma_R) / n)
  z        <- qnorm(1 - alpha / 2)
  c(LCI = S_signed - z * se, UCI = S_signed + z * se)
}

#' Unsigned truncated normal CI (anova, m1 >= 1)
#' @noRd
.resi_ci_normal_unsigned <- function(contrast, alpha = 0.05) {
  Sigma_R <- contrast$Sigma_R   # m1 x m1
  R_beta  <- contrast$R_beta    # m1
  Stilde  <- contrast$Stilde
  n       <- contrast$n
  m1      <- contrast$m1

  # df-corrected center: Shat = sqrt(max(0, Stilde^2 - m1/n))
  # Aligns CI center with chisq2S point estimator; same correction as QF CI.
  Shat <- sqrt(max(0, Stilde^2 - m1 / n))

  # Null upper bound: (1-alpha) quantile of Stilde under H0.
  # Under H0, sqrt(n)*||R_hat|| ~ ||N(0, I_m1)|| = chi(m1) because Sigma_R -> I
  # asymptotically. So P(Stilde >= UCI | S=0) = alpha gives
  # UCI = sqrt(qchisq(1-alpha, df=m1) / n).
  # For m1=1: sqrt(qchisq(1-alpha, 1)/n) = qnorm(1-alpha/2)/sqrt(n).
  UCI_null <- sqrt(qchisq(1 - alpha, df = m1) / n)

  if (Shat <= 0) {
    return(c(LCI = 0, UCI = UCI_null))
  }

  u       <- R_beta / Stilde                              # unit direction
  sigma2S <- as.numeric(t(u) %*% Sigma_R %*% u)         # scalar variance
  se      <- sqrt(pmax(sigma2S, 0) / n)

  # truncate_ci handles boundary; centered at Shat
  bounds <- .resi_truncate_ci(Shat, se, n, m1, alpha)

  # When LCI clips to 0, the directional delta-method UCI can underestimate
  # spread for m1>1 (single-direction 1D normal vs chi(m1) null distribution).
  # Use the chi(m1) null upper as a floor.
  if (bounds[1] == 0) {
    bounds[2] <- max(bounds[2], UCI_null)
  }

  c(LCI = bounds[1], UCI = bounds[2])
}

#' Quadratic-form Imhof CI (anova, m1 >= 1)
#' @noRd
.resi_ci_qf <- function(contrast, alpha = 0.05) {
  Sigma_R <- contrast$Sigma_R
  R_beta  <- contrast$R_beta
  Stilde  <- contrast$Stilde
  n       <- contrast$n
  m1      <- contrast$m1

  # T2_obs  <- max(0, n * Stilde^2 - m1)    # df-corrected test statistic (Shat^2 * n)
  # Shat    <- sqrt(T2_obs / n)              # df-corrected point estimate (for search bounds only)
  T2_obs <- n * Stilde^2                  # raw test statistic
  Shat   <- sqrt(max(0, Stilde^2 - m1/n)) # (for search bounds only)

  # EVD of Sigma_R
  eigR       <- eigen(Sigma_R, symmetric = TRUE)
  lambda     <- pmax(eigR$values, .Machine$double.eps)
  U          <- eigR$vectors

  # Estimated non-centrality direction (unit length in U-basis)
  delta_unit <- if (Stilde > 0) drop(t(U) %*% R_beta) / Stilde else
    rep(0, length(lambda))

  # Normal-approx SE of Stilde (for search interval sizing).
  # Delta-method: se_S = sqrt(u^T Sigma_R u / n), u = R_beta/Stilde.
  # Under H0, n*Stilde^2 ~ chi^2_{m1}, so Var(Stilde) -> m1/(2n) * (2/m1) = 1/n
  # giving SD(Stilde) ~ sqrt(m1/n) / sqrt(2) ... more precisely the SD of
  # chi(m1)/sqrt(n) is sqrt(m1/n) (since Var(chi(m1)) ~ m1 for large m1).
  # Use sqrt(m1/n) as a floor so the bracket is never degenerate near null.
  u_dir        <- if (Stilde > 0) R_beta / Stilde else eigR$vectors[, 1]
  sigma2S      <- as.numeric(t(u_dir) %*% Sigma_R %*% u_dir)
  se_S_delta   <- sqrt(pmax(sigma2S / n, .Machine$double.eps))
  se_S_floor   <- sqrt(m1 / n)
  se_S         <- max(se_S_delta, se_S_floor)

  # nc_k = n * S_b^2 * delta_unit_k^2 / lambda_k  (non-centrality for chi^2_1 k-th term)
  nc_fun <- function(S_b) n * S_b^2 * delta_unit^2 / lambda

  # Probability function P(Q(S_b) >= T2_obs)
  # For m1=1: use pchisq (no imhof precision issues)
  # For m1>1: use imhof
  if (m1 == 1L) {
    q <- T2_obs / lambda[1]
    # For large q or ncp, use the exact identity chi^2(1,lambda) = (Z + sqrt(lambda))^2:
    #   P(chi^2(1,lambda) >= q) = pnorm(sqrt(lambda) - sqrt(q)) + pnorm(-sqrt(lambda) - sqrt(q))
    # This avoids pnchisq convergence failures (which arise when q and ncp are both large,
    # e.g. ~1e12, as happens in ill-conditioned GLMs at small n).  The formula is
    # algebraically exact and pnorm is numerically stable at arbitrarily large arguments.
    prob_fn <- function(S_b) {
      nc <- nc_fun(S_b)
      if (!is.finite(nc)) return(1.0)
      if (q > 1e4 || nc > 1e4) {
        sq  <- sqrt(max(q,  0))
        snc <- sqrt(max(nc, 0))
        return(pnorm(snc - sq, lower.tail = FALSE) + pnorm(-snc - sq))
      }
      pchisq(q, df = 1, ncp = nc, lower.tail = FALSE)
    }
  } else {
    prob_fn <- function(S_b) {
      nc <- nc_fun(S_b)
      # Cap individual non-centrality components to avoid imhof overflow
      nc <- pmin(nc, 1e10 * max(T2_obs, 1))
      suppressWarnings(CompQuadForm::imhof(T2_obs, lambda = lambda, delta = nc)$Qq)
    }
  }

  # ---- Lower bound ----
  # P(Q(0)) = P(central weighted chi-sq >= T2_obs)
  p0 <- prob_fn(0)
  if (is.na(p0) || p0 < 0) p0 <- 0   # imhof can return negative for extreme cases

  if (p0 >= alpha / 2) {
    LCI <- 0
  } else {
    lower_candidates <- Shat * c(0.9, 0.7, 0.5, 0.3, 0.1)
    bracket_lo <- NA_real_
    for (lo in lower_candidates) {
      plo <- tryCatch(prob_fn(lo), error = function(e) NA_real_)
      if (!is.na(plo) && is.finite(plo) && plo >= 0 && plo < alpha / 2) {
        bracket_lo <- lo
        break
      }
    }
    if (is.na(bracket_lo)) bracket_lo <- 0

    LCI <- tryCatch(
      uniroot(function(s) prob_fn(s) - alpha / 2,
              lower = bracket_lo, upper = Shat,
              tol = se_S * 1e-3, extendInt = "upX")$root,
      error = function(e) {
        warning("QF lower CI bound search failed; using 0.")
        0
      }
    )
    LCI <- max(LCI, 0)
  }

  # ---- Upper bound ----
  # When p0 = P(Q(0) >= T2_obs) >= 1-alpha/2, T2_obs falls in the lower alpha/2
  # tail of the null distribution. Strict test inversion gives an empty upper set
  # because prob_fn is increasing and never drops to 1-alpha/2. Fallback: use the
  # (1-alpha) quantile of Stilde under H0: n*Stilde^2 ~ chi^2_{m1}, so
  # UCI_null = sqrt(qchisq(1-alpha, df=m1) / n). Same formula as the normal CI.
  # The same fallback is used when the bracket search fails.
  # Null upper bound: (1-alpha) quantile of Stilde under H0.
  # Under H0, n*Stilde^2 ~ chi^2_{m1} asymptotically, so
  # UCI_null = sqrt(qchisq(1-alpha, df=m1) / n).
  UCI_null <- sqrt(qchisq(1 - alpha, df = m1) / n)

  if (p0 >= 1 - alpha / 2) {
    # T2_obs is in the lower alpha/2 tail of the null distribution;
    # strict test inversion gives an empty upper set. Use the null upper bound.
    UCI <- UCI_null
  } else {
    # Normal case: uniroot in (Shat, bracket_up) where prob_fn crosses 1-alpha/2
    upper_search <- Shat + se_S * c(3, 6, 12, 25, 50, 100) * qnorm(1 - alpha / 2)
    bracket_up   <- NA_real_
    for (up in upper_search) {
      pup <- tryCatch(prob_fn(up), error = function(e) NA_real_)
      if (is.na(pup) || !is.finite(pup) || pup >= 1 - alpha / 2) {
        bracket_up <- up; break
      }
    }
    if (is.na(bracket_up)) {
      UCI <- UCI_null
    } else {
      UCI <- tryCatch(
        uniroot(function(s) {
          p <- prob_fn(s)
          if (is.na(p) || !is.finite(p)) return(1 - (1 - alpha / 2))
          p - (1 - alpha / 2)
        },
                lower = Shat, upper = bracket_up,
                tol = se_S * 1e-3)$root,
        error = function(e) UCI_null
      )
    }
  }

  c(LCI = max(LCI, 0), UCI = UCI)
}


# ============================================================
#  Contrast-matrix helpers
# ============================================================

#' @noRd
.get_L_coef <- function(model, coef_name) {
  coefs <- names(coef(model))
  coefs <- coefs[!is.na(coef(model))]
  idx   <- which(coefs == coef_name)
  if (length(idx) == 0) stop("Coefficient not found: ", coef_name)
  L <- matrix(0, 1, length(coefs))
  L[1, idx] <- 1
  L
}

#' Type-2 Anova contrast matrix for a single term (in beta space)
#' Adapted from get_L_anova2 in functions.R
#' @noRd
.get_L_anova_term <- function(model, term, vcovmat) {
  coefs       <- names(coef(model))
  not_aliased <- !is.na(coef(model))
  names_terms <- labels(terms(model))
  which_term  <- which(term == names_terms)
  factors     <- attr(terms(model), "factors")
  asgn        <- attr(model.matrix(model), "assign")
  asgn[!not_aliased] <- NA

  subs_term <- which(asgn == which_term)

  # relatives: terms that contain 'term' as a lower-order effect
  relatives <- setdiff(seq_along(names_terms), which_term)
  relatives <- relatives[sapply(names_terms[relatives], function(t2)
    all(factors[, term] <= factors[, t2]))]
  subs_rel  <- unlist(lapply(relatives, function(r) which(asgn == r)))

  Ip <- diag(sum(not_aliased))
  hyp1 <- Ip[subs_rel, , drop = FALSE]
  hyp2 <- Ip[c(subs_rel, subs_term), , drop = FALSE]

  if (nrow(hyp1) == 0) {
    L_out <- hyp2
  } else {
    # Orthogonalize: columns of hyp2 orthogonal to span(hyp1) under vcovmat metric
    L_out <- t(.resi_conjcomp(t(hyp1), t(hyp2), vcovmat))
  }
  L_out <- L_out[!apply(L_out, 1, function(x) all(x == 0)), , drop = FALSE]
  L_out
}

#' @noRd
.resi_conjcomp <- function(X, Z, ip = diag(nrow(X))) {
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0) return(Z)
  Z %*% qr.Q(xq, complete = TRUE)[, -(seq_len(xq$rank)), drop = FALSE]
}


# ============================================================
#  Main exported function
# ============================================================

#' Robust Effect Size Index with Asymptotic Confidence Intervals
#'
#' Computes RESI point estimates and asymptotic confidence intervals using
#' either a normal approximation (Zhang et al., 2025) or a quadratic-form
#' (Imhof/Davies) approach.
#'
#' @param model.full Fitted \code{lm} or \code{glm} model object.
#' @param data Data frame of model data.
#' @param vcovfunc Variance estimator for RESI point estimates. Default:
#'   \code{sandwich::vcovHC}.
#' @param coefficients Logical; include coefficient table. Default \code{TRUE}.
#' @param anova Logical; include anova table. Default \code{TRUE}.
#' @param alpha Numeric; significance level. Default \code{0.05}.
#' @param ci.method Character; \code{"normal"} (truncated normal) or
#'   \code{"qf"} (quadratic-form Imhof). Default \code{"normal"}.
#' @param type Character; HC type for the M-estimator sandwich used in CI
#'   construction. Default \code{"HC3"}.
#' @param unbiased Logical; use bias-corrected RESI point estimate. Default
#'   \code{TRUE}.
#' @param Anova.args List; additional arguments passed to \code{car::Anova}.
#' @param vcov.args List; additional arguments passed to \code{vcovfunc}.
#' @param ... Ignored.
#' @return A list of class \code{"resi"} with \code{coefficients} and/or
#'   \code{anova} tables containing RESI point estimates and CIs.
#' @importFrom sandwich vcovHC estfun
#' @importFrom CompQuadForm imhof
#' @importFrom car Anova
#' @importFrom lmtest coeftest
#' @importFrom stats coef hatvalues residuals fitted weights qnorm qchisq pchisq
#'   model.matrix terms uniroot
#' @export
resi_pe_asymptotic <- function(model.full,
                                data,
                                vcovfunc  = sandwich::vcovHC,
                                coefficients = TRUE,
                                anova        = TRUE,
                                alpha        = 0.05,
                                ci.method    = c("normal", "qf"),
                                type         = "HC3",
                                unbiased     = TRUE,
                                Anova.args   = list(),
                                vcov.args    = list(),
                                ...) {

  ci.method <- match.arg(ci.method)
  model     <- model.full
  n         <- nrow(model.matrix(model))

  if (missing(data)) {
    data <- if (!is.null(model$model)) model$model else
      stop("data argument required")
  }

  # ---- vcovfunc with extra args ----
  if (length(vcov.args) > 0) {
    vcovfunc2 <- function(x) {
      args <- c(list(x), vcov.args)
      do.call(vcovfunc, args)
    }
  } else {
    vcovfunc2 <- vcovfunc
  }

  # ---- model-level precomputation ----
  precomp <- .resi_precompute(model, type = type)

  # ---- vcov matrix for point estimates (same as resi_pe) ----
  vcovmat <- tryCatch(vcovfunc2(model), error = function(e)
    stop("vcovfunc failed: ", conditionMessage(e)))

  output <- list()

  # ============================================================
  #  Coefficients table  (signed normal CI always)
  # ============================================================
  if (coefficients) {
    coef_names <- names(coef(model))[!is.na(coef(model))]
    is_lm      <- inherits(model, "lm") && !inherits(model, "glm")
    rdf        <- if (is_lm) model$df.residual else NULL

    # point estimates via coeftest
    ctab <- tryCatch(
      lmtest::coeftest(model, vcov. = vcovmat),
      error = function(e) stop("coeftest failed")
    )
    torZ <- if ("z value" %in% colnames(ctab)) "z" else "t"

    coef_df <- data.frame(
      Estimate    = ctab[, "Estimate"],
      `Std. Error` = ctab[, "Std. Error"],
      Statistic   = ctab[, paste(torZ, "value")],
      `p-value`   = ctab[, paste0("Pr(>|", torZ, "|)")],
      check.names = FALSE,
      row.names   = rownames(ctab)
    )
    colnames(coef_df)[3:4] <- c(paste(torZ, "value"),
                                  paste0("Pr(>|", torZ, "|)"))
    if (torZ == "z") {
      coef_df$RESI <- suppressWarnings(z2S(coef_df[[3]], n, unbiased))
    } else {
      coef_df$RESI <- suppressWarnings(t2S(coef_df[[3]], rdf, n, unbiased))
    }

    # asymptotic signed CIs
    ci_mat <- do.call(rbind, lapply(coef_names, function(cn) {
      L_mod  <- .get_L_coef(model, cn)
      contr  <- tryCatch(.resi_contrast(precomp, L_mod, vcovmat_n = n * vcovmat),
                         error = function(e) NULL)
      if (is.null(contr)) return(c(LCI = NA_real_, UCI = NA_real_))
      .resi_ci_normal_signed(contr, alpha = alpha)
    }))

    coef_df[, paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")] <- ci_mat
    output$coefficients <- coef_df
  }

  # ============================================================
  #  Anova table  (unsigned CI; method depends on ci.method)
  # ============================================================
  if (anova) {
    is_lm <- inherits(model, "lm") && !inherits(model, "glm")

    # point estimates via car::Anova (matches resi_pe behaviour)
    if (is_lm) {
      anova_tab <- tryCatch(
        suppressMessages(do.call(car::Anova,
          c(list(mod = model, vcov. = vcovmat), Anova.args))),
        error = function(e) stop("car::Anova failed")
      )
      anova_tab <- anova_tab[rownames(anova_tab) != "Residuals", , drop = FALSE]
      rdf        <- model$df.residual
      anova_tab$RESI <- f2S(anova_tab[, "F"], anova_tab[, "Df"], rdf, n)
    } else {
      anova_tab <- tryCatch(
        suppressMessages(do.call(car::Anova,
          c(list(mod = model, test.statistic = "Wald",
                 vcov. = vcovmat), Anova.args))),
        error = function(e) stop("car::Anova failed")
      )
      anova_tab <- anova_tab[rownames(anova_tab) != "Residuals", , drop = FALSE]
      anova_tab$RESI <- chisq2S(anova_tab[, "Chisq"], anova_tab[, "Df"], n)
    }

    term_names <- rownames(anova_tab)

    # CI for each anova term
    ci_mat <- do.call(rbind, lapply(term_names, function(term) {
      L_mod <- tryCatch(
        .get_L_anova_term(model, term, vcovmat),
        error = function(e) NULL
      )
      if (is.null(L_mod) || nrow(L_mod) == 0)
        return(c(LCI = NA_real_, UCI = NA_real_))

      contr <- tryCatch(.resi_contrast(precomp, L_mod, vcovmat_n = n * vcovmat),
                        error = function(e) NULL)
      if (is.null(contr)) return(c(LCI = NA_real_, UCI = NA_real_))

      if (ci.method == "qf") {
        tryCatch(.resi_ci_qf(contr, alpha = alpha),
                 error = function(e) c(LCI = NA_real_, UCI = NA_real_))
      } else {
        .resi_ci_normal_unsigned(contr, alpha = alpha)
      }
    }))

    anova_tab[, paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")] <- ci_mat
    class(anova_tab) <- c("anova_resi", class(anova_tab))
    output$anova <- anova_tab
  }

  output$alpha     <- alpha
  output$ci.method <- ci.method
  output$type      <- type
  class(output)    <- c("resi", "list")
  output
}


# ============================================================
#  Utility helpers (all unexported)
# ============================================================

#' @noRd
.resi_mean_list <- function(lst) {
  Reduce("+", lst) / length(lst)
}

#' @noRd
.resi_mean_outer <- function(psi_list, sqrtw, lm_phi_row = FALSE) {
  n <- length(psi_list)
  m <- ncol(psi_list[[1]])
  B <- matrix(0, m, m)
  for (i in seq_len(n)) {
    psi_i <- psi_list[[i]]   # 1 x m
    if (lm_phi_row) {
      # phi row unweighted, beta rows weighted
      psi_w          <- psi_i
      psi_w[1, 2:m]  <- psi_i[1, 2:m] * sqrtw[i]
    } else {
      psi_w <- psi_i * sqrtw[i]
    }
    B <- B + crossprod(psi_w)
  }
  B / n
}

#' @noRd
.resi_sym <- function(M) (M + t(M)) / 2

#' @noRd
.resi_safe_inv <- function(M) {
  tryCatch(chol2inv(chol(M)), error = function(e) solve(M))
}

#' Truncated CI for unsigned RESI (Algorithm 1, Zhang et al. 2025)
#' @noRd
.resi_truncate_ci <- function(Shat, se, n, m1, alpha = 0.05) {
  z1 <- qnorm(1 - alpha / 2)
  SL <- Shat - z1 * se
  SU <- Shat + z1 * se

  if (SL > 0) return(c(SL, SU))

  # probability mass at/below zero under H0: P(T^2 > n*Shat^2 | S=0)
  # Under H0, T^2 ~ chi^2_{m1}, so P(T^2 > m1 + n*Shat^2) (since Shat^2 = (T^2-m1)/n)
  gamma_p <- pchisq(m1 + n * Shat^2, df = m1, lower.tail = FALSE)

  if (gamma_p < alpha / 2) {
    SU <- Shat + qnorm(1 - (alpha - gamma_p)) * se
  } else {
    SU <- Shat + qnorm(1 - alpha) * se
  }
  c(0, SU)
}
