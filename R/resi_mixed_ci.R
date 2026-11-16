# Asymptotic RESI confidence intervals for lmerMod and nlme::lme fits with random effects
# nested in an outermost cluster (derivation: AIgams notes/resi_mixed_ci.qmd).
#
# Var(R_hat) is estimated from school-level (outermost-cluster) influence
# contributions: (a) outcome noise through beta_hat, (b) sampling of the
# per-cluster information I_s (random design), (c) sampling of the meat
# U_s U_s' (robust covariance only), plus, for the model-based covariance,
# (psi) the variance-component uncertainty. The resulting contrast lists are
# passed to the CI functions in resi_asymptotic.R.

# Marginal-covariance quantities for theta (lme4 relative covariance
# parameters) and sigma. V = sigma^2 (I + Z Lambda Lambda' Z'), handled by
# Woodbury with M = Lambda' Z' Z Lambda + I.
.mixedCiGls <- function(base, theta, sigma) {
  Lt <- base$Lambdat
  Lt@x <- theta[base$Lind]
  LZt <- Lt %*% base$Zt
  ch <- Matrix::Cholesky(Matrix::tcrossprod(LZt), perm = TRUE, LDL = FALSE, Imult = 1)
  vinv <- function(B) {
    as.matrix(B - Matrix::crossprod(LZt, Matrix::solve(ch, LZt %*% B, system = "A")))
  }
  VX <- vinv(base$X)
  XVX <- crossprod(base$X, VX)
  beta <- drop(solve(XVX, crossprod(VX, base$y)))
  r <- drop(base$y - base$X %*% beta)
  list(VX = VX / sigma^2, info = XVX / sigma^2, beta = beta, r = r,
       rVr = sum(r * drop(vinv(r))),
       logdetM = 2 * as.numeric(Matrix::determinant(ch, sqrt = TRUE)$modulus),
       logdetXVX = as.numeric(determinant(XVX)$modulus))
}

.mixedCiDeviance <- function(base, psi) {
  k <- length(psi)
  gls <- .mixedCiGls(base, psi[-k], psi[k])
  N <- length(base$y); p <- ncol(base$X); s2 <- psi[k]^2
  if (base$REML) {
    (N - p) * log(2 * pi * s2) + gls$logdetM + gls$logdetXVX + gls$rVr / s2
  } else {
    N * log(2 * pi * s2) + gls$logdetM + gls$rVr / s2
  }
}

.mixedCiBase <- function(fit, cluster) {
  base <- if (inherits(fit, "lme")) .mixedCiBaseLme(fit) else {
    list(X = lme4::getME(fit, "X"), y = lme4::getME(fit, "y"), Zt = lme4::getME(fit, "Zt"),
         Lambdat = lme4::getME(fit, "Lambdat"), Lind = lme4::getME(fit, "Lind"),
         theta = lme4::getME(fit, "theta"), sigma = stats::sigma(fit),
         REML = lme4::isREML(fit), beta = lme4::fixef(fit),
         flist = as.list(lme4::getME(fit, "flist")))
  }
  fl <- base$flist
  if (is.null(cluster)) cluster <- names(fl)[which.min(vapply(fl, nlevels, integer(1L)))]
  cl <- if (length(cluster) == 1L && cluster %in% names(fl)) fl[[cluster]] else cluster
  cl <- droplevels(as.factor(cl))
  zt <- methods::as(base$Zt, "TsparseMatrix")
  if (any(tapply(as.integer(cl)[zt@j + 1L], zt@i, function(x) length(unique(x))) > 1L)) {
    stop("Mixed-model RESI CIs require every random effect to be nested within `cluster`.")
  }
  base$cluster <- cl
  base
}

# Reasons an nlme::lme fit is outside the supported class (homoscedastic,
# independent within-group errors; unstructured, diagonal or identity
# random-effect covariances), or NULL when supported.
.mixedCiLmeUnsupported <- function(fit) {
  ms <- fit$modelStruct
  if (!is.null(ms$varStruct)) return("`weights` (variance functions) are not supported")
  if (!is.null(ms$corStruct)) return("`correlation` structures are not supported")
  ok <- vapply(ms$reStruct, function(pd) inherits(pd, c("pdLogChol", "pdSymm", "pdDiag", "pdIdent")),
               logical(1L))
  if (!all(ok)) return("random-effect covariance classes other than pdSymm/pdLogChol/pdDiag/pdIdent are not supported")
  NULL
}

# The lme4 representation of an lme fit: V = sigma^2 (I + Z Lambda Lambda' Z'),
# with Lambda block diagonal in the per-level relative Cholesky factors and
# theta their free entries (lower triangle for unstructured, diagonal for
# pdDiag, one shared value for pdIdent), so the deviance Hessian respects the
# model's covariance constraints.
.mixedCiBaseLme <- function(fit) {
  why <- .mixedCiLmeUnsupported(fit)
  if (!is.null(why)) stop("Mixed-model RESI CIs for lme: ", why, ".")
  design <- .resiMixedModelDesign(fit, "lme")
  data <- design$data
  X <- stats::model.matrix(stats::terms(stats::formula(fit)), data = data,
                           contrasts.arg = fit$contrasts)
  beta <- nlme::fixef(fit)
  if (!identical(colnames(X), names(beta))) stop("Could not rebuild the lme fixed-effect design.")
  y <- nlme::getResponse(fit)
  re <- fit$modelStruct$reStruct
  n <- nrow(X)
  zi <- zj <- integer(0); zx <- numeric(0)
  li <- lj <- lt <- integer(0)
  theta <- numeric(0)
  offset <- 0L
  # reStruct is stored innermost first; getGroups(level = k) and fit$groups
  # run outermost first, so levels are matched by name.
  for (k in seq_along(re)) {
    pd <- re[[k]]
    Zk <- stats::model.matrix(stats::formula(pd), data = data)
    q <- ncol(Zk)
    grp <- droplevels(as.factor(nlme::getGroups(fit, level = match(names(re)[k], names(fit$groups)))))
    J <- nlevels(grp)
    gi <- as.integer(grp)
    for (c in seq_len(q)) {
      zi <- c(zi, seq_len(n)); zj <- c(zj, offset + (gi - 1L) * q + c); zx <- c(zx, Zk[, c])
    }
    D <- as.matrix(pd)
    if (inherits(pd, c("pdDiag", "pdIdent"))) {
      Lk <- diag(sqrt(diag(D)), q)
      pos <- cbind(seq_len(q), seq_len(q))
      tid <- if (inherits(pd, "pdIdent")) rep(1L, q) else seq_len(q)
      th <- if (inherits(pd, "pdIdent")) Lk[1, 1] else diag(Lk)
    } else {
      Lk <- t(chol(D))
      pos <- which(lower.tri(Lk, diag = TRUE), arr.ind = TRUE)
      pos <- pos[order(pos[, 2], pos[, 1]), , drop = FALSE]
      tid <- seq_len(nrow(pos))
      th <- Lk[pos]
    }
    # Lambdat = t(Lambda): entry (r, c) of the lower factor sits at (c, r).
    for (j in seq_len(J)) {
      li <- c(li, offset + (j - 1L) * q + pos[, 2])
      lj <- c(lj, offset + (j - 1L) * q + pos[, 1])
      lt <- c(lt, length(theta) + tid)
    }
    theta <- c(theta, th)
    offset <- offset + J * q
  }
  Zt <- Matrix::t(Matrix::sparseMatrix(i = zi, j = zj, x = zx, dims = c(n, offset)))
  Lambdat <- Matrix::sparseMatrix(i = li, j = lj, x = seq_along(li), dims = c(offset, offset))
  Lind <- lt[Lambdat@x]
  Lambdat@x <- theta[Lind]
  list(X = X, y = as.numeric(y), Zt = Zt, Lambdat = Lambdat, Lind = Lind,
       theta = theta, sigma = fit$sigma, REML = identical(fit$method, "REML"),
       beta = beta, flist = rev(design$groups))
}

# Per-cluster information contributions I_s (g x p^2) and scores U_s (g x p).
.mixedCiClusterPieces <- function(base, gls) {
  X <- base$X; p <- ncol(X)
  prods <- do.call(cbind, lapply(seq_len(p), function(b) X * gls$VX[, b]))
  list(I_s = rowsum(prods, base$cluster),
       U_s = rowsum(gls$VX * gls$r, base$cluster))
}

# CR2 cluster-robust covariance (Pustejovsky & Tipton, 2018) with working
# covariance V, matching clubSandwich::vcovCR(fit, type = "CR2") for lmerMod and lme
# fits with random effects nested in `cluster`, without its n_s x n_s matrix
# roots. clubSandwich adjusts W_s e_s by A_s = V_s^{1/2} G_s^{-1/2} V_s^{1/2},
# G_s = V_s^{1/2} (V_s - X_s M X_s') V_s^{1/2}. With V_s = sigma^2 (I + Z_s Z_s')
# (Z_s = rows of Z Lambda) and P = V_s^{-1/2} X_s chol(M)', G_s = V_s (I - P P') V_s
# equals sigma^4 I off span(Z_s, P), so its inverse root is computed in that
# small invariant subspace.
.mixedCiCR2 <- function(fit, cluster = NULL) {
  base <- .mixedCiBase(fit, cluster)
  gls <- .mixedCiGls(base, base$theta, base$sigma)
  M <- solve(gls$info)
  Cm <- chol(M)
  s2 <- base$sigma^2
  LZ <- Matrix::t(base$Lambdat %*% base$Zt)
  rows <- split(seq_along(base$y), base$cluster)
  u <- vapply(rows, function(idx) {
    Zs <- LZ[idx, , drop = FALSE]
    Zs <- as.matrix(Zs[, Matrix::colSums(Zs != 0) > 0, drop = FALSE])
    # Thin left singular vectors of Zs via the small q_s x q_s eigenproblem
    # (LAPACK dgesdd occasionally fails to converge on these matrices).
    # (Zs has no columns when theta = 0, i.e. a singular fit.)
    if (ncol(Zs) == 0L) {
      sv <- list(d = numeric(0), u = matrix(0, length(idx), 0L))
    } else {
      ez <- eigen(crossprod(Zs), symmetric = TRUE)
      keep <- ez$values > max(ez$values) * 1e-12
      sv <- list(d = sqrt(ez$values[keep]))
      sv$u <- Zs %*% (ez$vectors[, keep, drop = FALSE] %*% diag(1 / sv$d, length(sv$d)))
    }
    vpow <- function(B, a) {
      base$sigma^(2 * a) * (B + sv$u %*% (((1 + sv$d^2)^a - 1) * crossprod(sv$u, B)))
    }
    VmhX <- vpow(base$X[idx, , drop = FALSE], -0.5)
    P <- VmhX %*% t(Cm)
    Tq <- qr.Q(qr(cbind(sv$u, P)))
    TU <- crossprod(Tq, sv$u)
    TVT <- s2 * (diag(ncol(Tq)) + TU %*% (sv$d^2 * t(TU)))
    TP <- crossprod(Tq, P)
    eg <- eigen(TVT %*% (diag(ncol(Tq)) - tcrossprod(TP)) %*% TVT, symmetric = TRUE)
    Gmh <- eg$vectors %*% (t(eg$vectors) / sqrt(pmax(eg$values, .Machine$double.eps)))
    h <- vpow(matrix(gls$r[idx]), 0.5)
    Th <- crossprod(Tq, h)
    drop(crossprod(VmhX, (h - Tq %*% Th) / s2^2 + Tq %*% (Gmh %*% Th)))
  }, numeric(ncol(base$X)))
  out <- M %*% tcrossprod(u) %*% M
  dimnames(out) <- list(colnames(base$X), colnames(base$X))
  out
}

# Directional derivative of R = Sigma^{-1/2} b with respect to Sigma in
# direction dSigma, via the Lyapunov solution for d(Sigma^{1/2}).
.mixedCiDR <- function(eig, R, dSigma) {
  sq <- sqrt(eig$values)
  V <- eig$vectors
  M <- crossprod(V, dSigma %*% V)
  dHalf <- V %*% (M / outer(sq, sq, `+`)) %*% t(V)
  # d(Sigma^{-1/2}) = -Sigma^{-1/2} d(Sigma^{1/2}) Sigma^{-1/2}
  ihalf <- V %*% (t(V) / sq)
  -drop(ihalf %*% dHalf %*% R)
}

.mixedCiContrast <- function(L, n_k, beta, Vhat, Itot, pieces, robust,
                             design, counts = NULL, Jpsi = NULL, Vpsi = NULL) {
  m1 <- nrow(L)
  Sigma <- n_k * L %*% Vhat %*% t(L)
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  ihalf <- eig$vectors %*% (t(eig$vectors) / sqrt(eig$values))
  R <- drop(ihalf %*% L %*% beta)
  p <- ncol(L); g <- nrow(pieces$U_s)
  Iinv <- solve(Itot)
  Ibar <- Itot / g
  Bbar <- crossprod(pieces$U_s) / g
  # n_k is itself a sum of per-cluster unit counts when it counts units below
  # the cluster level, so it moves with the cluster information.
  LVL <- L %*% Vhat %*% t(L)
  dcounts <- if (design && !is.null(counts)) counts - mean(counts) else rep(0, g)
  d <- t(vapply(seq_len(g), function(s) {
    U <- pieces$U_s[s, ]
    out <- if (robust) drop(ihalf %*% L %*% (Iinv %*% U)) else rep(0, m1)
    dV <- matrix(0, p, p)
    if (design) {
      dI <- matrix(pieces$I_s[s, ], p, p) - Ibar
      dV <- if (robust) -Iinv %*% dI %*% Vhat - Vhat %*% dI %*% Iinv else -Iinv %*% dI %*% Iinv
    }
    if (robust) dV <- dV + Iinv %*% (tcrossprod(U) - Bbar) %*% Iinv
    if (design || robust) {
      out <- out + .mixedCiDR(eig, R, n_k * L %*% dV %*% t(L) + dcounts[s] * LVL)
    }
    out
  }, numeric(m1)))
  if (m1 == 1L) d <- matrix(d, ncol = 1L)
  d <- d * sqrt(g / (g - 1))
  var_R <- crossprod(d)
  if (!robust) var_R <- var_R + diag(m1) / n_k
  if (!is.null(Jpsi)) var_R <- var_R + Jpsi %*% Vpsi %*% t(Jpsi)
  var_R <- (var_R + t(var_R)) / 2
  # Influence values for the signed Cornish-Fisher skewness (m1 = 1). For the
  # model-based covariance, the outcome-noise term uses the cluster scores.
  phi <- NULL
  if (m1 == 1L) {
    phi_d <- if (robust) drop(d) else
      drop(d) + sqrt(g / (g - 1)) * drop(pieces$U_s %*% (Iinv %*% t(ihalf %*% L)))
    phi <- g * phi_d
  }
  list(m1 = m1, n = n_k, R_beta = R, Stilde = sqrt(sum(R^2)),
       Sigma_R = n_k * var_R, var_R = var_R, phi_tilde = phi, g = g)
}

.mixedCiInterval <- function(contrast, ci.method, alpha, signed) {
  if (signed) {
    if (ci.method == "cf") {
      cf <- contrast
      cf$n <- contrast$g
      cf$Sigma_R <- contrast$g * contrast$var_R
      return(.resi_ci_cf_signed(cf, alpha))
    }
    return(.resi_ci_normal_signed(contrast, alpha))
  }
  switch(ci.method,
         normal = .resi_ci_normal_unsigned(contrast, alpha),
         qf = .resi_ci_qf(contrast, alpha),
         cf = .resi_ci_cf(contrast, alpha))
}

#' Mixed-model RESI confidence intervals (internal)
#'
#' Computes asymptotic CIs for the coefficient (signed) and ANOVA (unsigned)
#' RESI of an `lmerMod` or `lme` fit, for terms normalized by `n_coef`/`n_anova`.
#' @param fit `lmerMod` or `lme` fit whose random effects are nested in `cluster`.
#' @param Vhat Covariance of the fixed effects used for the point estimates.
#' @param robust Logical; `Vhat` is a cluster-robust sandwich.
#' @param n_coef,n_anova Named term-specific normalizing counts.
#' @param ci.method One of `"normal"`, `"qf"`, `"cf"` (coefficients use the
#'   signed normal CI for `"qf"`).
#' @param design Logical; include the design (cluster information) sampling
#'   term, i.e. treat the cluster covariate profiles as random.
#' @param cluster Name of the outermost grouping factor.
#' @param n_level_coef,n_level_anova Optional named grouping-factor names
#'   whose counts are `n_coef`/`n_anova` (the `n level` columns).
#' @return A list with, for each requested `ci.method`, `coefficients` and
#'   `anova` matrices of LCI and UCI, and `se`, the estimated SD of each
#'   one-degree-of-freedom RESI estimate.
#' @keywords internal
#' @noRd
.resi_ci_mixed <- function(fit, Vhat, robust, n_coef, n_anova = NULL,
                           ci.method = c("normal", "qf", "cf"), alpha = 0.05,
                           design = TRUE, cluster = NULL,
                           n_level_coef = NULL, n_level_anova = NULL) {
  ci.method <- match.arg(ci.method, several.ok = TRUE)
  base <- .mixedCiBase(fit, cluster)
  gls <- .mixedCiGls(base, base$theta, base$sigma)
  pieces <- .mixedCiClusterPieces(base, gls)
  beta <- base$beta
  Vhat <- as.matrix(Vhat)[names(beta), names(beta)]
  p <- length(beta)
  assign <- attr(base$X, "assign")
  terms <- .mixedFixedTermLabels(fit)
  fl <- base$flist
  unit_counts <- function(n_k, lev) {
    if (n_k == nlevels(base$cluster)) return(NULL)
    if (is.null(lev) || !(lev %in% names(fl))) {
      lev <- names(fl)[vapply(fl, nlevels, integer(1L)) == n_k]
      if (length(lev) != 1L) return(NULL)
    }
    as.numeric(tapply(as.character(fl[[lev]]), base$cluster,
                      function(x) length(unique(x))))
  }
  contrasts <- c(
    lapply(seq_len(p), function(j) {
      nm <- names(beta)[j]
      list(L = diag(p)[j, , drop = FALSE], n = n_coef[[nm]],
           counts = unit_counts(n_coef[[nm]], n_level_coef[nm]))
    }),
    if (!is.null(n_anova)) lapply(names(n_anova), function(t) {
      list(L = diag(p)[assign == match(t, terms), , drop = FALSE], n = n_anova[[t]],
           counts = unit_counts(n_anova[[t]], n_level_anova[t]))
    }))

  Jpsi_all <- NULL
  if (!robust) {
    psi_hat <- c(base$theta, base$sigma)
    Rfun <- function(psi) {
      k <- length(psi)
      V <- solve(.mixedCiGls(base, psi[-k], psi[k])$info)
      unlist(lapply(contrasts, function(cn) {
        S <- cn$n * cn$L %*% V %*% t(cn$L)
        e <- eigen((S + t(S)) / 2, symmetric = TRUE)
        drop(e$vectors %*% (t(e$vectors) / sqrt(e$values)) %*% cn$L %*% beta)
      }))
    }
    Jpsi_all <- numDeriv::jacobian(Rfun, psi_hat)
    H <- numDeriv::hessian(function(psi) .mixedCiDeviance(base, psi), psi_hat)
    eh <- eigen((H + t(H)) / 2, symmetric = TRUE)
    keep <- eh$values > max(eh$values) * 1e-8
    Vpsi <- 2 * eh$vectors[, keep, drop = FALSE] %*%
      (t(eh$vectors[, keep, drop = FALSE]) / eh$values[keep])
  }

  ends <- cumsum(vapply(contrasts, function(x) nrow(x$L), integer(1L)))
  cts <- lapply(seq_along(contrasts), function(i) {
    cn <- contrasts[[i]]
    idx <- (ends[i] - nrow(cn$L) + 1L):ends[i]
    Jpsi <- if (robust) NULL else Jpsi_all[idx, , drop = FALSE]
    .mixedCiContrast(cn$L, cn$n, beta, Vhat, gls$info, pieces, robust,
                     design, cn$counts, Jpsi, if (robust) NULL else Vpsi)
  })
  out <- lapply(setNames(ci.method, ci.method), function(method) {
    res <- lapply(seq_along(cts), function(i) {
      .mixedCiInterval(cts[[i]], method, alpha, signed = i <= p)
    })
    tab <- list(coefficients = do.call(rbind, res[seq_len(p)]))
    rownames(tab$coefficients) <- names(beta)
    if (!is.null(n_anova)) {
      tab$anova <- do.call(rbind, res[-seq_len(p)])
      rownames(tab$anova) <- names(n_anova)
    }
    tab
  })
  out$se <- vapply(cts, function(ct) if (ct$m1 == 1L) sqrt(ct$var_R) else NA_real_, numeric(1L))
  names(out$se) <- c(names(beta), names(n_anova))
  out
}
