#######
# Function to construct CI for
# RESI in generalized linear models
######

#' Function generating CIs for the RESI via various bootstrap methods
#' @export
#' @param model.full the full model. It should be a `glm` object.
#' @param model.reduced the reduced `glm()` model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param r numeric, the number of boostrap replicates. By default, 1000 bootstraps will be implemented.
#' @param robust.var default to TRUE, whether to use the robust (sandwich) variance estimator when construct the Wald test statistic. If `TRUE`, the variance of the estimator will be obtained by using `sandwich::vcovHC()`` and the HC3 will be applied.
#' @param sigma2 the true variance of error under homoskedasticity assumption (added for simulations).
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use. 1: resampling covariates along with residuals (default); 2: fixing covariates and only bootstrapping residulas; 3: resampling covariates and residuals independently w/ replacements; 4. no sampling, just multipliers
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used0
#' @param correct for the linear regression models (i.e., `family = 'gaussian'` in the `glm()` function) whether the residuals with bias correction will be used, by default, FALSE.
#' @param num.cores The number of CPU cores to be used for calculating bootstrapped CIs, by default, only 1 core will be used.
#' @param digits the number of decimal digits in the output ANOES table. By default, 3

boot.ci <- function(model.full, model.reduced = NULL, r = 1000, robust.var = TRUE, multi = 'none', boot.type = 1, alpha = 0.05, correct = FALSE, num.cores = 1, digits = 3){

  # data frame
  data = model.full$model
  # model forms
  form.full <- formula(model.full)
  if (is.null(model.reduced)){
    form.reduced = as.formula(paste0(all.vars(form.full)[1], " ~ 1"))
  } else {
    form.reduced <- formula(model.reduced)
  }

  # glm family
  family = model.full$family$family

  # covariance estimator
  if (robust.var) {
    vcovfunc = sandwich::vcovHC
    } else {
      vcovfunc = vcov
    }

  # residual estimates (corrected or not)
  # I think this correction is for linear models??
  if (correct == TRUE) {
    data$resid = residuals(model.full, type = "response")/(1-hatvalues(model.full))
  } else {
    data$resid = residuals(model.full, type = "response")
  }

  # variable names of interest
  # var.names = all.vars(form.full)[-1]
  var.names = names(coef(model.full))
  var.names = var.names[var.names != "(Intercept)"]

  # Bootstraps
  temp <- simplify2array(parallel::mclapply(1:r,
                                            function(ind, boot.data, form.full, form.reduced, multi) {

                                              temp.data = boot.data # the input data

                                              repeat {
                                                boot.data = temp.data
                                                # version 1: resample X along with residuals non-parametrically
                                                if (boot.type == 1) {
                                                  boot.ind <- sample(1:nrow(boot.data), replace = TRUE)
                                                  boot.data <- boot.data[boot.ind, ]
                                                }
                                                # version 2: fix X and only resample residuals (w/ replacement)
                                                if (boot.type == 2) {
                                                  boot.ind <- sample(1:nrow(boot.data), replace = TRUE)
                                                  resid <- boot.data[boot.ind, 'resid']
                                                  boot.data <- cbind(boot.data[, names(boot.data) != 'resid'], resid)
                                                }
                                                # version 3: resampling covariates and residuals independently w/ replacements
                                                if (boot.type == 3) {
                                                  # bootstrap indicator
                                                  boot.ind1 <- sample(1:nrow(boot.data), replace = TRUE)
                                                  boot.ind2 <- sample(1:nrow(boot.data), replace = TRUE)
                                                  resid <- boot.data[boot.ind2, 'resid'] # bootstrapped residuals
                                                  boot.data <- cbind(boot.data[boot.ind1, names(boot.data) != 'resid'],
                                                                     resid)
                                                }
                                                # version 4: no resampling, just multipliers
                                                if (boot.type == 4){
                                                  boot.data = boot.data
                                                }
                                                # to detect whether the re-sample covariates have constant value, if yes, skip to next loop
                                                f = function(x) {length(unique(x))}
                                                if (all(apply(X = boot.data, MARGIN = 2, FUN = f)  > 1) ) break
                                              }

                                              # obtain multiplers
                                              n = nrow(boot.data)
                                              if (multi == "none") w <- rep(1, n)
                                              if (multi == "rad") w <- sample(x = c(-1, 1), size = n, replace = TRUE)
                                              if (multi == "normal") w <- rnorm(n)

                                              # calculate the bootstrapped outcome values
                                              ## note: replace the observed y with wild-bootstrapped y (so that I don't need to re-specify the model formula)
                                              boot.data$y <- predict(model.full, newdata = boot.data, type = "response") + boot.data$resid * w

                                              # fit glm on the bootstrapped data
                                              ## full model
                                              boot.model.full <- glm(form.full, data = boot.data, family = family)
                                              ## reduced model
                                              boot.model.reduced <- glm(form.reduced, data = boot.data, family = family)

                                              ## Overall (Wald) test stat
                                              wald.test = lmtest::waldtest(boot.model.reduced, boot.model.full, vcov = vcovfunc, test = 'Chisq')
                                              stats = wald.test$Chisq[2]
                                              names(stats) <- "Overall"

                                              # when `model.reduced = NULL`, output the wald test stat for each effect
                                              if (is.null(model.reduced)){
                                                ## Individual (Wald) test stats
                                                ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
                                                ind.stats <- (coef(boot.model.full)[var.names])^2/diag(as.matrix(vcovfunc(boot.model.full)[var.names, var.names]))
                                                stats <- c(stats, ind.stats)
                                                names(stats) <- c("Overall", var.names)
                                              }
                                              stats
                                            } # end of function in `lapply`
                                            ,
                                            boot.data = data,
                                            form.full = form.full,
                                            form.reduced = form.reduced,
                                            multi = multi,
                                            mc.cores = num.cores)
  )

  temp2 = temp # temp still has the row names!
  dim(temp2) = c(length(temp2)/r, r)
  boot.stats = t(temp2)
  # use the row names of `temp`
  colnames(boot.stats) = if(is.null(rownames(temp))) "Overall" else {rownames(temp)}


  # Deriving the point estimates for RESI
  # fit glm on the data
  ## full model
  mod1 <- glm(form.full, data = data, family = family)
  ## reduced model
  mod0 <- glm(form.reduced, data = data, family = family)

  ## Overall (Wald) test stat
  wald.test = lmtest::waldtest(mod0, mod1, vcov = vcovfunc, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]

  overall.resi.hat = ifelse(robust.var, chisq2S(stats, overall.df, res.df), f2S(stats, overall.df, res.df))

  # set up ANOES table
  anoes.tab = matrix(rep(NA, 2*6), nrow = 2, ncol = 6)
  anoes.tab[1, c(1, 2, 4) ] = c(stats, overall.df, overall.resi.hat)
  anoes.tab[2, 2] = res.df
  rownames(anoes.tab) = c("Tested", "Residual")

  # when `model.reduced = NULL`, output the wald test stat for each effect
  if (is.null(model.reduced)){
    # set up ANOES table
    anoes.tab = matrix(rep(NA, 2*8), nrow = 2, ncol = 8)
    anoes.tab[1, c(3, 4, 6) ] = c(stats, overall.df, overall.resi.hat)
    anoes.tab[2, 4] = res.df
    ## Individual (Wald) test stats
    ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
    ind.est <- coef(model.full)[var.names] # coef estiamtes
    ind.rob.se <- sqrt(diag(as.matrix(vcovfunc(model.full)[var.names, var.names]))) # corresponding robust s.e.
    ind.stats <- (ind.est/ind.rob.se)^2 # wald test statistics
    ind.resi.hat <- if(robust.var) {chisq2S(ind.stats, 1, res.df)} else {f2S(ind.stats, 1, res.df)}
    rownames(anoes.tab) = c("Overall", "Residual")
    anoes.tab = rbind(cbind(ind.est, ind.rob.se, ind.stats, 1, NA, ind.resi.hat, NA, NA), anoes.tab)
    colnames(anoes.tab) = c("estimate", "robust s.e.", "Chi-squared", "df", "p-val", "RESI", "LL", "UL")
  }
  colnames(anoes.tab) = c("Chi-squared", "df", "p-val", "RESI", "LL", "UL")
  # calculate p-values form Chi-sq dist
  anoes.tab[, 'p-val'] = pchisq(anoes.tab[, "Chi-squared"], df = anoes.tab[, "df"], lower.tail = FALSE)

  # Deriving the bootstrap CI for the RESI
  boot.S <- boot.stats # define boot.S which has the same frame as `boot.stats`
  ## converting statistics to RESI
  for (i in 1:ncol(boot.stats)) {
    boot.S[, i] <- if(robust.var) {chisq2S(boot.stats[, i], 1, res.df)} else {f2S(boot.stats[ , i], 1, res.df)}
  } # end of loop `i`
  CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2))
  CIs = t(CIs)
  if (nrow(CIs) == 1) rownames(CIs) = "Tested"
  anoes.tab[rownames(CIs), c("LL", "UL")] = CIs

  output <- list(model.full = formula(model.full),
                 model.reduced = form.reduced,
                 model.family = family,
                 boot.type = boot.type,
                 multiplier = multi,
                 alpha = alpha,
                 `number of bootstraps` = r,
                 correct = correct,
                 ANOES = round(anoes.tab, digits = digits)
                 )
  return(output)
  }

