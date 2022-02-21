#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' @param model.full the full model. It should be a `glm` object.
#' @param model.reduced the reduced `glm()` model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for linear models with splines)
#' @param anova whether to compute an ANOES table based on the Anova function. By default = `TRUE`
#' @param summary whether to produce a summary output with the RESI columns added. By default = `TRUE`
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use. 1: resampling covariates along with residuals (default); 2: fixing covariates and only bootstrapping residulas; 3: resampling covariates and residuals independently w/ replacements; 4. no sampling, just multipliers
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used0
#' @param correct for the linear regression models (i.e., `family = 'gaussian'` in the `glm()` function) whether the residuals with bias correction will be used, by default, FALSE.
#' @param num.cores The number of CPU cores to be used for calculating bootstrapped CIs, by default, only 1 core will be used.
#' @param ... Other arguments to be passed to Anova function
#' @importFrom car Anova
#' @importFrom lmtest waldtest
#' @importFrom sandwich vcovHC
#' @importFrom stats coef formula glm hatvalues pf predict quantile residuals update vcov
#' @export
#' @return


anoes <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE, nboot = 1000, vcovfunc = sandwich::vcovHC, multi = 'none', boot.type = 1, alpha = 0.05, correct = FALSE, num.cores = 1, ...){
  #browser()
  options(warn = 1)

  # check to make sure anova or summary is specified
  if (!anova & !summary){
    warning('no RESI output table (summary or anova) specified')
  }

  # data.frame
  if (missing(data)){
    if (class(model.full)[1] == 'glm'){
      data = model.full$data
    }
    else{
      data = model.full$model
    }
  }


  # model forms
  form.full <- formula(model.full)
  if (is.null(model.reduced)){
    form.reduced = as.formula(paste0(all.vars(form.full)[1], " ~ 1"))
    mod.reduced <- update(model.full, formula = form.reduced, data = model.full$model)
  } else {
    form.reduced <- formula(model.reduced)
    mod.reduced <- model.reduced
  }

  # glm family
  family = model.full$family$family

  # residual estimates (corrected or not)
  # I think this correction is for linear models??
  if (correct == TRUE) {
    data$resid = residuals(model.full, type = "response")/(1-hatvalues(model.full))
  } else {
    data$resid = residuals(model.full, type = "response")
  }

  # variable names
  var.names = names(coef(model.full))

  # Bootstraps
  temp <- simplify2array(parallel::mclapply(1:nboot,
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
                                                if (any(apply(X = boot.data, MARGIN = 2, FUN = f)  > 1) ) break
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
                                              boot.model.full <- update(model.full, data = boot.data)
                                              ## reduced model
                                              if (is.null(model.reduced)){
                                                boot.model.reduced = update(model.full, formula = form.reduced,  data = boot.data)
                                              }
                                              else{
                                                boot.model.reduced = update(mod.reduced,  data = boot.data)
                                              }


                                              ## Overall (Wald) test stat
                                              wald.test = lmtest::waldtest(boot.model.reduced, boot.model.full, vcov = vcovfunc, test = 'Chisq')
                                              stats = wald.test$Chisq[2]
                                              names.stats <- "Overall"
                                              names(stats) <- names.stats

                                              # when `summary = TRUE`, output the wald test stat for each effect
                                              if (summary){
                                                ## Individual (Wald) test stats
                                                ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
                                                ind.stats <- (coef(boot.model.full)[var.names])^2/diag(as.matrix(vcovfunc(boot.model.full)[var.names, var.names]))
                                                stats <- c(stats, ind.stats)
                                                names.stats <- c("Overall", var.names)
                                                names(stats) <- names.stats

                                              }
                                              # if (`model.reduced = NULL`) & `anova = TRUE`, calculate Wald test stat for each term
                                              if (is.null(model.reduced) & anova){
                                                if (class(model.full)[1] == 'lm'){
                                                  suppressMessages(anova.tab <- Anova(boot.model.full, vcov. = vcovfunc, ...))
                                                  term.stats <- anova.tab[,'F']
                                                }
                                                else{
                                                  suppressMessages(anova.tab <- Anova(boot.model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
                                                  term.stats <- anova.tab[,'Chisq']
                                                }

                                                stats <- c(stats, term.stats)
                                                names(stats) <- c(names.stats, rownames(anova.tab))
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
  dim(temp2) = c(length(temp2)/nboot, nboot)
  boot.stats = t(temp2)
  # use the row names of `temp`
  colnames(boot.stats) = if(is.null(rownames(temp))) "Overall" else {rownames(temp)}


  # Deriving the point estimates for RESI

  ## Overall (Wald) test stat
  wald.test = lmtest::waldtest(mod.reduced, model.full, vcov = vcovfunc, test = 'Chisq')
  stats = wald.test$Chisq[2]
  overall.df = wald.test$Df[2]
  res.df = wald.test$Res.Df[2]

  # not sure if this works right with other specifications of robust variance.
  # overall.resi.hat = ifelse(isTRUE(all.equal(vcovfunc, sandwich::vcovHC)), chisq2S(stats, overall.df, res.df), f2S(stats, overall.df, res.df))
  overall.resi.hat = chisq2S(stats, overall.df, res.df)

  # overall RESI CI
  # overall.stats = if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))) {chisq2S(boot.stats[,1], overall.df, res.df)} else {f2S(boot.stats[,1], overall.df, res.df)}
  overall.stats = chisq2S(boot.stats[,1], overall.df, res.df)
  overall.ci = quantile(overall.stats, probs = c(alpha/2, 1-alpha/2))


  # compute ANOVA-style table if `anova = TRUE`
  if (is.null(model.reduced) & anova){
    # set up ANOES tab based on Anova
    if (class(model.full)[1] == 'lm'){
      suppressMessages(anoes.tab <- Anova(model.full, vcov. = vcovfunc, ...))
      terms.resi.hat = f2S(anoes.tab[,'F'], anoes.tab[,'Df'], res.df)
      ## converting statistics to RESI
      boot.S <- boot.stats[,c(((ncol(boot.stats) - nrow(anoes.tab) + 1):ncol(boot.stats)))]
      for (i in 1:ncol(boot.S)) {
        boot.S[, i] <- f2S(boot.S[ , i], anoes.tab[i, 'Df'], res.df)
      } # end of loop `i`
    }
    else{
      suppressMessages(anoes.tab <- Anova(model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
      terms.resi.hat = chisq2S(anoes.tab[,'Chisq'], anoes.tab[,'Df'], res.df)
      boot.S <- boot.stats[,c(((ncol(boot.stats) - nrow(anoes.tab) + 1):ncol(boot.stats)))]
      for (i in 1:ncol(boot.S)) {
        boot.S[, i] <- chisq2S(boot.S[ , i], anoes.tab[i, 'Df'], res.df)
      } # end of loop `i`
    }

    anoes.tab = rbind(cbind(anoes.tab, RESI = terms.resi.hat), Overall = c(overall.df, stats, pchisq(stats, df = overall.df, lower.tail = FALSE), overall.resi.hat))

    CIs = apply(boot.S[,which(colnames(boot.S) != "Residuals")], 2,  quantile, probs = c(alpha/2, 1-alpha/2))
    CIs = t(CIs)
    anoes.tab[rownames(CIs), c("LL", "UL")] = CIs
    anoes.tab['Overall', c("LL", "UL")] = overall.ci
  }

  # compute summaryES table if `summary = TRUE`
  if (summary){
    summaryES.tab = matrix(nrow = 0, ncol = 8)

    ## Individual (Wald) test stats
    ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
    ind.est <- coef(model.full)[var.names] # coef estiamtes
    ind.se <- sqrt(diag(as.matrix(vcovfunc(model.full)[var.names, var.names]))) # corresponding (robust) s.e.
    ind.stats <- (ind.est/ind.se)^2 # wald test statistics
    #ind.resi.hat <- if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))) {chisq2S(ind.stats, 1, res.df)} else {f2S(ind.stats, 1, res.df)}
    ind.resi.hat <- chisq2S(ind.stats, 1, res.df)

    summaryES.tab = rbind(cbind(ind.est, ind.se,ind.stats, 1, NA, ind.resi.hat, NA, NA), summaryES.tab)
    rownames(summaryES.tab) = c(var.names)
    colnames(summaryES.tab) = c("estimate", "robust s.e.", "Chi-squared", "df", "p-val", "RESI", "LL", "UL")
    if (is.null(model.reduced)){
      summaryES.tab = rbind(summaryES.tab, c(NA, NA, stats, overall.df, NA, overall.resi.hat, overall.ci))
      rownames(summaryES.tab) = c(var.names, "Overall")
    }

    # Deriving the bootstrap CI for the RESI
    if (ncol(boot.stats) > (1 + length(var.names))){
      boot.S <- boot.stats[, 2:(ncol(boot.stats) - nrow(anoes.tab) + 1)]
      ## converting statistics to RESI
      for (i in 1:ncol(boot.S)) {
        #boot.S[, i] <- if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))){chisq2S(boot.S[, i], 1, res.df)} else {f2S(boot.S[ , i], 1, res.df)}
        boot.S[,i] <- chisq2S(boot.S[, i], 1, res.df)
      } # end of loop `i`
      CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2))
    }
    else{
      boot.S <- boot.stats[, 2:ncol(boot.stats)] # define boot.S which has the same frame as `boot.stats`
      ## converting statistics to RESI
      for (i in 1:ncol(boot.S)) {
        #boot.S[, i] <- if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))) {chisq2S(boot.S[, i], 1, res.df)} else {f2S(boot.S[ , i], 1, res.df)}
        boot.S[,i] <- chisq2S(boot.S[, i], 1, res.df)
      } # end of loop `i`
      CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2))
    }
    CIs = t(CIs)
    summaryES.tab[rownames(CIs), c("LL", "UL")] = CIs
    summaryES.tab[, 'p-val'] = pchisq(summaryES.tab[, "Chi-squared"], df = summaryES.tab[, "df"], lower.tail = FALSE)
  }


  output <- list(model.full = formula(model.full),
                 model.reduced = form.reduced,
                 model.family = family,
                 boot.type = boot.type,
                 multiplier = multi,
                 alpha = alpha,
                 `number of bootstraps` = nboot,
                 correct = correct)

  if (summary){
    output$summary = summaryES.tab
  }

  if (is.null(model.reduced) & anova){
    output$anova = anoes.tab
  }

  if (!is.null(model.reduced)){
    wald.test[2, c("RESI", "LL", "UL")] = c(overall.resi.hat, overall.ci)
    output$waldtest = wald.test
  }

  class(output) = "anoes"

  return(output)
}
