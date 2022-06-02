#' Analysis of Effect Sizes (ANOES) based on the Robust Effect Size index (RESI) for (generalized) linear regression models
#' @param model.full the full model.
#' @param model.reduced the reduced model to compare with the full model. By default `NULL`, it's the same model as the full model but only having intercept.
#' @param data optional data frame or object coercible to data frame of model.full data (required for some model types)
#' @param anova whether to compute an ANOES table based on the Anova function. By default = `TRUE`
#' @param summary whether to produce a summary table with the RESI columns added. By default = `TRUE`
#' @param nboot numeric, the number of bootstrap replicates. By default, 1000 bootstraps will be implemented.
#' @param vcovfunc the variance estimator function for constructing the Wald test statistic. By default, sandwich::vcovHC (the robust (sandwich) variance estimator)
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use. 1: resampling covariates along with residuals (default); 2: fixing covariates and only bootstrapping residulas; 3: resampling covariates and residuals independently w/ replacements; 4. no sampling, just multipliers
#' @param alpha significance level of the constructed CIs. By default, 0.05 will be used.
#' @param correct for the linear regression models (i.e., `family = 'gaussian'` in the `glm()` function) whether the residuals with bias correction will be used, by default, FALSE.
#' @param num.cores The number of CPU cores to be used for calculating bootstrapped CIs, by default, only 1 core will be used.
#' @param ... Other arguments to be passed to Anova function
#' @export
#' @return Returns a list of arguments and RESI estimates and confidence intervals


anoes <- function(model.full, model.reduced = NULL, data, anova = TRUE, summary = TRUE, nboot = 1000, vcovfunc = sandwich::vcovHC, multi = 'none', boot.type = 1, alpha = 0.05, correct = FALSE, num.cores = 1, ...){
  browser()

  # check to make sure anova or summary is specified
  if (!anova & !summary){
    warning('no RESI output table (summary or anova) specified')
  }

  # check for Anova method
  Anova.method = ifelse(paste("Anova", class(model.full)[1], sep = ".") %in% methods(class = class(model.full)[1]), TRUE, FALSE)

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
  cl <- class(model.full)[1]
  form.full <- formula(model.full)
  if (is.null(model.reduced)){
    if (!(cl %in% c("survreg", "coxph"))){
      form.reduced = as.formula(paste0(all.vars(form.full)[1], " ~ 1"))
    }
    else{
      form.reduced = as.formula(paste(format(formula(model.full)[[2]]), "~ 1"))
    }
    if (cl != "nls"){
      mod.reduced <- update(model.full, formula = form.reduced, data = data)
    }
  } else {
    form.reduced <- formula(model.reduced)
    mod.reduced <- model.reduced
  }


  # glm family
  family = model.full$family$family

  # residual estimates (corrected or not)
  # I think this correction is for linear models??
  if (cl != "coxph"){
    if (correct == TRUE) {
      data$resid = residuals(model.full, type = "response")/(1-hatvalues(model.full))
    } else {
      res = residuals(model.full, type = "response")
      # check if can just make this [,'resid'] (thinking about missing data)
      if (!is.null(names(res))){
        data[names(res), "resid"] = res
      }
      else{
        data[,"resid"] = res
      }
    }
  }

  # variable names
  var.names = names(coef(model.full))

  # default arg for vcovfunc will not work with (all?) nls
  if (cl == "nls" & identical(vcovfunc, sandwich::vcovHC)){
    vcovfunc = regtools::nlshc
    message("Default vcov function not applicable for model type, vcovfunc set to regtools::nlshc")
  }

  if (cl %in% c("survreg", "coxph")){
    # robust option in model setup
    vcovfunc = vcov
  }

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
                                              if (cl != "coxph"){
                                                boot.data[,all.vars(form.full)[1]] <- predict(model.full, newdata = boot.data, type = "response") + boot.data$resid * w
                                                # colnames(boot.data)[colnames(boot.data) == all.vars(form.full)[1]] = "old y"
                                                # colnames(boot.data)[colnames(boot.data) == "y"] = all.vars(form.full)[1]
                                              }

                                              # fit model on the bootstrapped data
                                              ## full model
                                              boot.model.full <- update(model.full, data = boot.data)

                                              vcovmat = vcovfunc(boot.model.full)
                                              # nls name correction
                                              if (identical(vcovfunc, nlshc)){
                                                rownames(vcovmat) = var.names
                                                colnames(vcovmat) = var.names
                                              }

                                              ## reduced model
                                              if (cl %in% c("glm", "lm")){
                                                if (is.null(model.reduced)){
                                                  boot.model.reduced = update(model.full, formula = form.reduced,  data = boot.data)
                                                }
                                                else{
                                                  boot.model.reduced = update(mod.reduced,  data = boot.data)
                                                }
                                                ## Overall (Wald) test stat
                                                wald.test = lmtest::waldtest(boot.model.reduced, boot.model.full, vcov = vcovfunc, test = 'Chisq')
                                                stats = wald.test$Chisq[2]
                                              }
                                              else{
                                                wald.test = wald.test(vcovmat[var.names, var.names], coef(boot.model.full), Terms = 2:length(coef(boot.model.full)))
                                                stats = wald.test$result$chi2["chi2"]
                                              }
                                              names.stats <- "Overall"
                                              names(stats) <- names.stats

                                              term.stats = NULL
                                              anova.tab = NULL
                                              # if (`model.reduced = NULL`) & `anova = TRUE` & there is an Anova method, calculate Wald test stat for each term
                                              if (is.null(model.reduced) & anova & Anova.method){
                                                if (cl == 'lm'){
                                                  suppressMessages(anova.tab <- Anova(boot.model.full, vcov. = vcovfunc, ...))
                                                  term.stats <- anova.tab[,'F']
                                                }
                                                else{
                                                  suppressMessages(anova.tab <- Anova(boot.model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
                                                  term.stats <- anova.tab[,'Chisq']
                                                }}

                                              # when `summary = TRUE`, output the wald test stat for each effect
                                              if (summary){
                                                ## Individual (Wald) test stats
                                                ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
                                                ind.stats <- (coef(boot.model.full)[var.names])^2/diag(as.matrix(vcovmat[var.names, var.names]))
                                                stats <- c(stats, ind.stats)
                                                names.stats <- c("Overall", var.names)
                                                names(stats) <- names.stats
                                              }

                                              stats <- c(stats, term.stats)
                                              names(stats) <- c(names.stats, rownames(anova.tab))
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

  vcovmat = vcovfunc(model.full)
  if (identical(vcovfunc, nlshc)){
    rownames(vcovmat) = var.names
    colnames(vcovmat) = var.names
  }

  ## Overall (Wald) test stat (chi-sq)
  if (cl %in% c("glm", "lm")){
    wald.test = lmtest::waldtest(mod.reduced, model.full, vcov = vcovfunc, test = 'Chisq')
    stats = wald.test$Chisq[2]
    overall.df = wald.test$Df[2]
    res.df = wald.test$Res.Df[2]
  }
  else{
    wald.test = wald.test(vcovmat[var.names, var.names], coef(model.full), Terms = 2:length(coef(model.full)))
    stats = wald.test$result$chi2["chi2"]
    overall.df = wald.test$result$chi2["df"]
    res.df = nrow(data) - overall.df
  }
  # overall RESI
  overall.resi.hat = chisq2S(stats, overall.df, res.df)
  # overall RESI CI
  overall.stats = chisq2S(boot.stats[,1], overall.df, res.df)
  overall.ci = quantile(overall.stats, probs = c(alpha/2, 1-alpha/2))


  # compute ANOVA-style table if `anova = TRUE`
  if (is.null(model.reduced) & anova & Anova.method){
    # set up ANOES tab based on Anova
    if (cl == 'lm'){
      suppressMessages(anoes.tab <- Anova(model.full, vcov. = vcovfunc, ...))
      terms.resi.hat = f2S(anoes.tab[,'F'], anoes.tab[,'Df'], res.df)
      ## converting statistics to RESI
      boot.S <- boot.stats[,c(((ncol(boot.stats) - nrow(anoes.tab) + 1):ncol(boot.stats)))]
      boot.S <- as.matrix(boot.S)
      for (i in 1:ncol(boot.S)) {
        boot.S[, i] <- f2S(boot.S[ , i], anoes.tab[i, 'Df'], res.df)
      } # end of loop `i`

    }
    else{
      suppressMessages(anoes.tab <- Anova(model.full, test.statistic = 'Wald', vcov. = vcovfunc, ...))
      terms.resi.hat = chisq2S(anoes.tab[,'Chisq'], anoes.tab[,'Df'], res.df)
      boot.S <- boot.stats[,c(((ncol(boot.stats) - nrow(anoes.tab) + 1):ncol(boot.stats)))]
      boot.S <- as.matrix(boot.S)
      for (i in 1:ncol(boot.S)) {
        boot.S[, i] <- chisq2S(boot.S[ , i], anoes.tab[i, 'Df'], res.df)
      } # end of loop `i`
    }

    anoes.tab = cbind(anoes.tab, RESI = terms.resi.hat)

    CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
    CIs = t(CIs)
    anoes.tab[1:nrow(CIs), c("LL", "UL")] = CIs
  }

  # compute summaryES table if `summary = TRUE`
  if (summary){
    summaryES.tab = matrix(nrow = 0, ncol = 8)

    ## Individual (Wald) test stats
    ## note: the results from car:Anova look wield when using glm()...so I calculate the stats manually
    ind.est <- coef(model.full)[var.names] # coef estiamtes
    ind.se <- sqrt(diag(as.matrix(vcovmat)[var.names, var.names])) # corresponding (robust) s.e.
    ind.stats <- (ind.est/ind.se)^2 # wald test statistics
    #ind.resi.hat <- if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))) {chisq2S(ind.stats, 1, res.df)} else {f2S(ind.stats, 1, res.df)}
    ind.resi.hat <- chisq2S(ind.stats, 1, res.df)

    summaryES.tab = rbind(cbind(ind.est, ind.se,ind.stats, 1, NA, ind.resi.hat, NA, NA), summaryES.tab)
    rownames(summaryES.tab) = c(var.names)
    colnames(summaryES.tab) = c("estimate", "s.e.", "Chi-squared", "df", "p-val", "RESI", "LL", "UL")

    # Deriving the bootstrap CI for the RESI
    if (ncol(boot.stats) > (1 + length(var.names))){
      boot.S <- boot.stats[, 2:(length(var.names) + 1)]
      boot.S <- as.matrix(boot.S)
      ## converting statistics to RESI
      for (i in 1:ncol(boot.S)) {
        #boot.S[, i] <- if(isTRUE(all.equal(vcovfunc, sandwich::vcovHC))){chisq2S(boot.S[, i], 1, res.df)} else {f2S(boot.S[ , i], 1, res.df)}
        boot.S[,i] <- chisq2S(boot.S[, i], 1, res.df)
      } # end of loop `i`
      CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2))
    }
    else{
      boot.S <- boot.stats[, 2:ncol(boot.stats)] # define boot.S which has the same frame as `boot.stats`
      boot.S <- as.matrix(boot.S)
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


  output <- list(model.full = form.full,
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

  if (is.null(model.reduced) & anova & Anova.method){
    output$anova = anoes.tab
  }

  if (cl %in% c("glm", "lm")){
    wald.test[2, c("RESI", "LL", "UL")] = c(overall.resi.hat, overall.ci)}
  else{
    wald.test$result$chi2[c("RESI", "LL", "UL")] = c(overall.resi.hat, overall.ci)
    wald.test = wald.test$result$chi2
  }
  output$overall = wald.test

  class(output) = "anoes"

  return(output)
}

