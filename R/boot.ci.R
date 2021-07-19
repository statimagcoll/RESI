#######
# Function to construct CI for
# RESI in linear models via Wild Bootstrap or non-parametric bootstrap .
######

#' function generating CIS via bootstrap
#' @export
#' @param model.full the full `lm()` model
#' @param model.reduced the reduced `lm()` model to compare with
#' @param r numeric, the number of boostrap replicates. By default, 1000 bootstraps will be implemented.
#' @param method method used to compute the statistics.
#' "F" corresponds to using var-cov estimator `vcov`, and the resulting statistic will be F-statistic;
#' "Chisq" corresponds to using "sandwich::vcovHC" and results a Chi-squared statistics;
#' "Known" corresponds to knowing the true error of the linear model with homoskedasticity. With method = "Known", sigma2 need to be specified.
#' @param sigma2 the true variance of error under homoskedasticity assumption.
#' @param multi the distribution from which the multipliers will be drawn: 'none' = the multipliers equal constant 1 (default); 'rad' = rademacher; 'normal' = Std Normal distribution
#' @param boot.type which type of bootstrap to use.
# 1: resampling covariates along with residuals;
# 2: fixing covariates and only bootstrapping residulas;
# 3: resampling covariates and residuals independently w/ replacements
# 4. no sampling, just multipliers
#' @param alpha significance level used to compute the CIs. By default, 0.05 will be used
#' @param correct whether the residuals with bias correction will be used, by default, TRUE.
#' @param num.cores The number of CPU cores to be used for calculating bootstrapped CIs, by default, only 1 core will be used.

# Note:
# 1. should we test the interaction terms?
# 2. Do we really need reduced model? The ANOVA gives the results that the additional variables account for xx effect size in addition to ....(other variables)
# 3. if we need a reduced model and it's not empty, how are the residual dfs defined? residual df of the full/reduced model?


boot.ci <- function(model.full, model.reduced = NULL, r = 1000, method = "F", multi = 'none', boot.type = 1, alpha = 0.05, correct = TRUE, num.cores = 1){

  # data frame
  data = model.full$model
  # model forms
  form.full <- formula(model.full)
  if (is.null(model.reduced)){
    form.reduced = NULL
  } else {
    form.reduced <- formula(model.reduced)
  }



  # check whether the two models are nested
  # if (!(length((setdiff(all.vars(form.full), all.vars(form.reduced)))) > 0 &
  # length(setdiff(all.vars(form.reduced), all.vars(form.full))) == 0 )) {
  #   stop("The models are not nested")
  # }

  # covariance estimator
  if (tolower(method) == 'f'){ # assuming homoskedasticity
    vcovfunc <- vcov
  } else {
    if (tolower(method) == "chisq"){ # using robust var-cov estimator (allowing heteroskedasticity)
      vcovfunc <- sandwich::vcovHC
    }  else {stop("Please correctly specify the method to be used to compute the var-cov estimate (i.e., 'F' or 'Chisq')")}
  }

  S_hat = NULL

  # residual estimates (corrected or not)
  if (correct == TRUE) {
    data$resid = residuals(model.full)/(1-hatvalues(model.full))
  } else {
    data$resid = residuals(model.full)
  }

  # # The variables of interest (difference between full & reduced models)
  # var.list <- setdiff(all.vars(form.full), all.vars(form.reduced))

  # bootstraps
  temp <- simplify2array(parallel::mclapply(1:r,
                       function(ind, boot.data, form.full, form.reduced, multi, vcovfunc) {

                         # needs further modification. for now, only works for 1 covariates whose name is `x`
                         repeat {
                           # version 1: resample X along with residuals non-parametrically
                           if (boot.type == 1) {
                             boot.ind <- sample(1:nrow(boot.data), replace = TRUE)
                             boot.data <- boot.data[boot.ind, ]
                           }
                           # version 2: fix X and only resample residuals (w/ replacement)
                           if (boot.type == 2) {
                             # bootstrap indicator
                             boot.ind <- sample(1:nrow(boot.data), replace = TRUE)
                             resid <- boot.data[boot.ind, 'resid']
                             boot.data <- cbind(boot.data[, names(boot.data) != 'resid'],
                                                resid)
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

                           if (unique(boot.data$x) > 1) break
                         }



                         # obtain multiplers
                         n = nrow(boot.data)
                         if (multi == "none") w <- rep(1, n)
                         if (multi == "rad") w <- sample(x = c(-1, 1), size = n, replace = TRUE)
                         if (multi == "normal") w <- rnorm(n)

                         # calculate the bootstrapped outcome values
                         ## note: replace the observed y with wild-bootstrapped y (so that I don't need to re-specify the model formula)
                         boot.data$y <- predict(model.full, newdata = boot.data) + boot.data$resid * w

                         # fit lm on the bootstrapped data
                         ## full model
                         boot.model.full <- lm(form.full, data = boot.data)
                         ## reduced model
                         if(!is.null(model.reduced)){
                           boot.model.reduced <- lm(form.reduced, data = boot.data)
                           # individual effects
                           # individual.anova <- suppressMessages(car::Anova(boot.model.full, test.statistics = "Wald", vcov. = vcovfunc))
                           # individual.stats <- individual.anova$F
                           # names(individual.stats) <- rownames(individual.anova)
                           #
                           # individual.stats <- individual.stats[var.list]
                           # overall effect
                           overall.stats <- lmtest::waldtest(boot.model.reduced, boot.model.full, vcov =  vcovfunc, test = "Chisq")$Chisq[2]
                           names(overall.stats) <- "Overall"
                           stats <- c(overall.stats)
                           stats
                         } else {
                           boot.model.reduced <- lm(as.formula(paste0(all.vars(form.full)[1], " ~ 1")), data = boot.data) # y ~ 1
                           # directly compute the statistics based on full model
                           # individual effects
                           individual.anova <- suppressMessages(car::Anova(boot.model.full, test.statistics = "Wald", vcov. = vcovfunc))
                           individual.stats <- individual.anova$F
                           names(individual.stats) <- rownames(individual.anova)
                           individual.stats <- individual.stats[names(individual.stats) != "Residuals"]
                           # overall effect
                           overall.stats <- lmtest::waldtest(boot.model.reduced, boot.model.full, vcov = vcovfunc(boot.model.full), test = "Chisq")$Chisq[2]
                           names(overall.stats) <- "Overall"
                           stats <- c(individual.stats, overall.stats)
                           stats
                         }
                       } # end of function in `lapply`
                       , boot.data = data, form.full = form.full, form.reduced = form.reduced, multi = multi, vcovfunc = vcovfunc, mc.cores = num.cores)
                      )
  temp2 = temp # temp still has the row names!
  dim(temp2) = c(length(temp2)/r, r)
  boot.stats = as.data.frame(t(temp2))
  # use the row names of `temp`
  colnames(boot.stats) = if(is.null(rownames(temp))) "Overall" else {rownames(temp)}

  # calculate the statistics for the data
  if (!is.null(form.reduced)){
    # ### individual effects
    # original.anova <- suppressMessages(car::Anova(model.full, test.statistics = "Wald", vcov. = vcovfunc))
    ### overall effect
    original.overall <- lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc, test = "Chisq")
    ### d.f
    dfs <- c(original.overall$Df[2], original.overall$Res.Df[2])
    names(dfs) <- c("Overall", "Residuals")
    # names(dfs) <- rownames(original.overall)
    # overall.df <- sum(dfs[names(dfs) != "Residuals"])
    # dfs <- c(dfs, Overall = overall.df)
    anova.tab = data.frame(Df = dfs,
                           Stat = c(original.overall$Chisq[2], NA),
                           'P-value' = c(original.overall$`Pr(>Chisq)`[2], NA))
  }  else {
    ### individual effects
    original.anova <- suppressMessages(car::Anova(model.full, test.statistics = "Wald", vcov. = vcovfunc))
    ### overall effect
    model.reduced <- lm(as.formula(paste0(all.vars(form.full)[1], " ~ 1")), data = data)
    original.overall <- lmtest::waldtest(model.reduced, model.full, vcov = vcovfunc, test = "Chisq")
    ### d.f
    dfs <- original.anova$Df
    names(dfs) <- rownames(original.anova)
    overall.df <- original.overall$Df[2]
    dfs <- c(dfs, Overall = overall.df)

    anova.tab = data.frame(Df = dfs,
                           Stat = c(original.anova$F, original.overall$Chisq[2]),
                           'P-value' = c(original.anova$`Pr(>F)`, original.overall$`Pr(>Chisq)`[2]))

    # list of effects of interest
    effect.list <- all.vars(form.full)[-1]
    anova.tab = anova.tab[c(effect.list, "Overall", "Residuals"), ]
  }

  ## conversion from stats to RESI
  ### define conversion function
  if (tolower(method) == 'chisq') { # when they are Chisq statistics
    stat2S <- function (chisq, df, rdf) {
      S = (chisq - df)/rdf
      sqrt(ifelse(S<0, 0, S))
    }
  } else { # when it's F statistic
    stat2S <- function (f, df, rdf) {
      S = (f*df*(rdf-2)/rdf - df)/rdf
      sqrt(ifelse(S<0, 0, S))
    }
  }

  boot.S <- boot.stats # define boot.S which has the same frame as `boot.stats`

  for (i in 1:ncol(boot.stats)) {
    boot.S[, i] <- stat2S(boot.stats[, i], dfs[names(boot.stats)[i]], dfs['Residuals'])
  } # end of loop `i`

    # compute the quantiles for each column (i.e., for each effect)
    CIs = apply(boot.S, 2,  quantile, probs = c(alpha/2, 1-alpha/2))
    CIs = t(CIs)
    # Stats <- c(original.anova$F, original.overall$Chisq[2])
    # p_values <- c(original.anova$`Pr(>F)`, original.overall$`Pr(>Chisq)`[2])
    dfs = anova.tab$Df
    names(dfs) = rownames(anova.tab)
    Ss <- stat2S(anova.tab$Stat, dfs, dfs['Residuals'])
    anova.tab <- cbind(anova.tab, RESI = Ss, rbind(CIs,  NA))

    # colnames(anova.tab) = c("df", "Stat", "P value", "RESI", paste0(alpha/2*100, "% ", c("LB", "UB")) )
    output = list(model.full = form.full, model.reduced = form.reduced,
                  boot.type = boot.type, multiplier = multi, method = method,
                  alpha = alpha, correct = correct, anova = round(anova.tab, 3))
    return(output)
}









