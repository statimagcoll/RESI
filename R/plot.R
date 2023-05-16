#' Plotting RESI Estimates and CIs
#'
#' This function uses base graphics to plot robust effect size (RESI) estimates and confidence intervals
#' from `resi`, `summary_resi`, and `anova_resi` objects.
#' @param x Object of `resi`, `summary_resi`, or `anova_resi` class
#' @param alpha Numeric, desired alpha level for confidence intervals
#' @param ycex.axis Numeric, scale specifically for the variable name labels
#' @param yaxis.args List, other arguments to be passed to \code{\link{axis}} for the y-axis
#' @param automar Logical, whether to automatically adjust the plotting margins to accommodate variable names. Default = `TRUE`
#' @param ... Other graphical parameters passed to \code{\link{plot}} and \code{\link{lines}}
#' @details This function creaes a forest-like plot with RESI estimates for each variable or factor.
#' The size of the left margin will be automatically adjusted (and returned to original after plotting)
#' unless `automar = FALSE`. Additional graphics parameters will be passed to the main
#' plot function, the confidence intervals. Arguments specifically for the y-axis (variable names)
#' can be specified using `yaxis.args`. To manually adjust the size of the y-axis labels without
#' affecting the x-axis, the user can specify a value for `ycex.axis`.
#' @return Returns a plot of RESI point estimates
#' @examples
#' # create a resi object
#' resi_obj <- resi(lm(charges ~ region * age + bmi + sex, data = RESI::insurance),
#' nboot = 10)
#'
#' # plot coefficients table, changing size of labels for both axes in the usual way
#' plot(resi_obj, cex.alpha = 0.7)
#'
#' # plot ANOVA table, changing the size of just the y-axis
#' plot(resi_obj, ycex.alpha = 0.8)
#' @importFrom graphics abline axis lines par strwidth
#' @export
plot.resi = function(x, alpha = NULL, ycex.axis = NULL, yaxis.args = list(),
                     automar = TRUE, ...){
  dots = list(...)

  if (!(is.null(x$alpha))){
    if (is.null(alpha)){
      alpha = x$alpha[1]
    }
    else{
      if (length(alpha) > 1){
        warning("\nOnly first alpha will be plotted")
        alpha = alpha[1]
      }
      if (!(alpha %in% x$alpha)){
        if (is.null(x$boot.results)){
          stop("\nSpecified alpha not found in the resi object")
        }
        CIs = apply(x$boot.results$t[,2:(1+nrow(x$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
        CIs = t(CIs)
        x$coefficients[1:nrow(CIs), c(paste(alpha/2*100, "%", sep=""), paste((1-rev(alpha)/2)*100, "%", sep=""))] = CIs
      }
    }
    ll = paste(alpha/2*100, "%", sep = "")
    ul = paste((1-alpha/2)*100, "%", sep = "")
  } else{
    ll = ul = "RESI"
  }

  vnames = rownames(x$coefficients)
  # width of labels
  if ((!(is.null(ycex.axis))) | "cex.axis" %in% names(dots)){
    if (!(is.null(ycex.axis))){
      ycex = ycex.axis
    } else{
      ycex = dots$cex.axis
    }
    w = strwidth(vnames, "inches", cex = ycex)
  } else{
    w = strwidth(vnames, "inches")
    ycex = 1
  }

  # save original margins
  omar = par("mai")
  if (automar){
    # set margins to accomodate labels
    par("mai" = omar + c(0, max(w) - omar[2] + 0.25, 0, 0))
  }

  er = 0.2 * (length(vnames)/(length(vnames) + 1))

  plot(x = x$coefficients[,"RESI"], y = length(vnames):1,
       xlim = c(min(0, min(x$coefficients[,ll])), max(x$coefficients[,ul])),
       xlab = "RESI Estimate", yaxt = "n", ylab = "",
       main = paste("Coefficient RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""),...)
  if (!(is.null(x$alpha))){
    for (i in length(vnames):1){
      lines(x = c(x$coefficients[-1*(i-length(vnames)) + 1,ll],
                  x$coefficients[-1*(i-length(vnames)) + 1,ul]), y = c(i,i),
            ...)
      lines(x = c(x$coefficients[-1*(i-length(vnames)) + 1,ll],
                  x$coefficients[-1*(i-length(vnames)) + 1,ll]),
            y = c(i-er,i+er),
            ...)
      lines(x = c(x$coefficients[-1*(i-length(vnames)) + 1,ul],
                  x$coefficients[-1*(i-length(vnames)) + 1,ul]),
            y = c(i-er,i+er),
            ...)}

  }
  do.call(axis, c(list(side = 2, at = length(vnames):1,
                       labels = vnames, las = 1,
                       cex.axis = ycex), yaxis.args))
  abline(v = 0, lty = 2)
  # return to original margins
  par("mai" = omar)
}

#' @export
plot.summary_resi = plot.resi

#' @export
plot.anova_resi = function(x, alpha = NULL, ycex.axis = NULL, yaxis.args = list(),
                           automar = TRUE, ...){
  dots = list(...)

  cols = grep("%", colnames(x))
  if (length(cols) != 0){
    if (is.null(alpha)){
      alpha = gsub("%", "", colnames(x)[cols[1]])
      alpha = as.numeric(alpha)*2/100
    }
    else{
      if (length(alpha) > 1){
        warning("\nOnly first alpha will be plotted")
        alpha = alpha[1]
      }
    }
    ll = paste(alpha/2*100, "%", sep="")
    ul = paste((1-alpha/2)*100, "%", sep="")
    if (!(ll %in% colnames(x))){
      stop("\nSpecified alpha not found in the resi object")
    }} else{
      ll = ul = "RESI"
    }

  vnames = rownames(x)[which(rownames(x) != "Residuals")]
  # width of labels
  if ((!(is.null(ycex.axis))) | "cex.axis" %in% names(dots)){
    if (!(is.null(ycex.axis))){
      ycex = ycex.axis
    } else{
      ycex = dots$cex.axis
    }
    w = strwidth(vnames, "inches", cex = ycex)
  } else{
    w = strwidth(vnames, "inches")
    ycex = 1
  }

  # save original margins
  omar = par("mai")
  if (automar){
    # set margins to accomodate labels
    par("mai" = omar + c(0, max(w) - omar[2] + 0.25, 0, 0))
  }

  er = 0.2 * (length(vnames)/(length(vnames) + 1))

  plot(x = x[which(rownames(x) != "Residuals"), "RESI"], y = length(vnames):1,
       xlim = c(0, max(x[,ul], na.rm = TRUE)), xlab = "RESI Estimate",
       yaxt = "n", ylab = "", main = paste("Anova RESI Estimates and ",
                                           (1-alpha)*100, "%", " CIs", sep=""),...)
  if (length(cols) != 0){
    for (i in length(vnames):1){
      lines(x = c(x[-1*(i-length(vnames)) + 1,ll],
                  x[-1*(i-length(vnames)) + 1,ul]), y = c(i,i), ...)
      lines(x = c(x[-1*(i-length(vnames)) + 1,ll],
                  x[-1*(i-length(vnames)) + 1,ll]),
            y = c(i- er,i+er), ...)
      lines(x = c(x[-1*(i-length(vnames)) + 1,ul],
                  x[-1*(i-length(vnames)) + 1,ul]),
            y = c(i-er,i+er), ...)

    }}
  do.call(axis, c(list(side = 2, at = length(vnames):1,
                       labels = vnames, las = 1,
                       cex.axis = ycex), yaxis.args))
  abline(v = 0, lty = 2)
  # return to original margins
  par("mai" = omar)
}
