#' @importFrom graphics abline axis lines
#' @export
plot.resi = function(x, alpha = NULL, ycex.axis = NULL, yaxis.args = list(),
                     automar = TRUE, ...){
  dots = list(...)

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
      CIs = apply(x$boot.results[,2:(1+nrow(x$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      CIs = t(CIs)
      x$coefficients[1:nrow(CIs), c(paste(alpha/2*100, "%", sep=""), paste((1-rev(alpha)/2)*100, "%", sep=""))] = CIs
    }
  }
  ll = paste(alpha/2*100, "%", sep = "")
  ul = paste((1-alpha/2)*100, "%", sep = "")

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
          ...)

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
  for (i in length(vnames):1){
    lines(x = c(x[-1*(i-length(vnames)) + 1,ll],
                x[-1*(i-length(vnames)) + 1,ul]), y = c(i,i), ...)
    lines(x = c(x[-1*(i-length(vnames)) + 1,ll],
                x[-1*(i-length(vnames)) + 1,ll]),
          y = c(i- er,i+er), ...)
    lines(x = c(x[-1*(i-length(vnames)) + 1,ul],
                x[-1*(i-length(vnames)) + 1,ul]),
          y = c(i-er,i+er), ...)

  }
  do.call(axis, c(list(side = 2, at = length(vnames):1,
                       labels = vnames, las = 1,
                       cex.axis = ycex), yaxis.args))
  abline(v = 0, lty = 2)
  # return to original margins
  par("mai" = omar)
}
