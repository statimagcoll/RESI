#' @importFrom graphics abline axis lines
#' @export
plot.resi <- function(x, ycex.axis = 1, alpha = NULL, ...){
  if (is.null(alpha)){
    alpha = x$alpha[1]
  }
  else{
    if (length(alpha) > 1){
      warning('\nOnly first alpha will be plotted')
      alpha = alpha[1]
    }
    if (!(alpha %in% x$alpha)){
      if (is.null(x$boot.results)){
        stop('\nSpecified alpha not found in the resi object')
      }
      CIs = apply(x$boot.results[,2:(1+nrow(x$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      CIs = t(CIs)
      x$coefficients[1:nrow(CIs), c(paste(alpha/2*100, '%', sep=''), paste((1-rev(alpha)/2)*100, '%', sep=''))] = CIs
    }
  }
  ll = paste(alpha/2*100, "%", sep = "")
  ul = paste((1-alpha/2)*100, "%", sep = "")
  plot(x = x$coefficients[,"RESI"], y = length(rownames(x$coefficients)):1,
       xlim = c(min(0, min(x$coefficients[,ll])), max(x$coefficients[,ul])),
       xlab = "RESI Estimate", yaxt = "n", ylab = "",
       main = paste("Coefficient RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""),...)
  for (i in length(rownames(x$coefficients)):1){
    lines(x = c(x$coefficients[-1*(i-length(rownames(x$coefficients))) + 1,ll],
                x$coefficients[-1*(i-length(rownames(x$coefficients))) + 1,ul]), y = c(i,i))
  }
  axis(2, length(rownames(x$coefficients)):1, rownames(x$coefficients), las = 1, cex.axis = ycex.axis)
  abline(v = 0, lty = 2)
}

#' @export
plot.summary_resi <- plot.resi

#' @export
plot.anova_resi <- function(x, alpha = NULL, ycex.axis = 1,...){
  cols = grep('%', colnames(x))
  if (is.null(alpha)){
    alpha = gsub("%", "", colnames(x)[cols[1]])
    alpha = as.numeric(alpha)*2/100
  }
  else{
    if (length(alpha) > 1){
      warning('\nOnly first alpha will be plotted')
      alpha = alpha[1]
    }
  }
  ll = paste(alpha/2*100, "%", sep="")
  ul = paste((1-alpha/2)*100, "%", sep="")
  if (!(ll %in% colnames(x))){
    stop('\nSpecified alpha not found in the resi object')
  }

  plot(x = x[which(rownames(x) != "Residuals"), "RESI"], y = length(which(rownames(x) != "Residuals")):1,
       xlim = c(0, max(x[,ul], na.rm = TRUE)), xlab = "RESI Estimate",
       yaxt = "n", ylab = "", main = paste("Anova RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""),...)
  for (i in length(which(rownames(x) != "Residuals")):1){
    lines(x = c(x[-1*(i-length(which(rownames(x) != "Residuals"))) + 1,ll],
                x[-1*(i-length(which(rownames(x) != "Residuals"))) + 1,ul]), y = c(i,i))
  }
  axis(2, length(which(rownames(x) != "Residuals")):1,
       rownames(x)[which(rownames(x) != "Residuals")], las = 1, cex.axis = ycex.axis)
  abline(v = 0, lty = 2)
}
