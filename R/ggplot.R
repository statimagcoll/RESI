#' @import ggplot2
#' @export
ggplot.resi = function(x, alpha = NULL, error.bars = TRUE, ...){
  #browser()
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

  dat = x$coefficients
  dat$ggpos = length(vnames):1
  dat = dat[order(dat$ggpos),]
  p = ggplot(dat, aes(x = RESI, y = ggpos)) +
    geom_point() +
    sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                      aes(x = get(ll),
                                                          xend = get(ul),
                                                          y = ggpos, yend = ggpos,
                                                         ))) +
    scale_y_continuous(name = NULL, breaks = 1:length(vnames), labels = rev(vnames)) +
    geom_vline(xintercept = 0, linetype="dotted") +
    xlab("RESI Estimate") +
    ggtitle(paste("Coefficient RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""))
  if (error.bars) {
    er = 0.2 * (length(vnames)/(length(vnames)))
    p = p + sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                              aes(x = get(ll),
                                                                  xend = get(ll),
                                                                  y = ggpos - er, yend = ggpos + er,
                                                              ))) +
      sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                          aes(x = get(ul),
                                                              xend = get(ul),
                                                              y = ggpos - er, yend = ggpos + er,
                                                          )))
  }
  p
}

#' @export
ggplot.summary_resi = ggplot.resi

#' @export
ggplot.anova_resi = function(x, alpha = NULL, error.bars = TRUE, ...){
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

  dat = as.data.frame(x[which(rownames(x) != "Residuals"),])
  dat$ggpos = length(vnames):1
  dat = dat[order(dat$ggpos),]
  p = ggplot(dat, aes(x = RESI, y = ggpos)) +
    geom_point() +
    sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                      aes(x = get(ll),
                                                          xend = get(ul),
                                                          y = ggpos, yend = ggpos,
                                                      ))) +
    scale_y_continuous(name = NULL, breaks = 1:length(vnames), labels = rev(vnames)) +
    xlab("RESI Estimate") +
    geom_vline(xintercept = 0, linetype="dotted") +
    ggtitle(paste("ANOVA RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep=""))
  if (error.bars) {
    er = 0.2 * (length(vnames)/(length(vnames)))
    p = p + sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                              aes(x = get(ll),
                                                                  xend = get(ll),
                                                                  y = ggpos - er, yend = ggpos + er,
                                                              ))) +
      sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                        aes(x = get(ul),
                                                            xend = get(ul),
                                                            y = ggpos - er, yend = ggpos + er,
                                                        )))
  }
  p
}





