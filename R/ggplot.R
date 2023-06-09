#' Plotting RESI Estimates and CIs
#'
#' This function uses ggplot2 graphics to plot robust effect size (RESI) estimates and confidence intervals
#' from `resi`, `summary_resi`, and `anova_resi` objects.
#' @param data Object of `resi`, `summary_resi`, or `anova_resi` class
#' @param mapping Ignored, included for consistency with `ggplot` generic
#' @param alpha Numeric, desired alpha level for confidence intervals
#' @param error.bars Logical, whether to include end caps on the confidence intervals. Default = `TRUE`
#' @param ... Ignored
#' @param environment Ignored, included for consistency with `ggplot` generic
#' @return Returns a ggplot of RESI point estimates
#' @examples
#' # create a resi object
#' resi_obj <- resi(lm(charges ~ region * age + bmi + sex, data = RESI::insurance),
#' nboot = 10)
#'
#' # plot ANOVA table
#' ggplot2::ggplot(anova(resi_obj))
#' @import ggplot2
#' @export
ggplot.resi = function(data, mapping, alpha = NULL, error.bars = TRUE, ...,
                       environment){
  dots = list(...)

  if (!is.null(data$alpha)){
    if (is.null(alpha)){
      alpha = data$alpha[1]
    }
    else{
      if (length(alpha) > 1){
        warning("\nOnly first alpha will be plotted")
        alpha = alpha[1]
      }
      if (!(alpha %in% data$alpha)){
        if (is.null(data$boot.results)){
          stop("\nSpecified alpha not found in the resi object")
        }
        CIs = apply(data$boot.results$t[,2:(1+nrow(data$coefficients))], 2,  quantile, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
        CIs = t(CIs)
        data$coefficients[1:nrow(CIs), c(paste(alpha/2*100, "%", sep=""), paste((1-rev(alpha)/2)*100, "%", sep=""))] = CIs
      }
    }
    ll = paste(alpha/2*100, "%", sep = "")
    ul = paste((1-alpha/2)*100, "%", sep = "")} else{
      ll = ul = "RESI"
    }

  vnames = rownames(data$coefficients)

  dat = data$coefficients
  ggpos = length(vnames):1
  dat$ggpos = ggpos
  dat = dat[order(dat$ggpos),]
  RESI = dat$RESI
  p = ggplot(dat, aes(data = RESI, y = ggpos)) +
    geom_point() +
    scale_y_continuous(name = NULL, breaks = 1:length(vnames), labels = rev(vnames)) +
    geom_vline(xintercept = 0, linetype="dotted") +
    xlab("RESI Estimate") +
    ggtitle(ifelse(is.null(data$alpha), "Coefficient RESI Estimates",
                   paste("Coefficient RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep="")))
  if (!is.null(data$alpha)){
    p = p + sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                              aes(x = get(ll),
                                                                  xend = get(ul),
                                                                  y = ggpos, yend = ggpos,
                                                              )))
    if (error.bars) {
      er = 0.2 * (length(vnames)/(length(vnames) + 1))
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
    }}
  p
}

#' @export
ggplot.summary_resi = ggplot.resi

#' @export
ggplot.anova_resi = function(data, mapping, alpha = NULL, error.bars = TRUE, ...,
                             environment){
  dots = list(...)
  cols = grep("%", colnames(data))
  if (length(cols) != 0){
    if (is.null(alpha)){
      alpha = gsub("%", "", colnames(data)[cols[1]])
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
    if (!(ll %in% colnames(data))){
      stop("\nSpecified alpha not found in the resi object")
    }} else{
      ll = ul = "RESI"
    }

  vnames = rownames(data)[which(rownames(data) != "Residuals")]

  dat = as.data.frame(data[which(rownames(data) != "Residuals"),])
  ggpos = length(vnames):1
  dat$ggpos = ggpos
  dat = dat[order(dat$ggpos),]
  RESI = dat$RESI
  p = ggplot(dat, aes(x = RESI, y = ggpos)) +
    geom_point()  +
    scale_y_continuous(name = NULL, breaks = 1:length(vnames), labels = rev(vnames)) +
    xlab("RESI Estimate") +
    geom_vline(xintercept = 0, linetype="dotted") +
    ggtitle(ifelse(length(cols) == 0, "ANOVA RESI Estimates",
                   paste("ANOVA RESI Estimates and ", (1-alpha)*100, "%", " CIs", sep="")))

  if (length(cols) != 0){
    p = p + sapply(length(vnames):1, function(i) geom_segment(data = dat[-1*(i-length(vnames)) + 1,],
                                                      aes(x = get(ll),
                                                          xend = get(ul),
                                                          y = ggpos, yend = ggpos,
                                                      )))
    if (error.bars) {
      er = 0.2 * (length(vnames)/(length(vnames) + 1))
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
    }}
  p
}





