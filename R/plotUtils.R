#' Plot the solution to an L0 segmentation problem
#' @param lzeroFit output from running a segmentation
#' @param trueSpikes if true spikes are known and provided this adds tick marks at specified index locations
#' @param xTimes timeseries corresponding to the time that each index of the fluorescence data correspond to. Include to plot the x-axis on this scale
#' @param spikeTimes times, on the same scale as xTimes, that spikes occur in the true (ground truth) data
#'
#' @export
#'
plotSegmentation <- function(lzeroFit, trueSpikes = NULL, xTimes = NULL, spikeTimes = NULL) {

  if (is.null(xTimes)) {
    xlab <- "Index"
    ind <- 1:length(lzeroFit$dat)
  } else {
    ind <- xTimes
    xlab <- "Time"
  }

  rng <- range(lzeroFit$dat)
  ylims <- c(floor(rng[1]), ceiling(rng[2]))
  plot(ind, lzeroFit$dat, cex = 0.5, pch = 20, col = "darkgrey", ylim = ylims, ylab = "", xlab = xlab)
  lines(ind, lzeroFit$fittedValues, col = "blue", lwd = 2)

  for (spike in lzeroFit$spikes) abline(v = ind[spike], col = "blue", lwd = 2)

  if (!is.null(trueSpikes)) {
    hh <- 0.05 * diff(ylims)
    for (spike in trueSpikes) segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
                                       y1 = hh + ylims[1], col = "black", lwd = 2)
  }

  if (!is.null(spikeTimes)) {
    hh <- 0.05 * diff(ylims)
    for (spike in spikeTimes) segments(x0 = spike, x1 = spike, y0 = ylims[1] - hh, y1 = hh +
                                         ylims[1], col = "black", lwd = 2)
  }
}

#' Plot mean squared error vs. tuning parameter from the cross-validation output
#' @param cvOut output list from cross validation procedure
#' @export
plotCV <- function(cvOut)
{
  plot(cvOut$lambdas, cvOut$cvError, cex = 0.5, pch = 20, xlab = "Penalty", ylab = "CV error", log = "xy")
  abline(v = cvOut$lambdaMin, col = "red")
  abline(v = cvOut$lambda1SE, col = "blue")
  legend("bottomright", col = c("red", "blue"), c("Min. value", "1se rule"), pch = 20)
}
