#' Plot the solution to an L0 segmentation problem
#' @param lzeroFit output from running estimatedSpikes
#' @param xlims optional parameter to specify the x-axis limits
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{estimateSpikes}},
#' \code{\link{print.estimatedSpikes}},
#' \code{\link{plot.estimatedSpikes}}.
#'
#' \strong{Cross validation:}
#' \code{\link{cv.estimateSpikes}},
#' \code{\link{print.cvSpike}},
#' \code{\link{plot.cvSpike}}.
#'
#' \strong{Simulation:}
#' \code{\link{simulateAR1}},
#' \code{\link{plot.simdata}}.
#' @export
#'
plot.estimatedSpikes <- function(lzeroFit, xlims = NULL) {
  ind <- 1:length(lzeroFit$dat)
  rng <- range(c(lzeroFit$dat, lzeroFit$fittedValues))
  ylims <- rng #c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(ind, lzeroFit$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlab = "Index")
  } else {
  plot(ind, lzeroFit$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims, xlab = "Time")
  }
  lines(ind, lzeroFit$fittedValues, col = "blue", lwd = 2)

  hh <- 0.01 * diff(ylims)
  for (spike in lzeroFit$spikes)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
                                          y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
}

#' Plot simulated data
#' @param dat output data from simulateAR1 or simulateDexp
#' @param xlims optional parameter to specify the x-axis limits
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{estimateSpikes}},
#' \code{\link{print.estimatedSpikes}},
#' \code{\link{plot.estimatedSpikes}}.
#'
#' \strong{Cross validation:}
#' \code{\link{cv.estimateSpikes}},
#' \code{\link{print.cvSpike}},
#' \code{\link{plot.cvSpike}}.
#'
#' \strong{Simulation:}
#' \code{\link{simulateAR1}},
#' \code{\link{plot.simdata}}.
#' @return Plot with simulated fluorescence (dark grey circles), calcium concentration (dark green line) and spikes (dark green tick marks on x-axis)
#'
#' @examples
#'
#' sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' @export
#'
plot.simdata <- function(dat, xlims = NULL) {
  rng <- range(dat$fl)
  ylims <- c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(dat$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims)
  } else {
    plot(dat$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims)
  }
  lines(dat$conc, col = "darkgreen", lwd = 2)

  hh <- 0.01 * diff(ylims)
  for (spike in dat$spikes)
  {
    segments(x0 = spike, x1 = spike, y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "darkgreen", lwd = 1)
  }
}


#' Print simulated data
#' @param dat simulated data
print.simdata <- function(dat){
  cat("\n Call: \n")
  dput(dat$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(dat$spikes), "\n")

  cat("\n Settings: \n")
  cat("Data length \t\t", length(dat$fl), "\n")
  cat("Model type \t\t", dat$type, "\n")
  cat("Gamma \t\t\t", dat$gam, "\n")
  cat("Pois mean \t\t", dat$poisMean, "\n")
  cat("SD \t\t\t", dat$sd, "\n")
}

#' Plot mean squared error vs. tuning parameter from the cross-validation output
#' @param cvOut output from cross validation procedure
#'
#' @export
plot.cvSpike <- function(cvOut)
{
  plot(cvOut$lambdas, cvOut$cvError, cex = 0.5, pch = 20, xlab = "Lambda", ylab = "CV error (Average MSE over folds)", log = "xy")
  abline(v = cvOut$lambdaMin, col = "red")
  abline(v = cvOut$lambda1SE, col = "blue")
  legend("topleft", col = c("red", "blue"), c("Min. value", "1se rule"), pch = 20)
}

#' Print estimated spikes
#'
#' @param lzeroFit estimated spikes
#'
#' @export
print.estimatedSpikes <- function(lzeroFit)
{
  cat("\n Call: \n")
  dput(lzeroFit$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(lzeroFit$spikes), "\n")

  cat("\n Settings: \n")
  cat("Data length \t\t", length(lzeroFit$dat), "\n")
  cat("Model type \t\t", lzeroFit$type, "\n")
  cat("Gamma \t\t\t", lzeroFit$gam, "\n")
  cat("Lambda \t\t\t", lzeroFit$lambda, "\n")
}

#' Print CV results
#'
#' @param cvOut output from CV
#'
#' @export
print.cvSpike <- function(cvOut) {
  cat("\n Call: \n")
  dput(cvOut$call)
  cat("\n Settings: \n")
  cat("Model type \t\t\t", cvOut$type, "\n")
  cat("Largest lambda \t\t\t", max(cvOut$lambdas), "\n")
  cat("Smallest lambda \t\t", min(cvOut$lambdas), "\n")

  cat("\n Output: \n")
  cat("Lambda to min CV Error \t\t", cvOut$lambdaMin, "\n")
  cat("Lambda 1SE for CV Error \t", cvOut$lambda1SE, "\n")
  if (cvOut$optimized) {
    cat("Optimal gamma at Lambda min \t",
        cvOut$optimalGam[cvOut$indexMin, ], "\n")
    cat("Optimal gamma at Lambda 1SE \t",
        cvOut$optimalGam[cvOut$index1SE, ], "\n")
  } else {
    cat("Gamma \t\t\t\t", as.numeric(cvOut$optimalGam), "\n")
  }


}
