#' Plot the solution to an L0 segmentation problem
#' @param x output from running estimatedSpikes
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
#'
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
plot.estimatedSpikes <- function(x, xlims = NULL, ...) {
  ind <- 1:length(x$dat)
  rng <- range(c(x$dat, x$fittedValues))
  ylims <- rng #c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlab = "Index")
  } else {
  plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims, xlab = "Time")
  }
  lines(ind, x$fittedValues, col = "blue", lwd = 2)

  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
                                          y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
}

#' Plot simulated data
#' @param x output data from simulateAR1
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
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
plot.simdata <- function(x, xlims = NULL, ...) {
  rng <- range(x$fl)
  ylims <- c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims)
  } else {
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims)
  }
  lines(x$conc, col = "darkgreen", lwd = 2)

  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = spike, x1 = spike, y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "darkgreen", lwd = 1)
  }
}


#' Print simulated data
#' @param x simulated data
#' @param ... arguments to be passed to methods
print.simdata <- function(x, ...){
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")

  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$fl), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Pois mean \t\t", x$poisMean, "\n")
  cat("SD \t\t\t", x$sd, "\n")
}

#' Plot mean squared error vs. tuning parameter from the cross-validation output
#' @param x output from cross validation procedure
#' @param ... arguments to be passed to methods
#' @export
plot.cvSpike <- function(x, ...)
{
  plot(x$lambdas, x$cvError, cex = 0.5, pch = 20, xlab = "Lambda", ylab = "CV error (Average MSE over folds)", log = "xy")
  abline(v = x$lambdaMin, col = "red")
  abline(v = x$lambda1SE, col = "blue")
  legend("topleft", col = c("red", "blue"), c("Min. value", "1se rule"), pch = 20)
}

#' Print estimated spikes
#'
#' @param x estimated spikes
#' @param ... arguments to be passed to methods
#' @export
print.estimatedSpikes <- function(x, ...)
{
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")

  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$dat), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Lambda \t\t\t", x$lambda, "\n")
}

#' Print CV results
#'
#' @param x output from CV
#' @param ... arguments to be passed to methods
#'
#' @export
print.cvSpike <- function(x, ...) {
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Settings: \n")
  cat("Model type \t\t\t", x$type, "\n")
  cat("Largest lambda \t\t\t", max(x$lambdas), "\n")
  cat("Smallest lambda \t\t", min(x$lambdas), "\n")

  cat("\n Output: \n")
  cat("Lambda to min CV Error \t\t", x$lambdaMin, "\n")
  cat("Lambda 1SE for CV Error \t", x$lambda1SE, "\n")
  if (x$optimized) {
    cat("Optimal gamma at Lambda min \t",
        x$optimalGam[x$indexMin, ], "\n")
    cat("Optimal gamma at Lambda 1SE \t",
        x$optimalGam[x$index1SE, ], "\n")
  } else {
    cat("Gamma \t\t\t\t", as.numeric(x$optimalGam), "\n")
  }


}
