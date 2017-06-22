#' Simulate fluorescence trace based on simple AR(1) generative model
#'
#' @details
#' Simulate fluorescence trace based on simple AR(1) generative model
#'
#' y_t = c_t + eps, eps ~ N(0, sd)
#'
#' c_t = gam * c_{t-1} + s_t
#'
#' s_t ~ Pois(poisMean)
#'
#' @param n number of timesteps
#' @param seed random seed
#' @param gam AR(1) decay rate
#' @param poisMean mean for Poisson distributed spikes
#' @param sd standard deviation
#'
#' @return spikes, fluorescence, and calcium concentration
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
#'
#' @examples
#' sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#' @import stats
#' @export
simulateAR1 <- function(n, gam, poisMean, sd, seed)
{
  set.seed(seed)
  eta <- numeric(n)
  c <- numeric(n)
  f <- numeric(n)
  for (i in 1:n)
  {
    eta[i] <- rpois(1, poisMean)
    if (i > 1)
      c[i] <- gam * c[i - 1] + eta[i]
    else
      c[i] <- 1 + eta[i]

    f[i] <- c[i] + rnorm(n = 1, mean = 0, sd = sd)
  }

  spikesOut <- unique((eta > 0) * (1:n))
  out <- list(spikes = spikesOut[-1], fl = f, conc = c, call = match.call(),
              gam = gam, poisMean = poisMean, type = "ar1",
              sd = sd, seed = seed)
  class(out) <- "simdata"
  return(out)
}
