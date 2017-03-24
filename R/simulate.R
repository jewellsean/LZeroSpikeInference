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
#' \code{\link{simulateDexp}},
#' \code{\link{plot.simdata}}.
#'
#' @examples
#' sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
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

  spikesOut <- unique(eta * (1:n))
  out <- list(spikes = spikesOut[-1], fl = f, conc = c, call = match.call(),
              gam = gam, poisMean = poisMean, type = "ar1",
              sd = sd, seed = seed)
  class(out) <- "simdata"
  return(out)
}


#' Simulate fluorescence trace based on difference of exponential
#' calcium dynamics
#'
#' @details
#' Simulate fluorescence trace based on difference of exponential
#' calcium dynamics
#'
#' y_t =  r_t + eps, eps ~ N(0, sd)
#'
#' c_t = gammaC * c_{t-1} + s_t
#'
#' d_t = gammaD * d_{t-1} + s_t
#'
#' r_t = (c_t - d_t)
#'
#' s_t ~ Pois(poisMean)
#'
#' @return list of spikes, fluorescence, and calcium concentration
#'
#' @param n number of timesteps
#' @param seed random seed
#' @param gams vector (gammaC, gammaD)
#' @param poisMean mean for Poisson distributed spikes
#' @param sd standard deviation
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
#' \code{\link{simulateDexp}},
#' \code{\link{plot.simdata}}.
#' @examples
#' sim <- simulateDexp(n = 500, gams = c(0.998, 0.7), poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#' @export
simulateDexp <- function(n, gams, poisMean, sd, seed)
{
  set.seed(seed)
  eta <- numeric(n)
  c <- numeric(n)
  d <- numeric(n)
  r <- numeric(n)
  f <- numeric(n)
  for (i in 1:n)
  {
    eta[i] <- rpois(1, poisMean)
    if (i > 1) {
        c[i] <- gams[1] * c[i - 1] + eta[i]
        d[i] <- gams[2] * d[i - 1] + eta[i]
    } else {
      c[i] <- 0
      d[i] <- 0
    }

    r[i] <- (c[i] - d[i])

    f[i] <- r[i] + rnorm(n = 1, mean = 0, sd = sd)
  }

  spikesOut <- unique(eta * (1:n))
  out <- list(spikes = spikesOut[-1], fl = f, conc = r, call = match.call(),
              gam = gams, poisMean = poisMean, type = "dexp",
              sd = sd, seed = seed)
  class(out) <- "simdata"
  return(out)
}
