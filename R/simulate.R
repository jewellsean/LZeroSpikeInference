#' Simulate flourscence trace based on generative model in paper
#' @param n number of timesteps
#' @param seed random seed
#' @param poisRate Poisson mean spike
#' @param decay ar(1) decay rate
#' @param sd standard deviation for white noise process
#'
#' @return list of spikes, flourscence, and concentration
#'
#' @examples
#' simulateAR1(n = 100, seed = 1, poisRate = 0.05, decay = 0.998, sd = 0.05)
#'
#' @export
simulateAR1 <- function(n, seed, poisRate, decay, sd)
{
  set.seed(seed)
  eta <- numeric(n)
  c <- numeric(n)
  f <- numeric(n)
  for (i in 1:n)
  {
    eta[i] <- rpois(1, poisRate)
    if (i > 1)
      c[i] <- decay * c[i - 1] + eta[i]
    else
      c[i] <- 1 + eta[i]

    f[i] <- c[i] + rnorm(n = 1, mean = 0, sd = sd)
  }

  spikesOut <- unique(eta * (1:n))
  return(list(spikes = spikesOut[-1], fl = f, conc = c))
}
