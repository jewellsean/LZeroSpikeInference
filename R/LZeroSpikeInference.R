#' LZeroSpikeInference: LZeroSpikeInference: A package for estimating spike
#' times from calcium imaging data using an L0 penalty
#'
#' This package implements an algorithm for deconvolving calcium imaging data
#' for a single neuron in order to estimate the times at which the neuron
#' spikes.
#'
#' @section LZeroSpikeInference functions:
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
#'
#' @details
#'
#' This package implements an algorithm for deconvolving calcium imaging data
#' for a single neuron in order to estimate the times at which the neuron
#' spikes. This algorithm solves the optimization problems
#'
#' \strong{AR(1)-model:}
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1} }
#' for the global optimum, where $y_t$ is the observed fluorescence at the tth
#' timepoint.
#'
#' \strong{AR(1) with intercept:}
#' minimize_{c1,...,cT,b1,...,bT} 0.5 sum_{t=1}^T (y_t - c_t - b_t)^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1}, b_t neq b_{t-1} }
#' where the indicator variable 1_{(A,B)} equals 1 if the event A cup B holds, and equals zero otherwise.
#'
#' \strong{Difference of Exponentials:}
#' minimize_{c1,...,cT,d1,...,dT} 0.5 sum_{t=1}^{T} ( y_t - (c_t - d_t) )^2 + lambda sum_{t =2}^{T} 1_{ c_t neq gamma _c c_{t-1}, d_t neq gamma_d d_{t-1} }.
#'
#' See Jewell and Witten, Exact Spike Train Inference Via L0 Optimization (2017)
#'
#' @examples
#' # To reproduce Figure 1 of Jewell and Witten, Exact Spike Train Inference
#' # Via L0 Optimization (2017)
#'
#' sampleData <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.15, seed = 8)
#' fit <- estimateSpikes(sampleData$fl, gam = 0.998, lambda = 8, type = "ar1")
#' plot(fit)
#'
#' @docType package
#' @name LZeroSpikeInference
NULL
