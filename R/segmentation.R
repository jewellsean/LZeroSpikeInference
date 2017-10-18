computeCost <- function(dat, optimalFits, ind, n, params) {
  if (n == 1 && optimalFits$hardThreshold == FALSE)
        return(0)

    ## Extract matrix from optimalFits object + AR(1) decay parameter
    sufficientStats <- optimalFits$activeRowSufficientStats
    if (optimalFits$type == "ar1") {
        ## Calculate sumGamma2 = \sum_{t=a}^b \gamma^(2(t-a)) and regression coefficient Ca =
        ## \sum_{t=a}^b y_t \gamma^{t-a} / sumGamma2
        sumGamma2 <- (1 - params ^ (2 * n)) / (1 - params^2)

        Ca <- sufficientStats[ind, 3] / sumGamma2
        if (optimalFits$hardThreshold == T) {
          if (Ca < 0) {
            Ca <- 0
          }
        }

        ## Calculate segment cost \sum_{t=a}^b y_t^2 / 2 - Ca \sum_{t=a}^b y_t \gamma^{t-a} + Ca^2 *
        ## sumGamma2
        segmentCost <- sufficientStats[ind, 2] - Ca * sufficientStats[ind, 3] +
          (Ca ^ 2) * (sumGamma2 / 2)
    }

    if (optimalFits$type == "intercept") {
        sumGammaC2 <- (1 - params^(2 * n))/(1 - params^2)
        sumGammaC <- (1 - params^n)/(1 - (params))

        rescalingFactor <- n * sumGammaC2 - sumGammaC^2

        Ca <- (1/rescalingFactor) * (n * sufficientStats[ind, 3] - sumGammaC * sufficientStats[ind,
            4])
        beta0 <- (1/rescalingFactor) * (sumGammaC2 * sufficientStats[ind, 4] - sumGammaC * sufficientStats[ind,
            3])

        segmentCost <- sufficientStats[ind, 2] - Ca * sufficientStats[ind, 3] - beta0 * sufficientStats[ind,
            4] + (Ca^2) * (sumGammaC2/2) + (beta0^2) * (n/2) + (Ca * beta0) * sumGammaC
    }


    return(segmentCost)
}

updateCostTable <- function(optimalFits, newDatPt, t, params) {

    ## Extract matrix from optimalFits object
    sufficientStats <- optimalFits$activeRowSufficientStats

    if (optimalFits$type == "ar1") {
        ## Sufficient statistics to store in a matrix with columns Inital point a \sum_{t=a}^b (y_t ^ 2)
        ## /2 \sum_{t=a}^b y_t \gamma^{t-a}

        if (is.null(sufficientStats)) {
            optimalFits$activeRowSufficientStats <- matrix(c(0, (newDatPt^2)/2, newDatPt), nrow = 1,
                ncol = 3)
            return(optimalFits)
        }

        sufficientStats[, 2] <- sufficientStats[, 2] + 0.5 * (newDatPt^2)
        sufficientStats[, 3] <- sufficientStats[, 3] + newDatPt * params^(t - (sufficientStats[,
            1] + 1))
    }

    if (optimalFits$type == "intercept") {
        ## Sufficient statistics to store in a matrix with columns Inital point a \sum_{t=a}^b (y_t ^ 2)
        ## /2 \sum_{t=a}^b y_t \gamma^{t-a} \sum_{t=a}^b y_t

        if (is.null(sufficientStats)) {
            optimalFits$activeRowSufficientStats <- matrix(c(0, (newDatPt^2)/2, newDatPt, newDatPt),
                nrow = 1, ncol = 4)
            return(optimalFits)
        }

        sufficientStats[, 2] <- sufficientStats[, 2] + 0.5 * (newDatPt^2)
        sufficientStats[, 3] <- sufficientStats[, 3] + newDatPt * params^(t - (sufficientStats[,
            1] + 1))
        sufficientStats[, 4] <- sufficientStats[, 4] + newDatPt
    }


    optimalFits$activeRowSufficientStats <- sufficientStats
    return(optimalFits)
}

addNewTimeOptimalFits <- function(optimalFits, indicesPtsKeep, t) {
    if (optimalFits$type == "ar1") {
        optimalFits$activeRowSufficientStats <- rbind(optimalFits$activeRowSufficientStats[indicesPtsKeep,
            ], matrix(c(t, 0, 0), nrow = 1, ncol = 3))
    }
    if (optimalFits$type %in% c("dexp", "intercept")) {
        optimalFits$activeRowSufficientStats <- rbind(optimalFits$activeRowSufficientStats[indicesPtsKeep,
            ], matrix(c(t, 0, 0, 0), nrow = 1, ncol = 4))
    }
    return(optimalFits)
}

computeSegmentation <- function(dat, params, penalty, type, hardThreshold) {
    n <- length(dat)
    ## rows index 0...t, cols index: t, F(t), \tau'_t, # taus
    table <- matrix(0, nrow = n + 1, ncol = 4)
    table[1, 2] <- -penalty
    R = c(0)  ## restricted set for mins

    optimalFits <- list(type = type, activeRowSufficientStats = NULL, hardThreshold = hardThreshold)
    keepChgPts <- list()
    for (t in 1:n) {
        nR <- length(R)
        minimizers <- numeric(nR)  ## minimize over R
        ## update cost table for this t
        optimalFits <- updateCostTable(optimalFits, dat[t], t, params)

        for (ind in 1:nR) {
            minimizers[ind] <- table[R[ind] + 1, 2] +
              computeCost(dat[(R[ind] + 1):t], optimalFits, ind, t - R[ind], params)
        }
        table[t + 1, 1] <- t
        table[t + 1, 2] <- min(minimizers) + penalty
        table[t + 1, 3] <- R[which.min(minimizers)]

        indicesPtsKeep <- (minimizers < table[t + 1, 2])
        # print(paste0("Cost F(r): ", table[t + 1, 2]))
        # print(minimizers)
        R <- c(R[indicesPtsKeep], t)
        keepChgPts[[t]] <- R

        table[t + 1, 4] <- length(R)

        if (t < n)
            optimalFits <- addNewTimeOptimalFits(optimalFits, indicesPtsKeep, t)
    }
    return(list(table = table, keepChgPts = keepChgPts))
}

findChangePts <- function(vecChgPts) {
    ind <- length(vecChgPts)
    changePts <- c()
    while (ind > 1) {
        ind <- vecChgPts[ind] + 1
        changePts <- c(ind, changePts)
    }
    changePts <- changePts - 1

    return(changePts)
}

computeFittedValues <- function(dat, changePts, params, type, hardThreshold = FALSE) {
    n <- length(dat)
    nSegments <- length(changePts)
    changePts <- c(changePts, n)
    if (type == "ar1") {
        X <- matrix(0, nrow = n, ncol = nSegments)
        for (k in 1:nSegments) {
            X[(changePts[k] + 1):changePts[k + 1], k] <- params^(0:(changePts[k + 1] - (changePts[k] +
                1)))
        }
        fit <- lm(dat ~ X - 1)

        if (hardThreshold == T) {
          return(pmax(fit$fitted.values, 0))
        } else {

          if (sum(fit$coefficients < 0) > 0) {
          warning("Check model fit carefully. In some segments calcium concentration may not 'decay' as expected. Most observed datapoints should be positive.")
        }
          return(fit$fitted.values)
        }


    }

    if (type == "intercept") {
        X <- matrix(0, nrow = n, ncol = 2 * nSegments)
        ind <- c(1, 2)
        for (k in 1:nSegments) {
            X[(changePts[k] + 1):changePts[k + 1], ind] <- cbind(params^(0:(changePts[k + 1] - (changePts[k] +
                1))), rep(1, changePts[k + 1] - changePts[k]))
            ind <- c(ind[2] + 1, ind[2] + 2)
        }

        fit <- lm(dat ~ X - 1)
        return(fit$fitted.values)
    }

}

#' Estimate spike train, underlying calcium concentration, and changepoints based on fluorescence
#' trace.
#'
#' @param dat fluorescence data
#' @param gam a scalar value for the AR(1)/AR(1) + intercept decay parameter.
#' @param lambda tuning parameter lambda
#' @param type type of model, must be one of AR(1) 'ar1', AR(1) + intercept 'intercept'.
#' @param calcFittedValues TRUE to calculate fitted values.
#' @param hardThreshold boolean specifying whether the calcium concentration must be non-negative (in the AR-1 problem)
#'
#' @return Returns a list with elements:
#' @return \code{changePts} the set of changepoints
#' @return \code{spikes} the set of spikes
#' @return \code{fittedValues} estimated calcium concentration
#'
#' @details
#'
#' This algorithm solves the optimization problems
#'
#' \strong{AR(1)-model:}
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1} }
#' for the global optimum, where $y_t$ is the observed fluorescence at the tth
#' timepoint.
#'
#' If hardThreshold = T then the additional constraint
#' c_t >= 0 is added to the optimzation problem above.
#'
#' \strong{AR(1) with intercept:}
#' minimize_{c1,...,cT,b1,...,bT} 0.5 sum_{t=1}^T (y_t - c_t - b_t)^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1}, b_t neq b_{t-1} }
#' where the indicator variable 1_{(A,B)} equals 1 if the event A cup B holds, and equals zero otherwise.
#'
#' See Jewell and Witten (2017) <arXiv:1703.08644>, section 2 and 5.
#'
#' Note that "changePts" and "spikes" differ by one index due to a quirk of the conventions used in the changepoint literature and the definition of a spike.
#'
#'
#' @examples
#'
#' sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' # AR(1) model
#'
#' fit <- estimateSpikes(sim$fl, gam = 0.998, lambda = 1, type = "ar1")
#' plot(fit)
#' print(fit)
#'
#' # AR(1) + intercept model
#' fit <- estimateSpikes(sim$fl, gam = 0.998, lambda = 1, type = "intercept")
#' plot(fit)
#' print(fit)
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
#'
#' @export
estimateSpikes <- function(dat, gam, lambda,
                           type = "ar1", calcFittedValues = TRUE, hardThreshold = FALSE) {
  checkValidType(type)
  checkValidParameters(gam, type)
  checkData(dat)
  table <- computeSegmentation(dat, gam, lambda, type, hardThreshold)
  keepChgPts <- table$keepChgPts
  table <- table$table
  changePts <- findChangePts(table[, 3])
  spikes <- changePts[-1] + 1
  if (calcFittedValues) {
    fittedValues <- computeFittedValues(dat, changePts, gam, type, hardThreshold)
  } else {
    fittedValues <- NULL
  }
  out <- list(spikes = spikes, fittedValues = fittedValues,
              dat = dat, type = type, changePts = changePts,
              call = match.call(),
              gam = gam,
              lambda = lambda,
              cost = table[(2:dim(table)[[1]]), 2],
              nIntervals = table[(2:dim(table)[[1]]) ,4],
              keepChgPts = keepChgPts)
  class(out) <- "estimatedSpikes"
  return(out)
}

checkValidParameters <- function(params, type)
{
  if (type %in% c("ar1", "intercept")) {
    if(params >= 1 || params <= 0)
      stop("Decay parameter must satisfy 0 < gamma < 1")
  }
}

checkValidType <- function(type) {
  if (!(type %in% c("ar1", "intercept")))
    stop("Model not implemented. Type must be one of ar1 or intercept.")
}

checkData <- function(dat) {
  if(mean(dat < 0) > 0.5) warning("Most observed data points should be positive.")
}
