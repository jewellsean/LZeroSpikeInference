computeCost <- function(dat, optimalFits, ind, n, params) {
    if (n == 1)
        return(0)

    ## Extract matrix from optimalFits object + AR(1) decay parameter
    sufficientStats <- optimalFits$activeRowSufficientStats
    if (optimalFits$type == "ar1") {
        ## Calculate sumGamma2 = \sum_{t=a}^b \gamma^(2(t-a)) and regression coefficient Ca =
        ## \sum_{t=a}^b y_t \gamma^{t-a} / sumGamma2
        sumGamma2 <- (1 - params^(2 * n))/(1 - params^2)
        Ca <- sufficientStats[ind, 3]/sumGamma2

        ## Calculate segment cost \sum_{t=a}^b y_t^2 / 2 - Ca \sum_{t=a}^b y_t \gamma^{t-a} + Ca^2 *
        ## sumGamma2
        segmentCost <- sufficientStats[ind, 2] - Ca * sufficientStats[ind, 3] + (Ca^2) * (sumGamma2/2)
    }

    if (optimalFits$type == "dexp") {
        sumGammaC2 <- (1 - params$gammaC^(2 * n))/(1 - params$gammaC^2)
        sumGammaD2 <- (1 - params$gammaD^(2 * n))/(1 - params$gammaD^2)
        sumGammaCD <- (1 - (params$gammaD * params$gammaC)^n)/(1 - (params$gammaD * params$gammaC))

        rescalingFactor <- (sumGammaC2 * sumGammaD2 - sumGammaCD^2)

        Ca <- (1/rescalingFactor) * (sumGammaD2 * sufficientStats[ind, 3] - sumGammaCD * sufficientStats[ind,
            4])
        Da <- (1/rescalingFactor) * (-sumGammaC2 * sufficientStats[ind, 4] + sumGammaCD * sufficientStats[ind,
            3])

        segmentCost <- sufficientStats[ind, 2] - Ca * sufficientStats[ind, 3] + Da * sufficientStats[ind,
            4] + (Ca^2) * (sumGammaC2/2) + (Da^2) * (sumGammaD2/2) - (Ca * Da) * sumGammaCD
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
    if (optimalFits$type == "dexp") {
        ## Sufficient statistics to store in a matrix with columns Inital point a \sum_{t=a}^b (y_t ^ 2)
        ## /2 \sum_{t=a}^b y_t \gamma_C^{t-a} \sum_{t=a}^b y_t \gamma_D^{t-a}

        if (is.null(sufficientStats)) {
            optimalFits$activeRowSufficientStats <- matrix(c(0, (newDatPt^2)/2, newDatPt, newDatPt),
                nrow = 1, ncol = 4)
            return(optimalFits)
        }

        sufficientStats[, 2] <- sufficientStats[, 2] + 0.5 * (newDatPt^2)
        sufficientStats[, 3] <- sufficientStats[, 3] + newDatPt * params$gammaC^(t - (sufficientStats[,
            1] + 1))
        sufficientStats[, 4] <- sufficientStats[, 4] + newDatPt * params$gammaD^(t - (sufficientStats[,
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

computeSegmentation <- function(dat, params, penalty, type) {
    n <- length(dat)
    table <- matrix(0, nrow = n + 1, ncol = 3)  ## rows index 0...t, cols index: t, F(t), \tau'_t
    table[1, 2] <- -penalty
    R = c(0)  ## restricted set for mins

    optimalFits <- list(type = type, activeRowSufficientStats = NULL)

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
        R <- c(R[indicesPtsKeep], t)

        if (t < n)
            optimalFits <- addNewTimeOptimalFits(optimalFits, indicesPtsKeep, t)
    }
    return(table)
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

computeFittedValues <- function(dat, changePts, params, type) {
    n <- length(dat)
    nSegments <- length(changePts)
    changePts <- c(changePts, n)
    if (type == "ar1") {
        X <- matrix(0, nrow = n, ncol = nSegments)
        for (k in 1:nSegments) {
            X[(changePts[k] + 1):changePts[k + 1], k] <- params^(0:(changePts[k + 1] - (changePts[k] +
                1)))
        }
        return(lm(dat ~ X - 1)$fitted.values)
    }
    if (type == "dexp") {
        X <- matrix(0, nrow = n, ncol = 2 * nSegments)
        ind <- c(1, 2)
        for (k in 1:nSegments) {
            X[(changePts[k] + 1):changePts[k + 1], ind] <- cbind(params$gammaC^(0:(changePts[k +
                1] - (changePts[k] + 1))), -params$gammaD^(0:(changePts[k + 1] - (changePts[k] +
                1))))
            ind <- c(ind[2] + 1, ind[2] + 2)
        }
        return(lm(dat ~ X - 1)$fitted.values)
    }

    if (type == "intercept") {
        X <- matrix(0, nrow = n, ncol = 2 * nSegments)
        ind <- c(1, 2)
        for (k in 1:nSegments) {
            X[(changePts[k] + 1):changePts[k + 1], ind] <- cbind(params^(0:(changePts[k + 1] - (changePts[k] +
                1))), rep(1, changePts[k + 1] - changePts[k]))
            ind <- c(ind[2] + 1, ind[2] + 2)
        }
        return(lm(dat ~ X - 1)$fitted.values)
    }

}

#' Segments based on type of model and penalty
#' @param dat noisy flourscence data
#' @param params model parameters. For the AR(1) and AR(1) intercept models this is the scalar decay parameter; this is a
#' dataframe with two parameters gammaC and gammaD for the difference of exponentials model. That is,
#'  \code{params <- data.frame(gammaC = 0.98, gammaD = 0.818)}, for the difference of exponentials model.
#' @param penalty tuning parameter lambda
#' @param type type of model, must be one of AR(1) 'ar1', AR(1) + intercept 'intercept', or difference of exponentials 'dexp'
#' @param calcFittedValues  TRUE to calculate fitted values.
#'
#' @return Returns a list with elements:
#' @return \code{changePts} the set of changepoints
#' @return \code{spikes} the set of spikes
#' @return \code{fittedValues} estimated calcium concentration
#' @return \code{dat} the data that was used to fit the model
#'
#' @examples
#'
#' sampleData <- simulateAR1(n = 500, seed = 1, poisRate = 0.01, decay = 0.998, sd = 0.15)
#'
#' # AR(1) model
#'
#' fit <- segment(sampleData$fl, params = 0.998, penalty = 1, type = "ar1")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # AR(1) + intercept model
#'
#' fit <- segment(sampleData$fl, params = 0.998, penalty = 1, type = "intercept")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # Difference of exponentials model
#'
#' params <- data.frame(gammaC = 0.998, gammaD = 0.82)
#' fit <- segment(sampleData$fl, params, penalty = 1, type = "dexp")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' @export
segment <- function(dat, params, penalty, type = c("ar1", "intercept", "dexp"), calcFittedValues = TRUE) {
    if (type %in% c("ar1", "dexp", "intercept")) {
        table <- computeSegmentation(dat, params, penalty, type)
        changePts <- findChangePts(table[, 3])
        spikes <- changePts[-1] + 1
        if (calcFittedValues) {
          fittedValues <- computeFittedValues(dat, changePts, params, type)
        } else {
          fittedValues <- NULL
        }
        return(list(changePts = changePts, spikes = spikes, fittedValues = fittedValues, dat = dat))
    } else {
        stop("Model not implemented. Type must be one of ar1, dexp, or intercept.")
    }
}
