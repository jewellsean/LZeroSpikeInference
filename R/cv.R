#' Cross-validate and optimize model parameters
#' See Algorithm 3 of Jewell and Witten (2017)
#' Exact Spike Train Inference Via L0 Optimization
#'
#' @param dat fluorescence trace (a vector)
#' @param params model parameters. For the AR(1) and AR(1) intercept models this is the scalar decay parameter; this is a
#' dataframe with two parameters gammaC and gammaD for the difference of exponentials model. That is,
#' \code{params <- data.frame(gammaC = 0.98, gammaD = 0.818)}, for the difference of exponentials model.
#' @param type type of model, must be one of AR(1) 'ar1', AR(1) with intercept 'intercept', or difference of exponentials 'dexp'
#' @param optimizeParams whether the model parameters should be optimized (TRUE) or just CV of the tuning parameter (FALSE)
#' @param nLambdas number of tuning parameters to estimate the model (grid of values is automatically produced)
#' @param lambdas vector of tuning parameters to use in cross-validation
#'
#' @return A list of values corresponding to the 2-fold cross-validation:
#' @return \code{cvError} the MSE for each tuning parameter
#' @return \code{cvSE} the SE for the MSE for each tuning parameter
#' @return \code{lambdas} tuning parameters
#' @return \code{params} matrix of (optimized) parameters, rows correspond to tuning parameters, columns correspond to optimized parameter
#' @return \code{lambdaMin} tuning parameter that gives the smallest MSE
#' @return \code{lamnda1SE} 1 SE tuning parameter
#' @return \code{indexMin} the index corresponding to lambdaMin
#' @return \code{index1SE} the index corresponding to lambda1SE
#'
#' @examples
#'
#' sampleData <- simulateAR1(n = 500, seed = 1, poisRate = 0.01, decay = 0.998, sd = 0.15)
#'
#' # AR(1) model
#' # Select tuning parameter with 1 SE rule
#' cvOut <- cv(sampleData$fl, params = 0.998, type = "ar1", optimizeParams = FALSE)
#' plotCV(cvOut)
#' fit <- segment(sampleData$fl, params = 0.998, penalty = cvOut$lambda1SE, type = "ar1")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # AR(1) model
#' # Select tuning parameter with 1 SE rule, provide tuning parameters
#' cvOut <- cv(sampleData$fl, params = 0.998, type = "ar1", optimizeParams = FALSE,
#' lambdas = 10^seq(-1, 2, length.out = 10))
#' plotCV(cvOut)
#' fit <- segment(sampleData$fl, params = 0.998, penalty = cvOut$lambda1SE, type = "ar1")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # AR(1) model
#' # Select tuning parameter with 1 SE rule and optimize decay parameter
#' cvOut <- cv(sampleData$fl, params = 0.998, type = "ar1", optimizeParams = TRUE)
#' plotCV(cvOut)
#' fit <- segment(sampleData$fl, params = cvOut$params[cvOut$index1SE, 1],
#' penalty = cvOut$lambda1SE, type = "ar1")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # AR(1) + intercept model
#' # Select tuning parameter with 1 SE rule and optimize decay parameter
#' cvOut <- cv(sampleData$fl, params = 0.998, type = "intercept", optimizeParams = TRUE)
#' plotCV(cvOut)
#' fit <- segment(sampleData$fl, params = cvOut$params[cvOut$index1SE, 1],
#' penalty = cvOut$lambda1SE, type = "ar1")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#' # Difference of exponentials model
#' # Select tuning parameter with 1 SE rule and optimize both model parameters
#' params <- data.frame(gammaC = 0.998, gammaD = 0.818)
#' cvOut <- cv(sampleData$fl, params = params, type = "dexp", optimizeParams = TRUE)
#' plotCV(cvOut)
#' optimParams <- data.frame(t(cvOut$params[cvOut$index1SE, ]))
#' fit <- segment(sampleData$fl, params = optimParams, penalty = cvOut$lambda1SE, type = "dexp")
#' plotSegmentation(fit, trueSpikes = sampleData$spikes)
#'
#'
#' @export
#'
cv <- function(dat, params, type = "ar1", optimizeParams = TRUE, nLambdas = 10, lambdas = NULL) {
    k <- 2  ## number of folds
    n <- length(dat)

    if (is.null(lambdas)) {
      lambdas <- createLambdaSequence(n, nLambdas)
    }

    ## Modified parameters for CV
    paramsTilde <- modifyParams(params, type, "fwd")

    foldid <- rep(seq(1, k), n)[1:n]
    cvMSE <- matrix(0, nrow = nLambdas, ncol = k)

    ## store all intermediate values for model parameters
    nParams <- length(params)
    if (optimizeParams) {
        paramsOut <- list()
        for (i in 1:nParams) paramsOut[[i]] <- matrix(0, nrow = nLambdas, ncol = k)
    } else paramsOut <- NULL

    for (fold in 1:k) {
        trainInd <- which(foldid != fold)
        testInd <- which(foldid == fold)

        trainDat <- dat[trainInd]
        testDat <- dat[testInd]
        nn <- length(trainInd)

        for (lambdaInd in 1:nLambdas) {
            if (optimizeParams) {
                segments <- segment(trainDat, paramsTilde, lambdas[lambdaInd], type, calcFittedValues = FALSE)
                paramsTilde <- optimParams(paramsTilde, trainDat, segments$changePts, lambdas[lambdaInd],
                  type)
            }
            segments <- segment(trainDat, paramsTilde, lambdas[lambdaInd], type)
            yhatTrain <- segments$fittedValues
            yhatTest <- c(yhatTrain[1], 0.5 * (yhatTrain[1:(nn - 1)] + yhatTrain[2:nn]))
            cvMSE[lambdaInd, fold] <- mean( (yhatTest - testDat) ^ 2)

            if (optimizeParams) {
                for (i in 1:nParams) paramsOut[[i]][lambdaInd, fold] <- as.numeric(paramsTilde[i])
            }
        }
    }

    cvErr <- rowMeans(cvMSE)
    cvse <- apply(cvMSE, 1, sd) / sqrt(k)

    indexMin <- which.min(cvErr)
    lambdaMin <- lambdas[indexMin]
    lambda1SE <- max(lambdas[cvErr <= cvErr[indexMin] + cvse[indexMin]])
    index1SE <- which(lambdas == lambda1SE)

    if (optimizeParams) {
        paramsOutTmp <- matrix(0, ncol = nParams, nrow = nLambdas)
        for (i in 1:nParams) paramsOutTmp[, i] <- as.matrix(modifyParams(rowMeans(paramsOut[[i]]),
            type, "bck"))
        paramsOut <- paramsOutTmp
        colnames(paramsOut) <- colnames(params)
    }

    out <- list(cvError = cvErr, cvSE = cvse, lambdas = lambdas, params = paramsOut, lambdaMin = lambdaMin,
                lambda1SE = lambda1SE, indexMin = indexMin, index1SE = index1SE)

    return(out)

}

yhatMSE <- function(params, y, changePts, type) {
    if (type == "dexp") {
        params <- data.frame(gammaC = params[1], gammaD = params[2])
        return(mean( (y - computeFittedValues(y, changePts, params, type)) ^ 2))
    }
    return(mean( (y - computeFittedValues(y, changePts, params, type)) ^ 2))
}

createLambdaSequence <- function(n, nLambdas = 10) {
    return(10^(seq(-1, 1.5, length.out = nLambdas)))
}
#
# Optimize model parameters based on a fixed tuning parameter (and spike times).
# @param params model parameters
# @param dat flourscence
# @param changePts locations of changepoints
# @param penalty tuning parameter
# @param type type of model
optimParams <- function(params, dat, changePts, penalty, type) {
    if (type %in% c("ar1", "intercept")) {
        return(optimize(f = yhatMSE, interval = c(0.9, 1), y = dat, changePts = changePts, type = type)$minimum)
    }

    if (type == "dexp") {
        params <- t(as.matrix(params))
        out <- constrOptim(theta = params, f = yhatMSE, ui = matrix(c(1, -1, 0, 0, 0, 0, 1, -1),
            nrow = 4), grad = NULL, ci = matrix(c(0, -1, 0, -1), nrow = 4), y = dat, changePts = changePts,
            type = type)$par
        return(data.frame(gammaC = out[1], gammaD = out[2]))
    }
}

modifyParams <- function(params, type, direction) {
    if (type %in% c("ar1", "intercept")) {
        if (direction == "fwd")
            return(params^2)
        if (direction == "bck")
            return(sqrt(params))
    }
    if (type == "dexp") {
        if (direction == "fwd")
            return(data.frame(params^2))
        if (direction == "bck")
            return(data.frame(sqrt(params)))
    }
}


