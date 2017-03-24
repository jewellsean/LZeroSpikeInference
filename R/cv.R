#' Cross-validate and optimize model parameters
#'
#' @details
#' We perform cross-validation over a one-dimensional grid of \eqn{\lambda} values.
#'  For each value of \eqn{\lambda} in this grid, we solve the corresponding optimization problem, that is, one of
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
#' minimize_{c1,...,cT,d1,...,dT} 0.5 sum_{t=1}^{T} ( y_t - (c_t - d_t) )^2 + lambda sum_{t =2}^{T} 1_{ c_t neq gamma _c c_{t-1}, d_t neq gamma_d d_{t-1} }
#' on a training set using a candidate value for \eqn{\gamma}. Given the resulting set of changepoints, we solve a constrained optimization problem for \eqn{\gamma}. We then refit the optimization problem with the optimized value of \eqn{\gamma},
#' and then evaluate the mean squared error (MSE) on a hold-out set. Note that in the final output of the algorithm,
#' we take the square root of the optimal value of \eqn{\gamma} in order to address the fact that the cross-validation
#' scheme makes use of training and test sets consisting of alternately-spaced timesteps.
#'
#' If there is a tuning parameter lambdaT in the path [lambdaMin, lambdaMax] that produces a fit with
#' less than 1 spike per 10,000 timesteps the path is truncated to [lambdaMin, lambdaT] and a warning is produced.
#'
#' See Algorithm 3 of Jewell and Witten (2017)
#' Exact Spike Train Inference Via L0 Optimization
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
#'
#' sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' # AR(1) model
#' outAR1 <- cv.estimateSpikes(sim$fl, type = "ar1")
#' plot(outAR1)
#' print(outAR1)
#' fit <- estimateSpikes(sim$fl, gam = outAR1$optimalGam[outAR1$index1SE, 1],
#' lambda = outAR1$lambda1SE, type = "ar1")
#' plot(fit)
#' print(fit)
#'
#' # AR(1) + intercept model
#' outAR1Intercept <- cv.estimateSpikes(sim$fl, type = "intercept",
#' lambdas = seq(0.1, 5, length.out = 10))
#' plot(outAR1Intercept)
#' print(outAR1Intercept)
#' fit <- estimateSpikes(sim$fl, gam = outAR1Intercept$optimalGam[outAR1Intercept$index1SE, 1],
#' lambda = outAR1Intercept$lambda1SE, type = "intercept")
#' plot(fit)
#' print(fit)
#'
#' # Difference of exponentials model
#' sim <- simulateDexp(n = 500, gams = c(0.998, 0.7), poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#' outCVDexp <- cv.estimateSpikes(sim$fl, type = "dexp", gam = c(0.998, 0.7))
#' plot(outCVDexp)
#' fit <- estimateSpikes(sim$fl, gam = c(0.998, 0.7), lambda = outCVDexp$lambda1SE, type = "dexp")
#' plot(fit)
#' print(fit)
#' @param dat fluorescence trace (a vector)
#' @param type type of model, must be one of AR(1) 'ar1', AR(1) with intercept 'intercept', or difference of exponentials 'dexp'
#' @param gam a scalar value for the AR(1)/AR(1) + intercept decay parameter or a vector (gammaC, gammaD) for the dexp model. If no value is provided an optimal value is selected for each tuning parameter
#' @param nLambdas number of tuning parameters to estimate the model (grid of values is automatically produced)
#' @param lambdas vector of tuning parameters to use in cross-validation
#'
#' @return A list of values corresponding to the 2-fold cross-validation:
#' @return \code{cvError} the MSE for each tuning parameter
#' @return \code{cvSE} the SE for the MSE for each tuning parameter
#' @return \code{lambdas} tuning parameters
#' @return \code{optimalGam} matrix of (optimized) parameters, rows correspond to tuning parameters, columns correspond to optimized parameter
#' @return \code{lambdaMin} tuning parameter that gives the smallest MSE
#' @return \code{lambda1SE} 1 SE tuning parameter
#' @return \code{indexMin} the index corresponding to lambdaMin
#' @return \code{index1SE} the index corresponding to lambda1SE
#'
#' @export
#'
cv.estimateSpikes <- function(dat, type = "ar1", gam = NULL,
                              lambdas = NULL, nLambdas = 10) {
    k <- 2  ## number of folds
    n <- length(dat)

    if (is.null(lambdas)) {
      lambdas <- createLambdaSequence(n, nLambdas)
    }

    if (is.null(gam))
    {
      optimizeGams = TRUE
      ## Modified parameters for CV
      if (type == "dexp") {
        params <- c(0.99, 0.82)
      } else {
        params <- 0.998
      }
    } else {
      optimizeGams = FALSE
      params <- gam
    }

    paramsTilde <- modifyParams(params, type, "fwd")

    foldid <- rep(seq(1, k), n)[1:n]
    cvMSE <- matrix(0, nrow = nLambdas, ncol = k)

    ## store all intermediate values for model parameters
    nParams <- length(params)
    if (optimizeGams) {
        paramsOut <- list()
        for (i in 1:nParams) paramsOut[[i]] <-
            matrix(0, nrow = nLambdas, ncol = k)
    } else paramsOut <- gam

    for (fold in 1:k) {
        trainInd <- which(foldid != fold)
        testInd <- which(foldid == fold)

        trainDat <- dat[trainInd]
        testDat <- dat[testInd]
        nn <- length(trainInd)

        for (lambdaInd in 1:nLambdas) {
            if (optimizeGams) {
                segments <- estimateSpikes(trainDat, paramsTilde,
                                           lambdas[lambdaInd], type,
                                           calcFittedValues = FALSE)
                paramsTilde <- optimParams(paramsTilde, trainDat,
                                           segments$changePts,
                                           lambdas[lambdaInd], type)
            }
            segments <- estimateSpikes(trainDat, paramsTilde,
                                       lambdas[lambdaInd], type)
            yhatTrain <- segments$fittedValues
            yhatTest <- c(yhatTrain[1], 0.5 * (yhatTrain[1:(nn - 1)] +
                                                 yhatTrain[2:nn]))
            cvMSE[lambdaInd, fold] <- mean( (yhatTest - testDat) ^ 2)

            if (optimizeGams) {
                for (i in 1:nParams) paramsOut[[i]][lambdaInd, fold] <-
                    as.numeric(paramsTilde[i])
            }
            if (length(segments$spikes) < nn / 10000 &&
                lambdaInd < nLambdas) {
              warning("Cross validation stopped early as less than 1 spike per 10,000 timesteps estimated. Rerun with smaller lambdas.")
              lambdas <- lambdas[1:lambdaInd]
              nLambdas <- length(lambdas)
              break
            }

        }
    }




    cvErr <- rowMeans(cvMSE[1: lambdaInd, ])
    cvse <- apply(cvMSE[1: lambdaInd, ], 1, sd) / sqrt(k)

    indexMin <- which.min(cvErr)
    lambdaMin <- lambdas[indexMin]
    lambda1SE <- max(lambdas[cvErr <= cvErr[indexMin] + cvse[indexMin]])
    index1SE <- which(lambdas == lambda1SE)

    if (optimizeGams) {
        paramsOutTmp <- matrix(0, ncol = nParams, nrow = nLambdas)
        for (i in 1:nParams)
          paramsOutTmp[, i] <- as.matrix(modifyParams(rowMeans(paramsOut[[i]][1: lambdaInd, ]),
            type, "bck"))
        paramsOut <- paramsOutTmp
        colnames(paramsOut) <- colnames(params)
    }

    out <- list(cvError = cvErr, cvSE = cvse, lambdas = lambdas,
                optimalGam = paramsOut, lambdaMin = lambdaMin,
                lambda1SE = lambda1SE, indexMin = indexMin,
                index1SE = index1SE,
                call = match.call(),
                optimized = optimizeGams,
                type = type)
    class(out) <- "cvSpike"
    return(out)

}

yhatMSE <- function(params, y, changePts, type) {
  return(mean( (y - computeFittedValues(y, changePts, params, type)) ^ 2))
}

createLambdaSequence <- function(n, nLambdas = 10) {
  return(10^(seq(-1, 1, length.out = nLambdas)))
}

optimParams <- function(params, dat, changePts, penalty, type) {
    if (type %in% c("ar1", "intercept")) {
        return(optimize(f = yhatMSE, interval = c(0.9, 1), y = dat, changePts = changePts, type = type)$minimum)
    }

    if (type == "dexp") {
        params <- as.matrix(params)
        out <- constrOptim(theta = params, f = yhatMSE, ui = matrix(c(1, -1, 0, 0, 0, 0, 1, -1),
            nrow = 4), grad = NULL, ci = matrix(c(0, -1, 0, -1), nrow = 4), y = dat, changePts = changePts,
            type = type)$par
        return(out)
    }
}

modifyParams <- function(params, type, direction) {
    if (type %in% c("ar1", "intercept", "dexp")) {
        if (direction == "fwd")
            return(params^2)
        if (direction == "bck")
            return(sqrt(params))
    }
}