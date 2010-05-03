"mrdrm" <- function(
predictor, response, weights, type = c("bin", "con"), fct, respLev, reference, level,
robust, mixVec = seq(0, 1, by = 0.05), logex = FALSE, bwLower = 0, compact = TRUE) 
{
    type <- match.arg(type)

    ## Using log-transformed predictor?
    logPred <- (min(predictor) > 1e-6)
    if (logPred)
    {
        logScale <- TRUE
        predictor <- log(predictor)
    } else {
        logScale <- FALSE
    }
    if (logex)
    {
        logScale <- TRUE
        predictor <- logExtend(predictor)
        logPred <- TRUE
    }

    ## Setting range for bandwidth
    diffVec <- abs(diff(predictor))
    diffVec <- diffVec[diffVec > 1e-6]
    bandwSeq <- seq(max(c(bwLower, min(diffVec))), max(diffVec), length.out = 20)

    ## Setting non-parametric control
    npCon <- NP.control.lr.wls(
    response.type = type,
    optim = "bandwidth.grid",
    weight.function = "gaus",
    bandwidth.type = "fixed.width",
    bandwidth = bandwSeq)

    ## Setting parametric control ... not used
    pCon <- P.control.glm.binomial.ml(link = "logit")

    ## Setting model-robust control
    mrCon <- SP.control.mrr(
    response.type = type,
    optim = "mixing.grid",
    mixing = mixVec)

    ## Fitting model-robust model
    if (missing(weights))
    {
        mrrFit <- SP.mrr(formula = response~predictor,
        data = data.frame(predictor = predictor, response = response),
        NP.control = npCon, P.control = pCon, SP.control = mrCon, logScale = logScale,
        general = TRUE, fct = fct, robust = robust, compact = compact)
    } else {
        mrrFit <- SP.mrr(formula = response~predictor,
        data = data.frame(predictor = predictor, response = response), cases = weights,
        NP.control = npCon, P.control = pCon, SP.control = mrCon, logScale = logScale,
        general = TRUE, fct = fct, robust = robust, compact = compact)
    }

    ## Calculating ED values
    lenrl <- length(respLev)
    edMat <- matrix(NA, lenrl, 3)
    for (i in 1:lenrl)
    {
        edRes <- ED.H(mrrFit, respLev[i], reference, level)
        edMat[i, 1] <- edRes$ed.alpha
        edMat[i, 2] <- edRes$ed.alpha.lower
        edMat[i, 3] <- edRes$ed.alpha.upper
    }
    if (logScale) {edMat <- exp(edMat)}
    colnames(edMat) <- c("Estimate", "Lower", "Upper")
    rownames(edMat) <- respLev

    ## Calculating AIC (maybe it should be AICC?)
    response <- as.vector(na.omit(response))
    nVal <- length(response)
    rssVal <- sum( (response - mrrFit$predictObs)^2 )
    aicVal <- nVal * log(2*pi) + nVal * log(rssVal/nVal) + nVal + 2 * mrrFit$modelDf

    list(edMat = edMat, aic = aicVal, fit = mrrFit, mixing = mrrFit$mixing, bandwidth = mrrFit$NP$bandwidth,
    logPred = logPred)
}

"logExtend" <- function(x) {
    xOrder <- order(x)
    x <- x[xOrder]

    xoIndex <- (x > 1e-12)
    xOver <- x[xoIndex]
    increm <- abs(mean(diff(log(unique(xOver)))))

    returnX <- x
    lxo <- log(xOver)
    returnX[xoIndex] <- lxo
    returnX[!xoIndex] <- min(lxo) - increm

    returnX[xOrder] <- returnX

    returnX
}

"plotmr" <- function(object, log = "", ...) {
    plot(object$"fit", with.NP = FALSE, with.P = FALSE, lscale = log, ...)  # CI.H=FALSE, with.ED=FALSE, alpha=50
}
