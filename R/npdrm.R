"npdrm" <- function(predictor, response, respLev, level) {
    ## Using log-transformed predictor?
    logPred <- (min(predictor) > 1e-6)
    if (logPred)
    {
        logScale <- TRUE
        predictor <- log(predictor)
    } else {
        logScale <- FALSE
    }
#    print(logScale)

    ## Determining whether or not the curve is decreasing
    if (coef(lm(response~predictor))[2] > 0)
    {
        monoArg <- "increasing"
    } else {
        monoArg <- "decreasing"
#        respLev <- 100 - respLev  # not needed because type = "prob"
        # to adjust for how ED() in 'EffectiveDose'
    }

    lenrl <- length(respLev)
    edMat <- matrix(NA, lenrl, 3)
    listXY <- list(x = predictor, y = response)

#    bwd <- 0.5
    edFit <- EffectiveDose:::ED(listXY, alpha = respLev/100, mono = monoArg, type = "prob")
    edFit@call$mono <- monoArg  # a little fix

    edMat[, 1] <- edFit@ED
    for (i in 1:lenrl)
    {
#        edMat[i, 1] <- EffectiveDose:::ED(listXY, alpha = respLev[i]/100)@ED
        edMat[i, 2:3] <- as.vector(EffectiveDose:::Boot.CI(listXY, alpha = respLev[i]/100, level = level,
        mono = monoArg, type = "prob")@CI)
    }
    if (logScale) {edMat <- exp(edMat)}
    colnames(edMat) <- c("Estimate", "Lower", "Upper")
    rownames(edMat) <- respLev

    list(list(edMat = edMat,
    aic = as.vector(EffectiveDose:::aic.ED.locfit(edFit)[4]),
    fit = edFit, mixing = 1, bandwidth = NA, logPred = logPred))
}
