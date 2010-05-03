"plot.mrdrc" <- function(x, ...) {
    object <- x
    pns.plot1(object$model, ...)
}

"pns.plot1" <- function(model, log = "x", ...) 
{
    mobject <- model

    if (mobject$model == "semi-parametric")
    {
        plotFct <- function(fitsElement)
        {
            plot(fitsElement$fit, with.NP = FALSE, with.P = FALSE,
            xlab = mobject$pred, ylab= mobject$resp, main = paste("Mixing:", mobject$fits[[1]]$mixing),
            logPred = mobject$fits[[1]]$logPred, lscale = log)
        }
    }
    if (mobject$model == "parametric")
    {
        plotFct <- function(fitsElement) {plot(fitsElement$fit, broken = TRUE, log = log, ...)}
    }
    if (mobject$model == "non-parametric")
    {
        return(list(status = 1, error.msg = "No plot available", caption = mobject$"model"))
    }
    
    lenMob <- length(mobject$fct)
    noCR <- ceiling(sqrt(lenMob))
    if (lenMob > 2) 
    {
        par(mfrow = c(noCR, noCR))
    } else {
        par(mfrow = c(1, lenMob))
    }
    lapply(mobject$"fits", plotFct)
    par(mfrow = c(1, 1))

    return(list(status = 0, error.msg = "", caption = mobject$"model"))
}
