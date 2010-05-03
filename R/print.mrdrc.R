print.mrdrc <- function(x, ...) {
    cat(x$result)
}

"pns.print1" <- function(x, ...) {
    mobject <- x  # $model

    appendFct <- function(x, y) {paste(x, y, sep = "")}

    textStr <- "USER CHOICES:\n\n"

    textStr <-appendFct(textStr, paste("Data type: ", mobject$type, "\n\n"))

    modelArg <- mobject$model
    textStr <- appendFct(textStr, paste("Model: ", modelArg, "\n\n"))
    if (modelArg != "non-parametric")
    {
        textStr <- appendFct(textStr, "Parametric components: \n\n")
        lapFct1 <- function(fct){fct$name}  # could also be fct$text
        textStr <- appendFct(textStr, paste(as.vector(unlist(lapply(mobject$fct, lapFct1))), collapse = ", "))
        textStr <- appendFct(textStr, "\n\n")
        robustYesNo <- ifelse(mobject$robust, "yes", "no")
        textStr <- appendFct(textStr, paste("Robust estimation: ", robustYesNo, "\n\n"))
    }

    textStr <- appendFct(textStr, paste("ED values: ", paste(mobject$respLev, collapse = ", "), "\n"))
    referenceArg <- mobject$reference
    refText <- ifelse(is.null(referenceArg), "none", referenceArg)
    textStr <- appendFct(textStr, paste("Reference: ", refText, "\n"))
    textStr <- appendFct(textStr, paste("Confidence level: ", mobject$level, "\n\n\n"))

    if ( (modelArg != "non-parametric") && (length(mobject$fct) > 1) )
    {
        textStr <- appendFct(textStr, "Ranked models with estimated ED values:\n\n")
    } else {
        textStr <- appendFct(textStr, "Estimated ED values:\n\n")
    }

    ## Inserting model names and tables
    if (modelArg != "non-parametric")
    {
    for (i in 1:length(mobject$fct))
    {
        mfitsElement <- mobject$fits[[i]]
        tempStr1 <- paste(mfitsElement$text, " (", mfitsElement$name, ")", sep = "")
        if (modelArg == "semi-parametric")
        {
            tempStr1 <- paste(tempStr1, "\n", "(mixing: ", as.character(mfitsElement$mixing), ")", sep = "")
        }
        tempStr2 <- paste(tempStr1, "\n\n", "$XTABLE_", i, "$\n\n", sep = "")
        textStr <- appendFct(textStr, tempStr2)
    }
    } else {
        textStr <- appendFct(textStr, paste("\n\n", "$XTABLE_1", sep = ""))
    }

    ## Returning the result string
    textStr
}


"pns.print2" <- function(x, ...) {
    mobject <- x  # $model

    lapFct2 <- function(fitsElement)
    {
        fctName <- fitsElement$name
        formem <- format(fitsElement$edMat, digits = 4)
        formaic <- format(fitsElement$aic, digits = 4)
        formaicd <- format(fitsElement$aicDiff, digits = 4)

        lendf <- nrow(formem)
        aicVec <- c(formaic, paste("(", formaicd, ")"), rep("", lendf - 2))
        respVec <- rownames(formem)

        retDF <- data.frame(aicc = aicVec, level = as.numeric(respVec), ed = fitsElement$edMat)
        names(retDF) <- c("AICC", "Level", "Est", "Lower", "Upper")
        return(retDF)
    }
    return(lapply(mobject$fits, lapFct2))
}
