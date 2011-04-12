"pnsdrm" <-  function(predictor, response, weights, type = c("continuous", "binomial"),
model = c("semi-parametric", "non-parametric", "parametric"), fct = NULL, robust = FALSE, respLev = c(10, 20, 50),
reference = NULL, level = 0.95, logex = FALSE)
{
#    require(drc, quietly = TRUE)  # as package dependency

    type <- match.arg(type)
    type2 <- ifelse(type == "continuous", "con", "bin")  # used in mrdrm
    model <- match.arg(model)

    ## Capturing the variable names (before manipulations begin)
    respName <- deparse(substitute(response))
    predName <- deparse(substitute(predictor))
    if (type == "binomial")
    {
        respName <- paste(respName, "/", deparse(substitute(weights)), sep = "")
    }
   
    ## Checking argument values
    if ( (type == "binomial") && (missing(weights)) )
    {
        stop("Weights need to be supplied for binomial data")
    }

    ## Converting to proportions in the binomial case
    if (type == "binomial")
    {
        response <- response / weights
    }
    
    ## Converting argument 'fct' to a list
    if (!is.list(fct[[1]]))
    {
        fct <- list(fct)
    }
   
    ## Fitting model(s) and calculating ED estimates
    resList <- switch(model, 
    "non-parametric" = {
#        detach(package:drc)
#        require(EffectiveDose, quietly = TRUE)

##        npdrm(predictor, response, respLev, level)
        stop("Non-parametric approach currently not available")
    },
    "parametric" = {
#        fctList <- list()    
#        if (is.list(fct[[1]]))
#        {
            lenFct <- length(fct)
            fctList <- list()
#            nameVec <- rep("", lenFct)
            for (i in 1:lenFct)
            {
                if (type == "binomial")
                {
                    tempFit <- drm(response~predictor, weights = weights, fct = fct[[i]], type = type,
                    na.action = na.omit)
                }
                if (type == "continuous")
                {
                    tempFit <- drm(response~predictor, fct = fct[[i]], type = type, na.action = na.omit)
                }
                if (identical(substr(fct[[i]]$"name", 1, 3), "LL2"))
                {
                    edMat <- drc:::ED(tempFit, respLev, interval = "fls", level = level, display = FALSE)
                } else {
                    edMat <- drc:::ED(tempFit, respLev, interval = "delta", level = level, display = FALSE)[, c(1, 3, 4)]
                }
                
                fctList[[i]] <- list(fit = tempFit, 
                edMat = edMat,
                aic = AIC(tempFit), mixing = 0)
#                fctList[[i]]$fit <- tempFit
#                fctList[[i]]$edMat <- ED(tempFit, respLev, ci = "delta", level = level)[, c(1, 3, 4)]
#                fctList[[i]]$aic <- AIC(tempFit) 
#                fctList[[i]]$mixing <- 0
#                nameVec[i] <- fct[[i]]$name
            }
#            names(fctList) <- nameVec  # used?
            
#            fctList
#        } else {    
#            tempFit <- drm(response~predictor, fct = fct, type = type)
#            fctList <- list(list(fit = tempFit, 
#            edMat = drc:::ED(tempFit, respLev, ci = "delta", level = level, display = FALSE)[, c(1, 3, 4)],
#            aic = AIC(tempFit), mixing = 0)) 
##            fctList[[1]]$fit <- tempFit
##            fctList[[1]]$edMat <- ED(tempFit, respLev, ci = "delta", level = level)[, c(1, 3, 4)]
##            fctList[[1]]$aic <- AIC(tempFit) 
##            fctList[[1]]$mixing <- 0
##            names(fctList) <- fct$name
##           
##            fctList
#        }
        fctList        
    },
    "semi-parametric" = {
        fctList <- list()    
#        if (is.list(fct[[1]]))
#        {
            lenFct <- length(fct)
#            fctList <- list()
#            nameVec <- rep("", lenFct)
            for (i in 1:lenFct)
            {
                fctList[[i]] <- mrdrm(predictor, response, weights, type2, fct[[i]], 
                respLev, reference, level, robust, logex = logex)
#                nameVec[i] <- fct[[i]]$name
            }
#            names(fctList) <- nameVec
            
#            fctList
#        } else {
#            fctList[[1]] <- mrdrm(predictor, response, weights, type2, fct, respLev, reference, level, robust)
##            names(fctList) <- fct$name
#        }
        fctList
    }
    )
   
#    pnsList <- list(fits = resList, model = model, fct = fct, type = type, respLev = respLev, 
#    reference = reference, level = level, robust = robust,
#    status = 0, error.msg = "", headline = "Dose-response analysis")  # simply hardcoded at the moment   

    pnsList <- list(status = 0, error.msg = "", headline = "Dose-response analysis")
    mobject1 <- list(fits = resList, model = model, fct = fct, type = type, respLev = respLev, 
    reference = reference, level = level, robust = robust)

    ## Ordering the fits
    updateFits <- function(mobject)
    {
        mfits <- mobject$fits
        
        aicVec <- as.vector(unlist(lapply(mfits, function(x) {x$aic})))
        aicVec2 <- aicVec - min(aicVec)
        for (i in 1:length(mfits))
        {
            mfits[[i]]$aicDiff <- aicVec2[i] 
            if (model != "non-parametric") 
            {
                mfits[[i]]$name <- mobject$fct[[i]]$name
                mfits[[i]]$text <- mobject$fct[[i]]$text
            }
        }
        mobject$fits <- mfits[order(aicVec)]

        mobject
    } 
    mobject2 <- updateFits(mobject1)
    mobject2$resp <- respName  
    mobject2$pred <- predName
    pnsList$model <- mobject2

    # Creating output: text and tables    
    pnsList$result <- pns.print1(mobject2)  # pns.print1(pnsList)
    pnsList$xtables <- pns.print2(mobject2)  # pns.print2(pnsList)
  
    class(pnsList) <- "mrdrc"
    pnsList
}

## Defining function for the GUI
"pnsdrm.calc" <- pnsdrm









