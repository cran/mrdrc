"feedSP" <- function(data, fct, logScale = TRUE)
{
    if (logScale) {expFct <- exp} else {expFct <- function(x){x}} 
    cases <- data$cases
    covariate <- expFct(data$covariate)  # covariate is on log scale
    response <- data$response
    data2 <- data
    data2$covariate <- expFct(data$covariate)  # covariate is on log scale

    if (!is.null(cases))
    {
        drmFit <- drm(response ~ covariate, weights = cases, data = data2, fct = fct, type = "binomial")
    } else {
        drmFit <- drm(response ~ covariate, fct = fct, data = data2)
    }
    
    poVec <- predict(drmFit, se.fit = FALSE)
    leny <- length(response)
    if (!is.null(cases))
    {
        poVec[poVec < 0.01] <- 0.01  # to avoid problems in 'Wmat' below
        poVec[poVec > 0.99] <- 0.99

        Wmat <- diag(as.vector(cases / (poVec * (1 - poVec))))
        # diag() needs a vector, not a column matrix
    } else {    
        Wmat <- diag(1, leny)
    }
    predFct <- function(new.data) 
    {
        new.data$covariate <- expFct(new.data$covariate)  # covariate is on log scale
        as.vector(predict(drmFit, new.data, se.fit = FALSE))
    }
    
    looPredict <- function(i) 
    {
        as.vector(predict(update(drmFit, data = data2[-i, ]), 
        data.frame(covariate = covariate[i]), se.fit = FALSE))
    }
    loopVec <- sapply(1:leny, looPredict)

#    Dmat <- m1drc$deriv1
    Dmat <- drmFit$deriv1
    hmFct <- function(new.data)
    {
#        print(Dmat)
#        print(Wmat)
        x <- expFct(new.data$covariate)  # covariate is on log scale
        pMat0 <- t(drmFit$parmMat) 
        pMat <- matrix(pMat0, ncol = ncol(pMat0), nrow = length(x), byrow = TRUE)

        drmFit$fct$deriv1(x, pMat) %*% solve(t(Dmat)%*%Wmat%*%Dmat) %*% t(Dmat) %*% Wmat 
    }
    
    res <- list(
    hatMatrix = hmFct, 
    hatMatrixObs = hmFct(data.frame(covariate = data$covariate)),  
    leave.one.out.predict = function() {loopVec},
    leaveOneOutObs = loopVec,
    modelDf = leny - df.residual(drmFit), 
    predict = predFct, 
    predictObs = poVec)
    
    structure(res, class = "P")
}
