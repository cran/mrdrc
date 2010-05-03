# Parametric engine for binomial data to be used in model robust regression
# data (class data.frame, object): observed data points, must hold varaibles
#   (covariate,response,cases)
# control (class P.ctr, object): controls for this engine
# leave.one.out (boolean, atomic): it true do leave-one-out computations
P.glm.binomial.ml <- function (
  data,
  control=P.control.glm.binomial.ml(),
  leave.one.out = TRUE
) {

##  message ("\n+++ P.glm.binomial.ml: creating object")

  # Handle arguments
  if ( !inherits(control,"P.ctr")) {
    stop ("P.glm.binomial.ml: non 'P.ctr' class control")
  }
  if ( control$response.type!="binomial" ) {
    stop ("P.glm.binomial.ml: control is not for 'binomial' response")
  }
  if (
    length(
      match(
        c("covariate","response","cases"),
        names (data))
    ) != 3
  ) {
    stop ("P.glm.binomial.ml: not all variables 'covariate','response','cases' in data")
  }
  family <- switch (
    control$link,
    "logit" =,
    "probit" =,
    "cauchit" = {
      binomial (control$link)
    },
    stop ("P.glm.binomial.ml: invalid link")
  )

  # Fit engine
  fitter <- function (subset) {
    cl <- call ("glm")
    cl$formula <- response ~ covariate
    cl$family <- family
    cl$weights <- data$cases
    cl$data <- data
    if ( !missing(subset) ) {  cl$subset <- subset  }
    eval (cl)
  }

  # Fit for observed data
  binomial.fit <- fitter ()

  # Predict engine
  predictImpl <- function (new.data) {
##    message ("P.glm.binomial.ml: predict implementation")
    predictNew <- predict (binomial.fit, new.data, type="response")
    names (predictNew) <- NULL
    predictNew
  }
  
  ## Inserted 2/2 2008 by Christian Ritz 
  ##  returning predictions on the linear predictor scale
  ##  in contrast to predictions on the probability scale
  ##  in predictImpl() above
  predict2 <- function (new.data) {
    predictNew <- predict (binomial.fit, new.data)
    names (predictNew) <- NULL
    predictNew
  }
  predictObs <- predictImpl (data)

  # Leave-one-out analysis
  leave.one.out.predictImpl <- function () {
##    message ("P.glm.binomial.ml: leave one out predict implementation")
    # TODO: Handle too few data points to produce fit
    n <- nrow(data)
    is <- 1:n
    sapply (
      1:n,
      function (index) {
        leave.out.fit <- fitter (subset=is[-index])
        predict (
          leave.out.fit,
          data.frame(covariate=data$covariate[index]),
          type="response"
       )
      }
    )
  }

  # Fixed matrices for hat matrix
  d.linkinv <- switch (
    control$link,
    "logit" = {
       dlogis
    },
    "probit" = {
       dnorm
    },
    "cauchit" = {
      dcauchy
    },
    stop ("P.glm.binomial.ml: invalid link for inverse")
  )
  X <- cbind(1,data$covariate)
  # This is wrong - instead use GLM fit to get
  # var.Y.obs <- (data$response*(1-data$response)) / data$cases
  # var.Y.inv.obs <- 1 / var.Y.obs
  # var.Y.inv.obs[is.infinite(var.Y.inv.obs)] <- 0
  predict.se <- predict (binomial.fit, data, type="response", se.fit = TRUE)

# Line below modified 2/2 2008 by Christian Ritz
#  var.Y.obs <- predict.se$se.fit^2  # wrong!!!

  var.Y.obs <- predict.se$fit * (1 - predict.se$fit) / data$cases
  var.Y.inv.obs <- 1 / var.Y.obs

  # f in MMR note is inverse link

# Line below modified 2/2 2008 by Christian Ritz
#  d.f.obs <- d.linkinv (predictObs)  # argument need to be on predictor scale
  d.f.obs <- d.linkinv(predict2(data))
  
  W <- diag( d.f.obs^2 * var.Y.inv.obs )
  hat.base <- solve(t(X)%*%W%*%X) %*% t(X) %*% W %*% diag(1/d.f.obs)

  # Hat matix for data points of interest
  hatMatrixImpl <- function (new.data) {
##    message ("P.glm.binomial.ml: hat matrix implementation")
    x <- cbind (1, new.data$covariate)

# Line below modified 2/2 2008 by Christian Ritz
#    d.f.new <- d.linkinv (predictImpl(new.data))
    d.f.new <- d.linkinv (predict2(new.data))
    
    d.f.new * x %*% hat.base
  }
  hatMatrixObs <- hatMatrixImpl (data)

  # Confidence interval based on H-matrix
  CI.HImpl <- function (new.data) {
##    message ("P.glm.binomial.ml: confidence interval implementation")
    ys <- predictImpl ( new.data )
    # TODO: There might be discrepancy btw implementation here and note on MRR
    H <- hatMatrixImpl ( new.data )
    sd.mu <- sqrt ( diag ( H %*% diag(var.Y.obs) %*% t(H) ) )
    list (
      lower=ys-1.96*sd.mu,
      upper=ys+1.96*sd.mu
    )
  }

  printImpl <- function (P.obj, verbose = FALSE) {
##    message ("P.glm.binomial.ml: print implementation")
##    message ("theta")
##    print (P.obj$theta)
##    message ("Predicted at obs")
##    print (P.obj$predictObs)
    if (verbose) {
##      message ("Model df")
##      print (P.obj$modelDf)
##      message ("Binomial likelihood fit")
##      b.f <- binomial.fit
##      b.f$call <- NULL
##      print (b.f)
##      message ("Hat matrix at obs")
##      print (P.obj$hatMatrixObs)
##      message ("Leave out one predictions at obs")
##      print (P.obj$leaveOneOutObs)
    }
  }

  plotImpl <- function (
    P.obj,
    leave.one.out = FALSE,
    CI.H = FALSE,
    CI.glm = FALSE,
    n=101
  ) {

    xmin <- min (data$covariate)
    xmax <- max (data$covariate)
    xs <- seq (xmin, xmax, length.out=n)
    xlim <- c (xmin, xmax)

    ys <- predictImpl (new.data=data.frame(covariate=xs))
    if (leave.one.out) {
      all.ys <- c(ys,data$response,P.obj$leaveOneOutObs)
      ylim <- c(min(all.ys),max(all.ys))
    } else {
      ylim <- c(min(ys,data$response),max(ys,data$response))
    }

    plot (
      xs, ys,
      xlim=xlim, ylim=ylim,
      type="l",
      main=sprintf(
        "P: theta=(%s)",
        paste(sprintf ("%.2e", P.obj$theta),collapse=",")
      ),
      xlab="covariate",
      ylab="response"
    )
    points (data$covariate, data$response, col="green")
    if (leave.one.out) {
      points (data$covariate, P.obj$leaveOneOutObs, col="red")
    }
    if ( CI.H && CI.glm ) {
      x.off <- (xmax-xmin) / (6*n)
    } else {
      x.off <- 0
    }
    if ( CI.H ) {
      # Observe this CI is based on the variance approximation in MRR note -
      # it will not coincide witht that of GLM ...
      ci <- P.obj$CI.H(data.frame(covariate=xs))
      mapply (
        function (x,low,up) {
          lines ( c(x,x) , c(low,up) , col="blue" )
        },
        xs-x.off,
        ci$lower,
        ci$upper
      )
    }
    if ( CI.glm ) {
      # GLM based CI for checking CI based on hat matrix
      pred.mu.sd <- predict (
        binomial.fit,
        data.frame(covariate=xs),
        type="response",
        se.fit = TRUE
      )
      sd.mu <- pred.mu.sd$se.fit
      mapply (
        function (x,y,sd.y) {
          lines ( c(x,x) , c(y+1.96*sd.y,y-1.96*sd.y) , col="black" )
        },
        xs+x.off,
        ys,
        sd.mu
      )
    }

    # TODO: Legend

  }

  res <- list (
    data=data,
    predictObs=predictObs,
    hatMatrixObs=hatMatrixObs,
    leaveOneOutObs=if(leave.one.out){leave.one.out.predictImpl()}else{NULL},
    theta=coef(binomial.fit),
    modelDf=length(coef(binomial.fit)),
    predict=predictImpl,
    leave.one.out.predict=leave.one.out.predictImpl,
    hatMatrix=hatMatrixImpl,
    CI.H=CI.HImpl,
    print=printImpl,
    plot=plotImpl
  )

  structure (
    res,
    class="P"
  )

}
