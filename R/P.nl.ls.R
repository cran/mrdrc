# Parametric engine for continuous data to be used in model robust regression
# data (class data.frame, object): observed data points, must hold varaibles
#   (covariate,response)
# control (class P.ctr, object): controls for this engine
# leave.one.out (boolean, atomic): it true do leave-one-out computations
P.nl.ls <- function (
  data,
  control=P.control.nl.ls(),
  leave.one.out = TRUE
) {

##  message ("\n+++ P.nl.ls: creating object")

  # Handle arguments
  if ( !inherits(control,"P.ctr")) {
    stop ("P.nl.ls: non 'P.ctr' class control")
  }
  if ( control$response.type!="continuous" ) {
    stop ("P.nl.ls: control is not for 'continuous' response")
  }
  if (
    length(
      match(
        c("covariate","response"),
        names (data))
    ) != 2
  ) {
    stop ("P.nl.ls: not all variables 'covariate','response' in data")
  }
  # TODO: Yet another test that these are indeed functions?
  fun <- get(control$fun)
  dfun <- get(control$dfun)
  ssfun <- get(control$ssfun)

  # Fit engine
  fitter <- function (subset, theta.guess) {
    if ( missing(subset) ) {
      these.data <- data
    } else {
      these.data <- data[subset,]
    }

    # Residual sum of squares to minimize
    rss <- function (theta) {
      residual <- (these.data$response - fun(these.data$covariate,theta))
      residual %*% residual
    }

    #  Derivative of residual sum of squares
    drss <- function (theta) {
      residual <- these.data$response - fun(these.data$covariate,theta)
      df <- apply (these.data["covariate"], 1, function(x) { dfun (x,theta) })
      -2 * residual %*% df
    }

    # Starting value for parameter
    if ( missing(theta.guess) ) {
      theta.guess <- ssfun (data=these.data)
    }

    # The nonlinear regression fit
    nonlinear.fit <- optim (
      theta.guess,
      rss,
      drss
    )
    if (nonlinear.fit$convergence!=0) {
      stop("P.nl.ls: unsuccesful nonlinear regression fit")
    }

    nonlinear.fit
  }

  # Fit for observed data
  nonlinear.fit <- fitter ()
  sigma2.hat <- nonlinear.fit$value / (nrow(data)-length(nonlinear.fit$par))

  # Predict engine
  predictImpl <- function (new.data) {
##    message ("P.nl.ls: predict implementation")
    predictNew <- fun (new.data$covariate, nonlinear.fit$par)
    names (predictNew) <- NULL
    predictNew
  }
  predictObs <- predictImpl(data)

  # Leave-one-out analysis
  leave.one.out.predictImpl <- function () {
##    message ("P.nl.ls: leave one out predict implementation")
    # TODO: Handle to few data points to produce fit
    n <- nrow(data)
    is <- 1:n
    leave.out <- sapply (
      1:n,
      function (index) {
        leave.out.fit <- fitter (
          subset=is[-index],
          theta.guess=nonlinear.fit$par
        )
        fun (data$covariate[index], leave.out.fit$par)
      }
    )
    names(leave.out) <- NULL
    leave.out
  }

  # Fixed matrices for hat matrix
  Fdot <- dfun (data$covariate, nonlinear.fit$par)
  hat.base <- solve(t(Fdot)%*%Fdot) %*% t(Fdot)

  # Hat matix for data points of interest
  hatMatrixImpl <- function (new.data) {
##    message ("P.nl.ls: hat matrix implementation")
    dfun (new.data$covariate, nonlinear.fit$par) %*% hat.base
  }
  hatMatrixObs <- hatMatrixImpl(data)

  # Confidence interval based on H-matrix
  CI.HImpl <- function (new.data) {
##    message ("P.nl.ls: confidence interval implementation")
    ys <- predictImpl ( new.data )
    H <- hatMatrixImpl ( new.data )
    sd.mu <-  sqrt ( diag(H %*% t(H)) * sigma2.hat )
    list (
      lower=ys-1.96*sd.mu,
      upper=ys+1.96*sd.mu
    )
  }

  printImpl <- function (P.obj, verbose = FALSE) {
##    message ("P.nl.ls: print implementation")
##    message ("theta")
##    print (P.obj$theta)
##    message ("Predicted at obs")
##    print (P.obj$predictObs)
    if (verbose) {
##      message ("Model df")
##      print (P.obj$modelDf)
##      message ("Nonlinear least squares fit")
##      print (nonlinear.fit)
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
    n=101
  ) {

    xmin <- min (data$covariate)
    xmax <- max (data$covariate)
    xs <- seq (xmin,xmax, length.out=n)
    xlim <- c(xmin,xmax)

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
    if ( CI.H ) {
      ci <- P.obj$CI.H(data.frame(covariate=xs))
      mapply (
        function (x,low,up) {
          lines ( c(x,x) , c(low,up) , col="blue" )
        },
        xs,
        ci$lower,
        ci$upper
      )
    }
    # TODO: Legend

  }

  res <- list (
    data=data,
    predictObs=predictObs,
    hatMatrixObs=hatMatrixObs,
    leaveOneOutObs=if(leave.one.out){leave.one.out.predictImpl()}else{NULL},
    theta=nonlinear.fit$par,
    modelDf=length(nonlinear.fit$par),
    sigma2.hat=sigma2.hat,
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
