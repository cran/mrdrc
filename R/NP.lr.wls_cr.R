# Non-parametric engine  to be used in model robust regression
# data (class data.frame, object): observed data points, must hold varaibles
#   (covariate,response[,cases])
# control (class P.ctr, object): controls for this engine
# leave.one.out (boolean, atomic): it true do leave-one-out computations
NP.lr.wls <- function (
  data,
  control,
  leave.one.out = TRUE
) {

##  message ("\n+++ NP.lr.wls: creating object")

  # Handle arguments
  if ( !inherits(control,"NP.ctr")) {
    stop ("NP.lr.wls: non 'NP.ctr' class control")
  }
  if ( control$response.type=="binomial" ) {
    required.var <- c("covariate","response","cases")
  } else if ( control$response.type=="continuous" ) {
    required.var <- c("covariate","response")
  } else {
    stop (
      paste(
        "NP.lr.wls: response type '",control$response.type,"' unknown",
        sep=""
      )
    )
  }
  if (
    length(
      match(
        required.var,
        names (data))
    ) != length(required.var)
  ) {
    stop (
      paste (
        "NP.lr.wls: not all variables '",
        paste(required.var,collapse="', '"),
        "' in data",
        sep=""
      )
    )
  }

  # Leave-one-out analysis
  leave.one.out.predictImpl <- function (bandwidth=control$bandwidth) {
##    message ("NP.lr.wls: leave one out predict implementation")
    n <- nrow(data)
    is <- 1:n
    sapply (
      1:n,
      function (index) {
        leave.out.fit <-  loc (
          formula = response ~ covariate,
          data=data,
          subset=is[-index],
          new.data=data.frame(covariate=data$covariate[index]),
          predict.at.data=FALSE,
          with.hat=FALSE,
          weight.function=control$weight.function,
          bandwidth.type=control$bandwidth.type,
          bandwidth=bandwidth
        )
        leave.out.fit[,3]
      }
    )
  }

  # Hat matix for data points of interest
  hatMatrixImpl <- function (new.data, bandwidth = control$bandwidth, dataSet = NULL) 
  {
##    message ("NP.lr.wls: hat matrix implementation")
      if (is.null(dataSet)) 
      {
          dataSet <- data
      }
      
      lr.wls.fit <- loc(
      formula = response ~ covariate,
      data = dataSet,
      new.data = new.data,
      predict.at.data = FALSE,
      with.hat = TRUE,
      weight.function = control$weight.function,
      bandwidth.type = control$bandwidth.type,
      bandwidth = bandwidth)

      lr.wls.fit[,4:ncol(lr.wls.fit)]
  }

#  print(control$criterion)
  # Fit engine
  if ( control$optim!="none" ) {
    # Set up criterion to optimize over
    switch (

      control$criterion,

      "LOGLIKELIHOOD"={
        if ( control$response.type=="binomial" ) {

          criterion <- function (bandw, verbose=FALSE) {
            p.pred.leave.one.out <- leave.one.out.predictImpl (bandw)
            # TODO: Is this realy sensible?
## Changed by Christian Ritz 19/2 2008            
#            p.pred.leave.one.out[p.pred.leave.one.out<0] <- 0
            p.pred.leave.one.out[p.pred.leave.one.out < 0.001] <- NA  # 0.01
#            p.pred.leave.one.out[1<p.pred.leave.one.out] <- 1
            p.pred.leave.one.out[p.pred.leave.one.out > 0.999] <- NA  # 0.99
## Changed by Christian Ritz 20/2 2008
#            successes.obs <- data$cases - data$cases*data$response
#            failures.obs <- data$cases - successes.obs
            successes.obs <- data$cases * data$response
            failures.obs <- data$cases - successes.obs

            # TODO: After fix which obs do actually enter this sum!? And
            #       does it make sense?
            ll <- sum(
              successes.obs*log(p.pred.leave.one.out) +
              failures.obs*log(1-p.pred.leave.one.out),
              na.rm = TRUE
            )
            if (verbose) {
#              print (
#                sprintf("logl: bandw=%e logl=%e", bandw, ll)
#              )
            }
## Changed by Christian Ritz 20/2 2008                       
##            ll
            -ll
          }

        } else {

          paste("NP.lr.wls: '",control$criterion,"' criterion not implemented for non-binomial data",sep="")

        }

      }, # End LOGLIKELIHOOD

      "PRESS"={
        if ( control$response.type=="binomial" ) {
          var.Y.inv.obs <- data$cases / (data$response*(1-data$response))
          var.Y.inv.obs[is.infinite (var.Y.inv.obs)] <- 0
        }
        # Criterion function to optimize
        criterion <- function (bandw, verbose=FALSE) {
          pred <- leave.one.out.predictImpl (bandw)
          H <- hatMatrixImpl (data, bandw)
          residual <- data$response - pred
          if (control$response.type == "binomial") {
            n <- sum (residual*var.Y.inv.obs*residual)
          } else {
            n <- sum (residual*residual)
          }
          d <- nrow(data) - sum(diag(H))
          
          if (verbose) {
#            print ("-----------------------------------------------------")
#            str <- sprintf ("bw=%.3f, PRESS = % e/% e = % e",bandw,n,d,n/d)
##            message (str)
##            print ( rbind(data$response,pred,residual))
            #str <- sprintf ("%.3f % e % e % e",bandw,n,d,n/d)
            #cat (str, file="press.log", fill=TRUE, append=TRUE)
          }
          n/d
        }
      }, # End PRESS
      # TODO: If adding another criterion function here, be careful that the
      #       following switch is still valid
      stop (
        paste (
         "NP.lr.wls: unknown criterion '",
         control$criterion,
         "' for bandwidth optimization",
         sep=""
        )
      )
    ) # End criterion switch

    # Do the optimization
    switch (

      control$bandwidth.type,

      "fixed.width" = {
        # criterion continuous in bandwidth
        bw.min <- min (control$bandwidth)
        bw.max <- max (control$bandwidth)
        if ( bw.min==bw.max ) {
          control$bandwidth <- bw.min
        } else {
          switch (
            control$optim,
            "bandwidth.grid" = {
              grid.eval <- sapply (
                control$bandwidth,
                function (bandw) {
#                  criterion.val <- criterion (bandw)
                  criterion.val <- try(criterion(bandw), silent = TRUE)
                  if (inherits(criterion.val, "try-error"))
                  {
                      criterion.val <- NA
                  }
                  # criterion.val <- criterion (bandw, verbose=T)
                  c(bandw=bandw, criterion.val=criterion.val)
                }
              )
#              print(grid.eval)
              best.bandw.index <- order(grid.eval["criterion.val",])[1]
              control$bandwidth <- grid.eval["bandw", best.bandw.index]
            },
            "one.step" = {
              bw.optim <- optim (
                (bw.min+bw.max)/2,
                criterion,
                # function (x) { criterion(x,verbose=T) },
                method="L-BFGS-B",
                lower=bw.min,
                upper=bw.max
              )
              if ( bw.optim$convergence!=0 ) {
                stop ("NP.lr.wls: failed optimization over bandwidth")
              } else {
                control$bandwidth <- bw.optim$par
              }
           },
           stop ("NP.lr.wls: '",control$optim,"' optim not implemented for \"fixed.bandwidth\"",sep="")
          ) # End optim switch
        }
      },

      "adaptive.portion" = {
        # Criterion discrete in portion - grid evaluation only of relevance
        bw.min <- min (control$bandwidth)
        bw.max <- max (control$bandwidth)
        if ( bw.min==bw.max ) {
          control$bandwidth <- bw.min
        } else {
          bw.of.interest <- 0:nrow(data) / nrow(data)
          subset.of.interest<-(bw.min<=bw.of.interest)&(bw.of.interest<=bw.max)
          bw.of.interest <- bw.of.interest[subset.of.interest]
          cris <- sapply (
            bw.of.interest,
            function (bw) {
              c(bw,criterion(bw))
              # c(bw,criterion(bw,verbose=T))
            }
          )
          o.cris <- order(cris[2,])
          control$bandwidth <- cris[1,o.cris[1]]
        }
      },

      "adaptive.count" = {
        # Criterion discrete in count - grid evaluation only of relevance
        bw.min <- min (control$bandwidth)
        bw.max <- max (control$bandwidth)
        if ( bw.min==bw.max ) {
          control$bandwidth <- bw.min
        } else {
          bw.of.interest <- 1:length(unique(data$covariate))
          subset.of.interest<-(bw.min<=bw.of.interest)&(bw.of.interest<=bw.max)
          bw.of.interest <- bw.of.interest[subset.of.interest]
          cris <- sapply (
            bw.of.interest,
            function (bw) {
              c(bw,criterion(bw))
              # c(bw,criterion(bw,verbose=T))
            }
          )
          o.cris <- order(cris[2,])
          control$bandwidth <- cris[1,o.cris[1]]
        }
      },

      stop (
        paste (
         "NP.lr.wls: bandwidth type '",
         control$bandwidth.type,
         "' for bandwidth optimization",
         sep=""
        )
      )

    ) # End bandwidth.type switch

  }

  # Not until now are we sure on bandwidth -> these calls can not be situated
  # right after the hatMatrix function
  hatMatrixObs <- hatMatrixImpl(data)
  modelDf1 <- sum(diag(hatMatrixObs))
  modelDf2 <- sum(diag(t(hatMatrixObs)%*%hatMatrixObs))

  # Predict engine
  predictImpl <- function (new.data) {
##    message ("NP.lr.wls: predict implementation")
    lr.wls.fit <- loc (
      formula = response ~ covariate,
      data=data,
      new.data=new.data,
      predict.at.data=FALSE,
      with.hat=FALSE,
      weight.function=control$weight.function,
      bandwidth.type=control$bandwidth.type,
      bandwidth=control$bandwidth
    )
    lr.wls.fit[,3]
  }
  predictObs <- predictImpl(data)

  # Variance for CI
  switch (
    control$response.type,
    "continuous" = {
      # Refer to Loader (2.18)
      sigma2.hat <- sum((data$response - predictObs)^2) /
                    (nrow(data)-2*modelDf1 + modelDf2)
    },
    "binomial" = {
      var.Y.obs <- (data$response*(1-data$response)) / data$cases
    },
    stop ("NP.lr.wls: unknown data type for variance")
  )

  # Confidence interval based on H-matrix
  CI.HImpl <- function (new.data) {
##    message ("NP.lr.wls: confidence interval implementation")
    ys <- predictImpl ( new.data )
    H <- hatMatrixImpl ( new.data )
    switch (
      control$response.type,
      "continuous" = {
        sd.mu <-  sqrt ( diag(H %*% t(H)) * sigma2.hat )
      },
      "binomial" = {
        sd.mu <- sqrt ( diag ( H %*% diag(var.Y.obs) %*% t(H) ) )
      },
      stop ("NP.lr.wls: unknown data type for CI")
    )
    list (
      lower=ys-1.96*sd.mu,
      upper=ys+1.96*sd.mu
    )
  }

  printImpl <- function (NP.obj, verbose=FALSE) {
##    message ("NP.lr.wls: print implementation")
##    message ("Predicted at obs")
##    print (NP.obj$predictObs)
##    message ("bandwidth.type")
##    print (NP.obj$bandwidth.type)
##    message ("bandwidth")
##    print (NP.obj$bandwidth)
    if (verbose) {
##      message ("Model df1")
##      print (NP.obj$modelDf1)
##      message ("Model df2")
##      print (NP.obj$modelDf2)
##      message ("Hat matrix at obs")
##      print (NP.obj$hatMatrixObs)
##      message ("Leave out one predictions at obs")
##      print (NP.obj$leaveOneOutObs)
    }
  }

  plotImpl <- function (
    NP.obj,
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
      all.ys <- c(ys,data$response,NP.obj$leaveOneOutObs)
      ylim <- c(min(all.ys),max(all.ys))
    } else {
      ylim <- c(min(ys,data$response),max(ys,data$response))
    }

    plot (
      xs, ys,
      xlim=xlim, ylim=ylim,
      type="l",
      main=sprintf("NP: %s, bandwidth=%.2e", NP.obj$bandwidth.type, NP.obj$bandwidth),
      xlab="covariate",
      ylab="response"
    )
    points (data$covariate, data$response, col="green")
    if (leave.one.out) {
      points (data$covariate, NP.obj$leaveOneOutObs, col="red")
    }
    if ( CI.H ) {
      ci <- NP.obj$CI.H (data.frame(covariate=xs))
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
    modelDf1=modelDf1,
    modelDf2=modelDf2,
    predict=predictImpl,
    leave.one.out.predict=leave.one.out.predictImpl,
    hatMatrix=hatMatrixImpl,
    CI.H=CI.HImpl,
    print=printImpl,
    plot=plotImpl,
    bandwidth.type=control$bandwidth.type,
    bandwidth=control$bandwidth
  )

  structure (
    res,
    class="NP"
  )

}
