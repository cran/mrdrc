# Model roboust regression as described in TODO:CITE
# data: data frame with records on form (response,covariate[,cases])
# NP.control (class NP.ctr, object):control for non-parametric sub-model engine
# P.control (class P.ctr, object): control for parametric sub-model engine
# SP.control (class SP.ctr, object): control for this semi-pamatric engine
SP.mrr.raw <- function (
  data,
  NP.control,
  P.control,
  SP.control,
  leave.one.out = TRUE,
  logScale, general, fct, robust, compact
) {

  # Handle arguments
  if ( !(
    inherits(NP.control,"NP.ctr") &&
#    inherits (P.control,"P.ctr") &&
    inherits (SP.control,"SP.ctr")
  ) ) {
    stop ("SP.mrr.raw: one of the control parameters not of proper class")
  }

  # Check data  
  
## Modified 15/1 2008 by Christian Ritz  
#  if (
#    SP.control$response.type!=P.control$response.type ||
#    SP.control$response.type!=NP.control$response.type
#  ) {
#    stop ("SP.mrr.raw: different response types for SP/NP/P controls!")
#  }
  if ( SP.control$response.type=="binomial" ) {
    binomial.response <- TRUE
    required.var <- c("covariate","response","cases")
  } else if ( SP.control$response.type=="continuous" ) {
    binomial.response <- FALSE
    required.var <- c("covariate","response")
  } else {
    stop (
      paste(
        "SP.mrr.raw: response type '",SP.control$response.type,"' unknown",
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
        "SP.mrr.raw: not all variables '",
        paste(required.var,collapse="', '"),
        "' in data",
        sep=""
      )
    )
  }
  
  # Criterion function
  # A note on design: Dactory kind of design where calling the criterion.fac
  # will provide a criterion function object for optimization. That object
  # again may be a one or two parameter function depending on optimization
  # mode. Branching on optimization mode and the actual optimization is
  # handled seperately since the criterion function may be continuous or
  # discrete in the parameter(s) to optimize over.
  switch (

    SP.control$criterion,

    "LOGLIKELIHOOD" = {
      if ( binomial.response ) {

        criterion.fac <- function (np.arg, p.engine, verbose=FALSE) {
  
          if ( inherits(np.arg,"NP") ) {
  
            np.engine <- np.arg
            # One-parameter criterion function
            # TODO: What if shooting below 0 or over 1! Right now it goes into
            # na because a p.pred outside [0,1] will imply taking the log of a
            # negative number
            function ( mix ) {
              p.pred.leave.one.out <-
                mix*np.engine$leaveOneOutObs + (1-mix)*p.engine$leaveOneOutObs
              # TODO: Is this really sensible? We will now have a log(0) which
              # yields -Inf and in R multiplying this with non-zero number
              # preserves infinity whilst multiplying with 0 yields a NaN ->
              # in first case logL then becomes -Inf, in second case that
              # particular observation is dropped from the sum
##              p.pred.leave.one.out[p.pred.leave.one.out<0] <- 0
##              p.pred.leave.one.out[1<p.pred.leave.one.out] <- 1
## Changed by Christian Ritz 20/2 2008               
              p.pred.leave.one.out[p.pred.leave.one.out < 0.001] <- NA
              p.pred.leave.one.out[p.pred.leave.one.out > 0.999] <- NA   
## Changed by Christian Ritz 20/2 2008                         
#              successes.obs <- data$cases - data$cases*data$response
#              failures.obs <- data$cases - successes.obs
              successes.obs <- data$cases*data$response
              failures.obs <- data$cases - successes.obs
              
              # TODO: After fix which obs do actually enter this sum!? And
              #       does it make sense?
              ll <- sum(
                successes.obs*log(p.pred.leave.one.out) +
                failures.obs*log(1-p.pred.leave.one.out),
                na.rm = TRUE
              )
              if (verbose) {
                print (
                  sprintf("logl: mix=%.3f logl=%e", mix, ll)
                )
              }
## Changed by Christian Ritz 20/2 2008 
##              ll
              -ll
            }
  
          } else if ( inherits(np.arg,"NP.ctr") ) {
  
            # Two-parameter criterion function
            function ( mix.bandw ) {
              np.ctr <- np.arg
              np.ctr$bandwidth <- mix.bandw[2]
              np.ctr$optim <- FALSE
              np.engine <- NP.lr.wls (data, np.ctr)
              cri <- criterion.fac (np.engine, p.engine, verbose)
              if ( verbose ) {
                print ( sprintf ("logl: bandw=%.3f", np.ctr$bandwidth) )
              }
              cri ( mix.bandw[1] )
            }
 
          }
        }

      } else {

        paste("SP.mrr.raw: '",SP.control$criterion,"' criterion not implemented for non-binomial data",sep="")

      }

    }, # End LOGLIKELIHOOD

    "PRESS" = {
      if ( binomial.response ) {
        var.Y.inv.obs <- data$cases / (data$response*(1-data$response))
        # TODO: This is lame! Personal correspondance revealed that Nottingham
        # and Birch instead assigned a nominal value in place of infinite
        # weight
        var.Y.inv.obs[is.infinite (var.Y.inv.obs)] <- 0
      }
      criterion.fac <- function (np.arg, p.engine, verbose = FALSE) {

        if ( inherits(np.arg,"NP") ) {

          np.engine <- np.arg
          # One-parameter criterion function
          function ( mix ) {
            residual <- data$response - 
              (mix*np.engine$leaveOneOutObs + (1-mix)*p.engine$leaveOneOutObs)
            if ( binomial.response ) {
              w.rss <- sum( residual * var.Y.inv.obs * residual )
            } else {
              w.rss <- sum( residual * residual )
            }
            convex.df <- nrow(data) -
               (mix*np.engine$modelDf1  + (1-mix)*p.engine$modelDf)
            press <- w.rss / convex.df
            if (verbose) {
              print (
                sprintf(
                  "press: mix=%.3f w.rss=%e df=%5.2f press=%e",
                  mix, w.rss, convex.df, press
                )
              )
            }
            press
          }

        } else if ( inherits(np.arg,"NP.ctr") ) {

          # Two-parameter criterion function
          function ( mix.bandw ) {
            np.ctr <- np.arg
            np.ctr$bandwidth <- mix.bandw[2]
            np.ctr$optim <- FALSE
            np.engine <- NP.lr.wls (data, np.ctr)
            cri <- criterion.fac (np.engine, p.engine, verbose)
            if ( verbose ) {
              print ( sprintf ("press: bandw=%.3f", np.ctr$bandwidth) )
            }
            cri ( mix.bandw[1] )
          }

        }
      }

    }, # End PRESS

    stop (
      paste("SP.mrr.raw: '",SP.control$criterion,"' unknown criterion",sep="")
    )

  ) # End criterion switch

  # Optimization for resulting NP-P engines and mixing parameter
  NP.engine <- NULL
  P.engine <- NULL
  criterion.val <- NULL
  switch (

    SP.control$optim,

    "none" = {
      NP.engine <- NP.lr.wls (data, NP.control)
      NP.control$bandwidth <- NP.engine$bandwidth
      if ( binomial.response ) {
        P.engine <- P.glm.binomial.ml (data, P.control)
      } else {
        P.engine <- P.nl.ls (data, P.control)
      }
      criterion <- criterion.fac (NP.engine, P.engine)
      criterion.val <- criterion (SP.control$mixing)
    },

    "mixing.grid" = {
      # Seperate estimation of the non-parametric and parametric sub-models.
      # If no optim in this engine, then there is only optimization over mixing
      
## Changed 14/1 2008 by Christian Ritz      
      if (!general) {
          if ( binomial.response ) {
              P.engine <- P.glm.binomial.ml (data, P.control)
          } else {
              P.engine <- P.nl.ls (data, P.control)
          }
      } else {
          P.engine <- feedSP(data, fct = fct, logScale = logScale) 
      }


      if (compact) # Added 21/4 2008 by Christian Ritz
      {
          dataset2 <- data.frame(
          covariate = sort(unique(data$covariate)),
          response = with(data, tapply(response, covariate, mean)))
          
          if (!is.null(data$cases))
          {
              dataset2$cases <- with(data, tapply(cases, covariate, sum))
          }  
          
          NP.engine <- NP.lr.wls(dataset2, NP.control)
          NP.engine$predictObs <- NP.engine$predict(data.frame(covariate = data$covariate))
          NP.engine$hatMatrixObs <- NP.engine$hatMatrix(data.frame(covariate = data$covariate), dataSet = data)          
          NP.engine$leaveOneOutObs <- rep(NP.engine$leaveOneOutObs, with(data, tapply(covariate, covariate, length)))
                         
      } else if (robust) # Added 25/2 2008 by Christian Ritz ... not working the right way
      {
#          print(P.engine$predictObs)
          dataset2 <- data
          dataset2$response <- dataset2$response - P.engine$predictObs
#          print(dataset2$response)
          NP.engine <- NP.lr.wls(dataset2, NP.control)
      } else {
          NP.engine <- NP.lr.wls(data, NP.control)
      }
#      NP.engine <- NP.lr.wls (data, NP.control)
      NP.control$bandwidth <- NP.engine$bandwidth

      criterion <- criterion.fac (NP.engine, P.engine)
      grid.eval <- sapply (
        SP.control$mixing,
        function (mix) {
          criterion.val <- criterion (mix)
          c (mixing=mix, criterion.val=criterion.val)
        }
      )
#      print(grid.eval)
      # TODO: This is dubious. As an example with LOGLIKELIHOOD criterion
      # where one or more probabilites has been truncated to 0 or 1, the logL
      # criterion may assume -Inf over a range of the mixing parameter or have
      # discontinuities. Evaluating over the grid and naively picking mixing
      # parameter from a simple ordering of the criterion evaluated at the
      # grid point will handle such behaviour of the criterion in a most
      # probably unsatisfactory way
      best.mixing.index <- order (grid.eval["criterion.val",])[1]
      SP.control$mixing <- grid.eval["mixing", best.mixing.index]
      criterion.val <- grid.eval["criterion.val", best.mixing.index]
    },

    "two.step" = {
      # Seperate estimation of the non-parametric and parametric sub-models.
      # If no optim in this engine, then there is only optimization over mixing
      NP.engine <- NP.lr.wls (data, NP.control)
      NP.control$bandwidth <- NP.engine$bandwidth
      if ( binomial.response ) {
        P.engine <- P.glm.binomial.ml (data, P.control)
      } else {
        P.engine <- P.nl.ls (data, P.control)
      }
      criterion <- criterion.fac (NP.engine, P.engine)
      mix.min <- min (SP.control$mixing)
      mix.max <- max (SP.control$mixing)
      # TODO: Will this actually work - jumping number of observations in
      #       optimization criterion cause non-continuous func!
      mix.optim <- optim (
        (mix.min+mix.max) / 2,
        criterion,
        method="L-BFGS-B",
        lower=mix.min,
        upper=mix.max
      )
      if ( mix.optim$convergence!=0 ) {
        stop ("SP.mrr.raw: failed optimization over mixing")
      } else {
        SP.control$mixing <- mix.optim$par
        criterion.val <- mix.optim$val
      }
    },

    "joint" = {
      if ( binomial.response ) {
        P.engine <- P.glm.binomial.ml (data, P.control)
      } else {
        P.engine <- P.nl.ls (data, P.control)
      }
      criterion <- criterion.fac (NP.control, P.engine)
      mix.min <- min (SP.control$mixing)
      mix.max <- max (SP.control$mixing)
      switch (
        NP.control$bandwidth.type,
        "fixed.width"={
          # Criterion continuous in both mixing and bandwidth
           bw.min <- min (NP.control$bandwidth)
           bw.max <- max (NP.control$bandwidth)
           mix.bw.optim <- optim (
             c ( (mix.min+mix.max)/2, (bw.min+bw.max)/2 ),
             criterion,
             method="L-BFGS-B",
             lower=c(mix.min,bw.min),
             upper=c(mix.max,bw.max)
           )
          if ( mix.bw.optim$convergence!=0 ) {
            stop ("SP.mrr.raw: failed optimization over mixing-bandwidth")
          } else {
            NP.control$optim <- FALSE
            NP.control$bandwidth <- mix.bw.optim$par[2]
            NP.engine <- NP.lr.wls (data, NP.control)
            SP.control$mixing <- mix.bw.optim$par[1]
            criterion.val <- mix.bw.optim$val
          }
        },
        "adaptive.portion"={
          # Criterion continuous in mixing, discrete in bandwidth
          stop ("SP.mrr.raw: joint optimzation for adaptive.portion not implemented")
        },
        "adaptive.count"={
          # Criterion continuous in mixing, discrete in bandwidth
          stop ("SP.mrr.raw: joint optimzation for adaptive.count not implemented")
        },
        stop (
          paste("SP.mrr.raw: '",SP.control$optim,"' unknown optim mode",sep="")
        )
      ) # End bandiwth.type switch
    },

    stop (
      paste("SP.mrr.raw: '",SP.control$optim,"' unknown optim mode",sep="")
    )

  ) # End optim switch

##  message ("SP.mrr.raw: SP.control$mixing")
##  print (SP.control$mixing)
##  message ("SP.mrr.raw: NP.control$bandwidth")
##  print (NP.control$bandwidth)
##  message ("SP.mrr.raw: NP.engine$bandwidth")
##  print (NP.engine$bandwidth)

  # At observed data points
  predictObs <- SP.control$mixing*NP.engine$predictObs +
      (1-SP.control$mixing)*P.engine$predictObs
  hatMatrixObs <- SP.control$mixing*NP.engine$hatMatrixObs +
      (1-SP.control$mixing)*P.engine$hatMatrixObs
  modelDf <- SP.control$mixing*NP.engine$modelDf1 +
      (1-SP.control$mixing)*P.engine$modelDf

  # Predict engine
  # new.data (class data.frame, object): data frame with variable covariate
  predictImpl <- function (new.data) {
##    message ("SP.mrr.raw: predict implementation")
    SP.control$mixing*predict(NP.engine,new.data) +
      (1-SP.control$mixing)*predict(P.engine,new.data)
  }

  # Leave-one-out analysis
  leave.one.out.predictImpl <- function () {
##    message ("SP.mrr.raw: leave one out predict implementation")
    SP.control$mixing*NP.engine$leaveOneOutObs +
      (1-SP.control$mixing)*P.engine$leaveOneOutObs
  }

  # Hat matrix
  # new.data (class data.frame, object): data frame with variable covariate
  hatMatrixImpl <- function(new.data) 
  {
##    message ("SP.mrr.raw: hat matrix implementation")
# Commented out by Christian Ritz 21/4 2008
#    SP.control$mixing*hatMatrix(NP.engine, new.data, dataSet = data) +
#      (1-SP.control$mixing)*hatMatrix(P.engine,new.data)
      (SP.control$mixing*NP.engine$hatMatrix(new.data, dataSet = data) +
      (1-SP.control$mixing)*P.engine$hatMatrix(new.data))
  }

  # Variance for CI
  switch (
    SP.control$response.type,
    # TODO: This is just an idea ... but is it OK?
    "continuous" = {
      sigma2.hat <- sum((data$response - predictObs)^2) / (nrow(data)-modelDf)
    },
    # TODO: Use predicted proportions instead?
    "binomial" = {
      var.Y.obs <- (data$response*(1-data$response)) / data$cases
    },
    stop ("SP.mrr.raw: unknown response type for variance")
  )

  # CI for ED based on H based CI
  # alpha (numeric, atomic): percentage
  ED.HImpl <- function (alpha, reference, level) {
##    message ("SP.mrr.raw: ED_alpha estimation")
    lower.covariate <- min(data$covariate)
    upper.covariate <- max(data$covariate)
    response.lower.covariate <- 
      predictImpl (data.frame(covariate=lower.covariate))
    response.upper.covariate <- 
      predictImpl (data.frame(covariate=upper.covariate))
    is.increasing.curve <- response.upper.covariate > response.lower.covariate
    lower.response <- min (response.lower.covariate,response.upper.covariate)
    upper.response <- max (response.lower.covariate,response.upper.covariate)

    if (is.increasing.curve) {propor <- alpha/100} else {propor <- 1 - alpha/100}
    switch (
      SP.control$response.type,
      "continuous" = {
        # TODO: Is this appropriate target response for inversion?
## Changed 15/1 2008 by Christian Ritz        
#        target.response <- response.lower.covariate + 
#                (alpha/100)*(response.lower.covariate-response.upper.covariate)

         if (!is.null(reference)) {lrVal <- reference} else {lrVal <- lower.response} 
         target.response <- lrVal + propor * (upper.response - lrVal)
      },
      "binomial" = {
        if (!is.null(reference)) {lrVal <- reference} else {lrVal <- 0}     
        target.response <- lrVal + propor * (1 - lrVal)
      },
      stop ("SP.mrr.raw: unknown response type for ED")
    )

#    print(target.response)
#    print(response.lower.covariate)
#    print(response.upper.covariate)

    # Point estimate of ed.alpha
    if ( target.response < lower.response ) {
      # TODO: What really to do?
      ed.alpha.point.estimate <- NaN
    } else if (target.response > upper.response) {
      # TODO: What really to do?
      ed.alpha.point.estimate <- NaN
    } else {
      ed.alpha.point.estimate <- uniroot (
        function (x) {
          predictImpl(data.frame(covariate=x))-target.response
        },
        c(lower.covariate, upper.covariate)
      )$root
    }

    # Confidence interval based on H-matrix
    # new.data (class data.frame, object): data frame with variable covariate
    CI.HImpl <- function (new.data) {
##    message ("SP.mrr.raw: confidence interval implementation")
    ys <- predictImpl ( new.data )
    H <- hatMatrixImpl ( new.data )
    switch (
      SP.control$response.type,
      "continuous" = {
        sd.mu <-  sqrt ( diag(H %*% t(H)) * sigma2.hat )
      },
      "binomial" = {
        sd.mu <- sqrt ( diag ( H %*% diag(var.Y.obs) %*% t(H) ) )
#        print(c(ys, sd.mu))
      },
      stop ("SP.mrr.raw: unknown response type for CI")
      )
      CIquan <- qnorm(1 - (1 - level)/2)
      list (
        lower = ys - CIquan * sd.mu,
        upper = ys + CIquan * sd.mu
      )
    }    
        
    # Working on intersection with lower CI.H for mu  
    lower.CI.lower.covariate <-
      CI.HImpl (data.frame(covariate=lower.covariate))$lower
    lower.CI.upper.covariate <-
      CI.HImpl (data.frame(covariate=upper.covariate))$lower
    lower.lower.CI <- min (lower.CI.lower.covariate,lower.CI.upper.covariate)
    upper.lower.CI <- max (lower.CI.lower.covariate,lower.CI.upper.covariate)
    if ( target.response < lower.lower.CI ) {
      ed.CI.1 <- NaN
    } else if (target.response > upper.lower.CI) {
      ed.CI.1 <- NaN
    } else {
      ed.CI.1 <- uniroot (
        function (x) {
          CI.HImpl(data.frame(covariate=x))$lower - target.response
        },
        c(lower.covariate, upper.covariate)
      )$root
    }

    # Working on intersection with upper CI.H for mu
    upper.CI.lower.covariate <- 
      CI.HImpl (data.frame(covariate=lower.covariate))$upper
    upper.CI.upper.covariate <-
      CI.HImpl (data.frame(covariate=upper.covariate))$upper
    lower.upper.CI <- min (upper.CI.lower.covariate,upper.CI.upper.covariate)
    upper.upper.CI <- max (upper.CI.lower.covariate,upper.CI.upper.covariate)
    if ( target.response < lower.upper.CI ) {
      ed.CI.2 <- NaN
    } else if (target.response > upper.upper.CI) {
      ed.CI.2 <- NaN
    } else {
      ed.CI.2 <- uniroot (
        function (x) {
          CI.HImpl(data.frame(covariate=x))$upper - target.response
        },
        c(lower.covariate, upper.covariate)
      )$root
    }

    # Assemble CI for ED - truncate at range of observed covariate!
    # TODO: Is this appropriate truncation
    # TODO: In case ed.alpha.point.estimate is NaN these truncations are not
    #       necessarily meaningful
    if (is.increasing.curve) {
      if ( is.nan(ed.CI.1) ) {
        ed.alpha.upper <- upper.covariate
      } else {
        ed.alpha.upper <- ed.CI.1
      }
      if ( is.nan(ed.CI.2) ) {
        ed.alpha.lower <- lower.covariate
      } else {
        ed.alpha.lower <- ed.CI.2
      }
    } else {
      if ( is.nan(ed.CI.1) ) {
        ed.alpha.lower <- lower.covariate
      } else {
        ed.alpha.lower <- ed.CI.1
      }
      if ( is.nan(ed.CI.2) ) {
        ed.alpha.upper <- upper.covariate
      } else {
        ed.alpha.upper <- ed.CI.2
      }
    }
    if (is.nan(ed.alpha.point.estimate))
    {
        ed.alpha.point.estimate <- NA
        ed.alpha.lower <- NA
        ed.alpha.upper <- NA
    }
    list (
      ed.alpha=ed.alpha.point.estimate,
      ed.alpha.lower=ed.alpha.lower,
      ed.alpha.upper=ed.alpha.upper
    )
  }

  # Print the supplied SP object
  # SP.obj (class SP, object): SP object
  # verbose (boolean, atomic): verbose output or not
  printImpl <- function (SP.obj, verbose = FALSE) {
##    message ("SP.mrr.raw: print implementation")
##    message ("mixing")
##    print (SP.obj$mixing)
##    message ("Model df")
##    print (SP.obj$modelDf)
##    message ("Criterion value")
##    print (SP.obj$criterion.val)
    if (verbose) {
##      message ("NP")
##      print (SP.obj$NP)
##      message ("P")
##      print (SP.obj$P)
    }
  }

  # Plot the supplied SP object
  # SP.obj (class SP, object): SP object
  # leave.one.out (boolean, atomic): show leave-on-out predictions
  # n (numeric, atomic): number of equidistant covariate values to evalue at
  # with.SP (boolean, atomic): do sub-plot of SP
  # with.NP (boolean, atomic): do sub-plot of NP
  # with.P (boolean, atomic): do sub-plot of P
  # with.ED (boolean, atomic): do ED
  # alpha (numeric, atomic/vector): percentages
  plotImpl <- function (
    SP.obj,
    leave.one.out = FALSE,
    CI.H = FALSE,
    n = 101,
    with.SP = TRUE,
    with.NP = TRUE,
    with.P = TRUE,
    with.ED = FALSE,
    alpha = 50,
    main = "", xlab = "", ylab = "", logPred = FALSE, lscale, col = "blue"
  ) {
##    message ("SP.mrr.raw: plot implementation")
    xmin <- min(data$covariate)
    xmax <- max(data$covariate)
    xs <- seq(xmin, xmax, length.out = n)

    ## Inserted by Christian Ritz 5/5 2008
    ##  to avoid problems with zeros
    if (identical(lscale, "x"))
    {
        xsInd <- (abs(xs) < 1e-10)
        conLevel <- 0.1 * min(xs[!xsInd])
        xs[xsInd] <- conLevel
    }
#    xlim <- c(xmin, xmax)
    xlim <- range(xs)

    # Sub plot
    subplot <- function(M.obj) {

      ys <- predict (M.obj, data.frame(covariate = xs)) # TODO: Does this work?
      if (leave.one.out) {
        all.ys <- c(ys,data$response,M.obj$leaveOneOutObs)
        ylim <- c(min(all.ys),max(all.ys))
      } else {
        ylim <- c(min(ys,data$response),max(ys,data$response))
      }
## Commented out by Christian Ritz 14/2 2008  
#      if ( inherits(M.obj,"SP") ) {
#        main <- sprintf ("SP: mixing=%.3f", M.obj$mixing)
#      } else if ( inherits(M.obj,"NP") ) {
#        main <- sprintf("NP: %s, bandwidth=%.2e",
#                        M.obj$bandwidth.type, M.obj$bandwidth)
#      } else if ( inherits(M.obj,"P") ) {
#        main <- sprintf(
#          "P: theta=(%s)",
#          paste(sprintf ("%.2e", M.obj$theta),collapse=",")
#        )
#      }

      ## Getting back from log scale
      if (logPred) 
      {
          xs <- exp(xs)
          xlim <- exp(xlim)
          xVar <- exp(data$covariate)    
      } else {
          xVar <- data$covariate    
      }
      if (identical(lscale, "x"))
      {
        xvInd <- (abs(xVar) < 1e-10)
        xVar[xvInd] <- conLevel
      }
      
      ## Plotting
      plot (
        xs, ys,
        xlim = xlim, ylim = ylim,
        type = "l",
        main = main,
        xlab = xlab,  # "covariate",
        ylab = ylab,   # "response"
        log = lscale
      )
      if ( CI.H ) {
        ci <- M.obj$CI.H (data.frame(covariate=xs))
        mapply (
          function (x,low,up) {
            lines ( c(x,x) , c(low,up) , col="grey" )
          },
          xs,
          ci$lower,
          ci$upper
        )
      }
      if (leave.one.out) {
        points (data$covariate, M.obj$leaveOneOutObs, col="red")
      }
      if ( with.ED && inherits(M.obj,"SP") ) {
        sapply (
          alpha,
          function (a) {
            ed <- ED.H (M.obj, a)
            lines (
              c(ed$ed.alpha.lower, ed$ed.alpha.upper),
              c(a/100, a/100),
              col="blue"
            )
          }
        )
      }
      points (xVar, data$response, col = col)
      # TODO: Legend
    }

#    par.old <- par (mfcol=c(1,with.SP+with.NP+with.P))

    if (with.SP) { subplot (SP.obj)  }
    if (with.NP) { subplot (SP.obj$NP)  }
    if (with.P) { subplot (SP.obj$P)  }

#    par (par.old)

  }

  res <- list(
    data=data,
    predictObs=predictObs,
    hatMatrixObs=hatMatrixObs,
    leaveOneOutObs=if(leave.one.out){leave.one.out.predictImpl()}else{NULL},
    modelDf=modelDf,
    mixing=SP.control$mixing,
    NP=NP.engine,
    P=P.engine,
    criterion.val=criterion.val,
    predict=predictImpl,
    leave.one.out.predict=leave.one.out.predictImpl,
    hatMatrix=hatMatrixImpl,
#    CI.H=CI.HImpl,
    ED.H=ED.HImpl,
    print=printImpl,
    plot=plotImpl
  )

  structure (
    res,
    class="SP"
  )

}
