# Local linear regression fit engine at xy from data points xs,ys
# xy (numeric, atomic/vector): data point to estimate at. If atomic then just
#   the covariate value of interest. If 2-vector covariate and response.
# xs (numeric, vector): observed covariates
# ys (numeric, vector): observed responses
# ws (numeric, vector): weight that applies multiplacatively with the weight
#   by the weight function
# with.hat (boolean, atomic): if true evaluate hat matrix at predictions
# weight.function (class function, object): function that maps the real
#   numbers into [0;1] putting weights on the observations that enter the
#   weighted regression that local regression boild down to
# bandwidth.type (character, ?): the type of bandwidth, refer to the code
#   to figure out what exactly how they work ...
# bandwidth (numeric, atomic): the actual bandwidth
# bandwidth.trc (boolean, atomic): if true use truncation in bandwidth
loc.reg <- function (
  xy,
  xs,
  ys,
  ws=rep(1,length(xs)),
  with.hat=FALSE,
  weight.function,
  bandwidth.type=c("fixed.width", "adaptive.portion", "adaptive.count"),
  bandwidth,
  bandwidth.trc=FALSE,
  safe=TRUE
) {

  # Book-keeping of data point of interest
  if ( length(xy)==2 ) {
    x <- xy[1]
    y <- xy[2]
  } else if ( length(xy)==1 ) {
    x <- xy[1]
    y <- NA
  } else {
    stop ("loc.reg: 1st arg should be vector of length 1 (x) or 2 (x,y)")
  }

  # Valid bandwidth?
  if ( missing(bandwidth.type) ) {
    stop ("loc.reg: missing 'bandwidth.type'")
  } else {
    bandwidth.type <- match.fun.arg (bandwidth.type) 
  }
  if ( missing(bandwidth) ) {
    stop ("loc.reg: missing 'bandwidth'")
  }

  # Really safe test on arguments passed
  if ( safe ) {
    # Check data point of interest
    if ( !(is.numeric(xy) && is.vector(xy)) ) {
      stop ("loc.reg: xy must be numeric vector")
    }

    # Check observed data points
    if ( !is.numeric(xs) || !is.numeric(ys) || (length(xs)!=length(ys)) ) {
      stop ("loc.reg: 'xs' and 'ys' must be numeric and of equal length")
    }
  
    # Check that specified weight funtion is a function
    if ( !is.function(weight.function) ) {
      cl <- match.call(expand.dots = FALSE)
      m <- match ( c("weight.function"),names(cl),0 )
      mal.fun <- cl[[m]]
      stop (
        paste (
         "loc.reg: '",
         mal.fun,
         "' not a (weight) function",
         sep=""
        )
      )
    }

    # Check bandwidth
    if ( !is.numeric(bandwidth) ) {
      stop ("loc.reg: bandwidth must be numeric vector")
    }
  }

  # Window points to be used and their distance weights
  centered.covariates <- xs-x
  dists <- abs (centered.covariates)

  # This switch clause must cover all vals in formal arg
  switch (
    bandwidth.type,
    fixed.width = {
      # message ("fixed.width")
      the.bandwidth <- bandwidth
    },
    adaptive.portion = {
      # message ("adaptive.portion")
      if (bandwidth.trc) {
        window.n <- trunc(bandwidth*length(xs))
      } else {
        window.n <- round(bandwidth*length(xs))
      }
      o.dists <- order (dists)
      the.bandwidth <- dists[o.dists][window.n]
    },
    adaptive.count = {
      # message ("adaptive.count")
      # Adaptive smoothing window by number of unique covariate-values
      unique.centered.covariates <- unique(centered.covariates)
      n <- length(unique.centered.covariates)
      abs.unique.c.c <- abs(unique.centered.covariates)
      o.covariates <- order (abs.unique.c.c)
      w.c <- bandwidth
      if (w.c < n) {
        # Somewhat arbitrary using midpoint - log-scale difficulty
        # Observe two signed unique centered covariates may have same absolute
        # value. The mean of such two is the same absoluta value and with a
        # weight function of weight zero at the max distance then gives zero
        # weight. Therefore extra if-caluse needed. On the other hand no more
        # than two signed unique centred covariates may have same abs value.
        if ( abs.unique.c.c[o.covariates][w.c] !=
             abs.unique.c.c[o.covariates][w.c+1]  ) {
          the.bandwidth <- sum(abs.unique.c.c[o.covariates][c(w.c,w.c+1)])/2.0
        } else {
          if (w.c+2<=n) {
            the.bandwidth <-sum(abs.unique.c.c[o.covariates][c(w.c,w.c+2)])/2.0
          } else {
            the.bandwidth <- 1.2 * abs.unique.c.c[o.covariates][n]
          }
          # Deprecated
          # the.bandwidth <-sum(abs.unique.c.c[o.covariates][c(w.c,w.c+1)])/2.0
        }
      } else {
        # Absolutely ad hoc 1.2
        the.bandwidth <- 1.2 * abs.unique.c.c[o.covariates][n]
      }

    }
  )

  # The local regression weights
  lr.weights <- weight.function(
    centered.covariates,
    bandwidth=the.bandwidth
  )
  lr.weights <- lr.weights*ws

  # Safeguard on effective  local regression 
  zero.weights <- (lr.weights==0)
  if ( length(xs[!zero.weights]) < 1 ) {
    stop ("loc.reg: less than 1 data point in smoothing window")
  }
# Commented out 8/2 2008 by Christian Ritz
#   else if ( length(unique(xs[!zero.weights])) == 1 ) {
#    warning (sprintf("loc.reg: data points of unique x=%f",x))
#  }

  # Estimation - do the local linear regression
  loc.data <- data.frame (x=centered.covariates, y=ys)
  loc.lm <- lm(formula=y~x, data=loc.data, weights=lr.weights)
  loc.res <- c (
    covariate=as.vector(x),
    response=as.vector(y),
    pred.loc.reg=as.vector(coef(loc.lm)["(Intercept)"])
  )
  if ( with.hat ) {
    h <- hat.point (loc.lm, 0)
    names(h) <- sapply(1:length(h), function (x) sprintf("h%d",x))
    loc.res <- c(loc.res, h)
  }

  loc.res
}
class(loc.reg) <- c(class(loc.reg),"loc.fit.engine")


# TODO: Make a local likelihood fit engine


# One-response one-covariate no-interpolation local linear regression (or
# local likelihood whenever loc.lik implemented)
# formula (class formula, object): specfies response and covariate variables
#   incl possible maps 
# data (class data.frame, object): observed data points
# subset (.): refer to subset for model.frame
# weights (numeric, vector): ancillary weights that will enter the local fit
#   engine optimization criterion. Refer to specific engine.
# new.data (class data.frame, object): data points to predict/estimate at
# predict.at.data (boolean, atomic): if true prediction/estimation at observed
#   data points
# loc.fit (character): name of local regression/likelihood fit engine to
#   use. Such engine inherits 'function' and 'loc.fit.engine'
# with.hat (boolean, atomic): if true evaluate hat matrix at predictions
# weight.function (class function, object): weight function to use in local
#   fit engine.
# bandwidth.type (character, ?): bandwidth type. Refer to specific engine.
# bandwidth (numeric, atomic): actual bandwidth. Refer to specific engine.
# bandwidth.trc (boolean, atomic): if true use truncation in bandwidth
#   computations. Refer to specific engine.
loc <- function (
  formula,
  data,
  subset,
  weights,
  new.data,
  predict.at.data=missing(new.data),
  loc.fit="loc.reg", # Should be string in order to ease debugging
  with.hat,
  weight.function,
  bandwidth.type,
  bandwidth,
  bandwidth.trc
) {

  # Check that specified local regression/likelihood engine exists
  if ( !exists(loc.fit, mode="function") || 
       !inherits(get(loc.fit),"loc.fit.engine") ) {
    # TODO: Is there more elegant way of doing this entire function body?!
    cl <- match.call(expand.dots = FALSE)
    m <- match ( c("loc.fit"),names(cl),0 )
    if ( m>0 ) { # Passed engine
      mal.fun <- cl[[m]]
    } else {     # Default engine
      mal.fun <- formals()$loc.fit
    }
    stop (
      paste (
       "loc: '",
       mal.fun,
       "' not a function that inherits 'loc.fit.engine'",
       sep=""
      )
    )
  }

  # Fetch data to use
  cl <- match.call (expand.dots = FALSE)
  m <- match (c("formula", "data", "subset", "weights"), names(cl), 0)
  cl <- cl[c(1, m)]
  cl$drop.unused.levels <- TRUE
  cl[[1]] <- as.name("model.frame")
  mf <- eval(cl, parent.frame())

  # Check that there is one respons and one covariate
  mt <- attr(mf, "terms")
  attr (mt, "intercept") <- 0
  if ( is.empty.model(mt) ) {
    stop ("loc: empty model")
  }
  if ( attr(mt,"response")==0 ) {
    stop ("loc: no model response")
  }
  Ys <- model.response(mf, "any")
#  print(Ys)
#  print(is.numeric(Ys))
#  print(is.list(Ys))
#  print(is.vector(Ys))
#  print(mode(Ys))
#  Ys <- as.vector(Ys)
 
## Commented out by Christian Ritz 27/2 2008  
#  if ( !is.vector(Ys) ) {
#    stop ("loc: response not a vector")
#  }
## Inserted instead
  Ys <- as.vector(Ys)  
  
  Xs <- model.matrix (mt, mf)
  if ( !ncol(Xs)==1 ) {
    stop ("loc: covariate not a vector (one-column matrix)")
  }

  # Fetch new.data to use
  # TODO: Think it over if this will actually do the trick!
  #       Specify under what circumstances it does work ...
  if ( !missing(new.data) ) {
    cl$formula <- formula[c(1,3)]
    cl$data <- new.data
    cl$subset <- NULL
    cl$weights <- NULL
    new.mf <- eval(cl, parent.frame())
    new.mt <- attr(new.mf, "terms")
    attr (new.mt, "intercept") <- 0
    new.Xs <- model.matrix (new.mt, new.mf)
    if ( !ncol(new.Xs)==1 ) {
      stop ("loc: covariate not a vector (one-column matrix) in new data")
    }
  }

  # Set up call that delegates to loc fit engine
  # TODO: Could not make it work with match.call ...
  cl <- call (loc.fit)
  cl$xs <- Xs[,1]
  cl$ys <- Ys
  if ( !missing(with.hat) ) {  cl$with.hat <- with.hat  }
  if ( !missing(weight.function) ) {  cl$weight.function <- weight.function  }
  if ( !missing(bandwidth.type) ) {  cl$bandwidth.type <- bandwidth.type  }
  if ( !missing(bandwidth) ) {  cl$bandwidth <- bandwidth  }
  if ( !missing(bandwidth.trc) ) {  cl$bandwidth.trc <- bandwidth.trc  }
  if ( !missing(weights) ) {  cl$ws <- mf[,"(weights)"]  }

  # Dispatch at specified data points
  if ( predict.at.data ) {
    res.obs <- apply (
                 cbind(Xs[,1],Ys),
                 1,
                 function (xy) {
                   cl$xy <- xy
                   eval(cl)
                 }
               )
  }
  predict.at.new.data <- !missing(new.data)
  if ( predict.at.new.data ) {
    res.pred <- apply (
                  new.Xs,
                  1,
                  function (x) {
                    cl$xy <- x
                    eval(cl)
                  }
                )
  }

  if ( predict.at.data && predict.at.new.data ) {
    res <- data.frame(res.obs, res.pred)
  } else if ( predict.at.data ) {
    res <- data.frame(res.obs)
  } else if ( predict.at.new.data ) {
    res <- data.frame(res.pred)
  } else {
    res <- NULL
    return (res)
  }

  # TODO: Put some proper names on the rows that are not
  response.str <- deparse(formula[[2]])
  covariate.str <- deparse(formula[[3]])
  rownames(res)[1] <- covariate.str
  rownames(res)[2] <- response.str
  rownames(res)[3] <- paste("pred",response.str,sep=".")
  colnames(res) <- NULL
  t(res)

}
