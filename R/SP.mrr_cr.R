# R-type of interface with formula to specify data for model roboust
#   regression. The model roboust estimation engine itself is SP.mrr.raw.
# formula (class formula, object): specifies response and covariate incl
#   possible maps. For binomial data, the response can be
#     1) portion of successes
#     2) (successes, failures)
#     3) factor with lowest level failure and other levels success
# data (class data.frame, object): observed data points
# subset (.): refer to subset for model.frame
# cases (numeric, vector): for binomial response, number of cases (trials)
#   when response given as portion of successes
# experiment (numeric, factor): with binomial data with response coded a
#   factor that codes for experiment (repetition)
# NP.control (class NP.ctr, object):control for non-parametric sub-model engine
# P.control (class P.ctr, object): control for parametric sub-model engine
# SP.control (class SP.ctr, object): control for the semi-pamatric engine
#   SP.mrr.raw 
SP.mrr <- function (
  formula,
  data,
  subset,
  cases,
  experiment,
  NP.control,
  P.control,
  SP.control,
  logScale = TRUE,  
  general = TRUE,
  fct,
  robust = FALSE,
  compact = TRUE
) {

  # Handle arguments
  if ( !(
    inherits(NP.control,"NP.ctr") &&
#    inherits (P.control,"P.ctr") &&
    inherits (SP.control,"SP.ctr")
  ) ) {
    stop ("SP.mrr: one of the control parameters not of proper class")
  }

  # Decide if binomial or continuous response
## Modified 15/1 2008 by Christian Ritz
#  if ( NP.control$response.type!=P.control$response.type ) {
#    stop ("SP.mrr: different response types for NP and P models")
#  }

## Changed 21/4 by Christian Ritz
  if (SP.control$response.type == "binomial") 
  {
      binomial.response <- TRUE
  } else if (SP.control$response.type == "continuous") 
  {
      binomial.response <- FALSE
  } else {
      stop ("SP.mrr: unknown P.control$response.type")
  }

  # Fetch data for SP.mrr.raw and bring to format suitable for SP.mrr.raw.
  # That format is (row-wise)
  #   (covariate,response[,cases])
  mf <- match.call (expand.dots = FALSE)
  m <- match (c("formula","data","subset","cases","experiment"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE       # Handle response given as factor
  mf[[1]] <- as.name ("model.frame")
  mf <- eval (mf, parent.frame())

  # Y
  Y <- model.response(mf, "any")
  if ( is.null(Y) ) {
    stop ("SP.mrr: no model response")
  }

  # Cases
  cases <- model.extract (mf,"cases")

  # Experiment
  experiment <- model.extract (mf,"experiment")

  # X
  mt <- attr(mf, "terms")
  attr(mt,"intercept") <- 0           # Force no intercept
  if ( is.empty.model(mt) ) {
    stop ("SP.mrr: empty model")
  }
  X <- model.matrix(mt, mf)
  if ( ncol(X)!=1 ) {
    stop ("SP.mrr: more than one covariate variable")
  } else {
    if ( is.factor(Y) ) {
      if ( !binomial.response ) {
        stop ("SP.mrr: factor response for non-binomial response")
      }
      if ( is.factor(experiment) ) {
        # Response factor format: lowest level failure, other levels success
        # Experiment factor codes repeated experiments. This clause is
        # somewhat tricky, but it produces a set of logical vectors to subset
        # the observed data into unique (covariate,experiment) combinations so
        # that there can be aggregated accordingly. Also these actual unique
        # (covariate,experiment)
        # Unique (covariate,experiment) combinations
        experiment.covariate <- data.frame(covariate=X[,1],experiment)
        u.e.c <- unique(experiment.covariate)
        # Apply on u.e.c caused unwanted type converstion on covariate
        # into factor as experiment. Thus apply on integer index.
        aggregate.subsets <- sapply (
          1:nrow(u.e.c),
          function (i) {
            (experiment.covariate$experiment==u.e.c[i,"experiment"]) &
            (experiment.covariate$covariate==u.e.c[i,"covariate"])
          }
        )
        raw.data <- data.frame (covariate=u.e.c$covariate)
      } else {
        # Unique covariate
        u.c <- data.frame(covariate=unique(X[,1]))
        aggregate.subsets <- apply (
          u.c,
          1,
          function (x) {
            X[,1]==x
          }
        )
        raw.data <- u.c
      }
    } else {
      # Any other response format
      raw.data <- data.frame (covariate=X[,1])
    }
  }

  # Y ctd
  if ( binomial.response ) {
    if ( is.factor(Y) ) {
      # Factor formatted response
      if ( nlevels(Y)<2 ) {
        stop ("SP.mrr(binomial): less than two level response")
      }
      failure.level <- levels(Y)[1]
      failure.indices <- (Y==failure.level)
      response.cases <- apply (
        aggregate.subsets,
        2,
        function (aggregate.indices) {
          these.failures <- (aggregate.indices & failure.indices)
          these.successes <- (aggregate.indices & !failure.indices)
          if ( is.null(cases) ) {
            f <- sum (these.failures)
            s <- sum (these.successes)
          } else {
            f <- sum(cases[these.failures])
            s <- sum (cases[these.successes])
          }
          n <- f+s
          p <- s / n            
          c(response=p,cases=n)
        }
      )
      raw.data$response <- response.cases["response",]
      raw.data$cases <- response.cases["cases",]
    } else if ( is.vector(Y) ) {
      # Response univariate vector
      raw.data$response <- Y
    } else if ( is.array(Y) && ncol(Y)==2 ) {
      # Response pairs of (success,failure) counts
      raw.data$response <- Y[,1]/(Y[,1]+Y[,2])
    } else {
      # No can do
      stop ("SP.mrr(binomial): more than two response variable")
    }
  } else {
#    if ( !is.null(dim(Y)) ) {
#      stop ("SP.mrr(continuous): more than one response variable")
#    } else {
        raw.data$response <- as.vector(Y)  # Modified by Christian Ritz 21/4 2008
#    }
  }

  # Cases ctd
  if ( !is.null(cases) ) {
    if ( is.factor(Y) ) {
      # Handled in clause for Y
    } else {
      # Any other response format
      raw.data$cases=cases
    }
  }

  # TODO: Any other final checks on that data are OK?

  SP.mrr.raw (
    raw.data,
    NP.control,
    P.control,
    SP.control,
    logScale = logScale,    
    general = general,
    fct = fct,
    robust = robust,
    compact = compact
  )

}
