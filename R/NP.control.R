# Local regression control parameters
# response.type (cheracter, ?): response type
# criterion (cheracter, ?): criterion to optimize over wrt bandwidth
# optim (cheracter, ?): the optimization mode
# weight.function (cheracter, ?): the weight function for local regression
# bandwidth.type (cheracter, ?): the bandwidth type
# bandwidth (numeric, atomic/vector): the actual bandwidth, use a range or
#   grid if appropriate with the chocen optimization mode
NP.control.lr.wls <- function (
  response.type=c("binomial", "continuous"),
  criterion=c("LOGLIKELIHOOD", "PRESS"),
  optim=c("none", "bandwidth.grid", "one.step"),
  weight.function=c("gaussian.weight","tricubic.weight"),
  bandwidth.type=c("fixed.width", "adaptive.portion", "adaptive.count"),
  bandwidth=c(0.05,1.0)
) {

  # Data type
  if ( missing (response.type) ) {
    stop ("NP.control.lr.wls: missing 'response.type'")
  } else {
    response.type <- match.fun.arg (response.type)
  }

  # Check validity of criterion
  if ( missing (criterion) ) {
    criterion <- switch (
      response.type,
      "binomial" = { "LOGLIKELIHOOD" },
      "continuous" = { "PRESS" },
      stop (
        paste (
          "NP.control.lr.wls: can not choose defatult criterion for response.type '",
          response.type,
          "'.",
          sep=""
        )
      )
    )
  } else {
    criterion <- match.fun.arg (criterion)
  }

  # Valid optimization mode
  if ( missing(optim) ) {
    stop ("NP.control.lr.wls: missing 'optim'")
  } else {
    optim <- match.fun.arg (optim)
  }

  # Valid weight.function?
  if ( missing(weight.function) ) {
    stop ("NP.control.lr.wls: missing 'weight.function'")
  } else {
    weight.function <- match.fun.arg (weight.function)
  }
  if ( !exists(weight.function, mode="function") ) {
    stop ("NP.control.lr.wls: '",weight.function,"' not a weight function",sep="")
  }

  # Valid bandwidth?
  if ( missing(bandwidth.type) ) {
    stop ("NP.control.lr.wls: missing 'bandwidth.type'")
  } else {
    bandwidth.type <- match.fun.arg (bandwidth.type)
  }

  structure (
    list (
      response.type=response.type,
      criterion=criterion,
      optim=optim,
      weight.function=get(weight.function),
      bandwidth.type=bandwidth.type,
      bandwidth=bandwidth
    ),
    class="NP.ctr"
  )

}
