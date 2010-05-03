# Controllers for model roboust regression
# response.type (character, ?): the response type
# criterion (character, ?): criterion to optimize over wrt maxing and
#   optionally also bandwidth
# optim (character, ?): optimization mode
# mixing (numeric, atomic/vector): the actual mixing, use a range or
#   grid if appropriate with the chocen optimization mode
SP.control.mrr <- function (
  response.type=c("binomial", "continuous"),
  criterion=c("LOGLIKELIHOOD", "PRESS"),
  optim=c("none", "mixing.grid", "two.step", "joint"),
  mixing=c(0, 1)
) {

  # Check for validity of response.type parameter
  if ( missing(response.type) )  {
    stop ("SP.control.mrr: missing response.type")
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
          "SP.control.mrr: can not choose defatult criterion for response.type '",
          response.type,
          "'.",
          sep=""
        )
      )
    )
  } else {
    criterion <- match.fun.arg (criterion)
  }

  # Check for validity of optim parameter  
  if ( missing(optim) )  {
    # TODO: Is there concensus on this choice?
    # Default to joint
    optim <- "two.step"
  } else {
    optim <- match.fun.arg (optim)
  }

  # Check validity of mixing
  if (
    !(
      ( is.numeric(mixing) && is.vector(mixing) )
        ||
      ( mixing == FALSE )
    )
  ) {
    stop ("SP.control.mrr: mixing must be numeric vector of FALSE")
  }
  if ( is.numeric(mixing) ) {
    if ( !(0<=min(mixing) && max(mixing)<=1) ) {
      stop ("SP.control.mrr: mixing be subset of [0,1]")
    }
  }

  structure (
    list (response.type=response.type, criterion=criterion, optim=optim,
          mixing=mixing),
    class="SP.ctr"
  )

}
