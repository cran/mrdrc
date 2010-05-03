# Controller for binomial data parametric engine for model robust regression
# link (cheracter, ?): link to be used in glm
P.control.glm.binomial.ml <- function (
  link=c("logit","probit","cauchit")
) {

  # Check validity of link
  if ( missing(link) ) {
    link <- "logit"
  } else {
    link <- match.fun.arg (link)
  }

  structure (
    list (response.type="binomial", family="binomial", link=link),
    class="P.ctr"
  )

}

# Controller for continuous data parametric engine for model robust regression
# fun (cheracter, ?): mean structure
# dfun (cheracter, ?): derivative of fun with respect to parameter
# ssfun (cheracter, ?): selfstarter
P.control.nl.ls <- function (
  fun="logistic4",
  dfun=paste("d.",fun,sep=""),
  ssfun=paste("ss.",fun,sep="")
) {

  # Check existence of mean, derivative and selfstarter
  sapply (
    list(fun,dfun,ssfun),
    function (str) {
      if ( !exists(str,mode="function") ) {
        stop (paste("P.control.nl.ls: '",str,"' not a function",sep=""))
      }
    }
  )

  structure (
    list (response.type="continuous", fun=fun, dfun=dfun, ssfun=ssfun),
    class="P.ctr"
  )

}
