hatMatrix <- function (x, ...) UseMethod ("hatMatrix", x)
leaveOneOutPredict <- function (x, ...) UseMethod ("leaveOneOutPredict", x)
CI.H <- function (x, ...) UseMethod ("CI.H", x)
ED.H <- function (x, ...) UseMethod ("ED.H", x)

# Non-parametric
predict.NP <- function (NP.obj, new.data) {

##  message ("\n*** NP predict dispatcher")
  NP.obj$predict(new.data=new.data)

}

leaveOneOutPredict.NP <- function (NP.obj) {

##  message ("\n*** NP leave.one.out.predict dispatcher")
  NP.obj$leave.one.out.predict()

}

hatMatrix.NP <- function (NP.obj, new.data) {

##  message ("\n*** NP hatMatrix dispatcher")
  NP.obj$hatMatrix(new.data=new.data)

}

CI.H.NP <- function (NP.obj, new.data) {

##  message ("\n*** NP CI.H dispatcher")
  NP.obj$CI.H(new.data=new.data)

}

print.NP <- function (NP.obj, verbose=FALSE) {

##  message ("\n*** NP print dispatcher")
  NP.obj$print(NP.obj=NP.obj, verbose=verbose)

}

plot.NP <- function (NP.obj, ...) {

##  message ("\n*** NP plot dispatcher")
  NP.obj$plot(NP.obj=NP.obj, ...)

}

# Parametric
predict.P <- function (P.obj, new.data) {

##  message ("\n*** P predict dispatcher")
  P.obj$predict(new.data=new.data)

}

leaveOneOutPredict.P <- function (P.obj) {

##  message ("\n*** P leave.one.out.predict dispatcher")
  P.obj$leave.one.out.predict()

}

hatMatrix.P <- function (P.obj, new.data) {

##  message ("\n*** P hatMatrix dispatcher")
  P.obj$hatMatrix(new.data=new.data)

}

CI.H.P <- function (P.obj, new.data) {

##  message ("\n*** P CI.H dispatcher")
  P.obj$CI.H(new.data=new.data)

}

print.P <- function (P.obj, verbose=FALSE) {

##  message ("\n*** P print dispatcher")
  P.obj$print(P.obj=P.obj, verbose=verbose)

}

plot.P <- function (P.obj, ...) {

##  message ("\n*** P plot dispatcher")
  P.obj$plot(P.obj=P.obj, ...)

}

# Semi-parametric
predict.SP <- function (SP.obj, new.data) {

##  message ("\n*** SP predict dispatcher")
  SP.obj$predict(new.data=new.data)

}

leaveOneOutPredict.SP <- function (SP.obj) {

##  message ("\n*** SP leave.one.out.predict dispatcher")
  SP.obj$leave.one.out.predict()

}

hatMatrix.SP <- function (SP.obj, new.data) {

##  message ("\n*** SP hatMatrix dispatcher")
  SP.obj$hatMatrix(new.data=new.data)

}

CI.H.SP <- function (SP.obj, new.data) {

##  message ("\n*** SP CI.H dispatcher")
  SP.obj$CI.H(new.data=new.data)

}

## Modified 15/1 by Christian Ritz
ED.H.SP <- function (SP.obj, alpha, reference, level) {

##  message ("\n*** SP ED.H dispatcher")
  SP.obj$ED.H(alpha=alpha, reference = reference, level = level)

}

print.SP <- function (SP.obj, verbose=FALSE) {

##  message ("\n*** SP print dispatcher")
  SP.obj$print(SP.obj=SP.obj, verbose=verbose)

}

plot.SP <- function (SP.obj, ...) {

##  message ("\n*** SP plot dispatcher")
  SP.obj$plot(SP.obj=SP.obj, ...)

}
