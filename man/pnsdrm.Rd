\name{pnsdrm}

\alias{pnsdrm}
\alias{pnsdrm.calc}

\alias{pns.plot1}

\title{Parametric, non-parametric or semi-parametric dose-response modelling}

\description{
  Parametric, non-parametric or semi-parametric dose-response modelling of both continuous and quantal data.
}

\usage{
pnsdrm(predictor, response, weights, type = c("continuous", "binomial"), 
model = c("semi-parametric", "non-parametric", "parametric"), 
fct = NULL, robust = FALSE, respLev = c(10, 20, 50), 
reference = NULL, level = 0.95, logex = FALSE)

pnsdrm.calc(predictor, response, weights, type = c("continuous", "binomial"), 
model = c("semi-parametric", "non-parametric", "parametric"), 
fct = NULL, robust = FALSE, respLev = c(10, 20, 50), 
reference = NULL, level = 0.95, logex = FALSE)
}

\arguments{
  \item{predictor}{numeric vector of concentrations/doses.}
  \item{response}{numeric vector of response values (proportions in case of quantal data).}
  \item{weights}{numeric vector of weights needed for quantal data.}
  \item{type}{character string specifying the type of response.}
  \item{model}{character string specifying the model to be fit.}
  \item{fct}{a built-in function or a list of built-in functions from the package 'drc'.}
  \item{robust}{logical specifying whether or not a robust approach should be used. Only for the
  semi-parametric approach.}
  \item{respLev}{numeric vector of requested ED level.}
  \item{reference}{optional reference value for the lower limit.}
  \item{level}{numeric specifying the confidence level.}
  \item{logex}{logical indicating whether or not a logarithmic x axis should be used.}
}

\details{
  The parametric estimation is based on the model fitting function \code{\link[drc]{drm}} in the package 'drc'. 
  The non-parametric estimation relies on the 'locfit' package.
  
  The semi-parametric approach is mainly based on the development in Nottingham and Birch (2000), whereas the
  non-parametric approach uses on the package 'EffectiveDose' which implements the method introduced in 
  Dette \emph{et al} (2004).
  
  \code{plot} and \code{print} methods are available. 
}

\value{
  A list containing the requested ED values and additional information about the underlying 
  model fit(s).
}

\references{
  Dette, H., Neumeyer, N. and Pilz, K. F. (2004) A Note on Nonparametric Estimation of the Effective Dose
  in Quantal Bioassay, \emph{J. Amer. Statist. Assoc.}, \bold{100}, 503--510.

  Nottingham, Q. and Birch, J. B. (2000) A Semiparametric Approach to Analysing Dose-Response Data,
  \emph{Statist. Med.}, \bold{19}, 389--404.
}

\author{
  Christian Ritz (wrapper functions)
  Mads Jeppe Tarp-Johansen (internal functions)
}

%\note{
%  The implementation of this function as well as all other functions in the package 'mrdrc' has been funded by
%  European Centre for the Validation of Alternative Methods, EU Joint Research Centre under lot 3 of the 
%  project "Quality assessment and novel statistical analysis techniques for toxicological data".
%}

%\seealso{
%  More examples are found in the help pages for \code{\link{bin.mat}} and \code{\link{exp.a}}.
%}

\examples{

## Analysing deguelin (in the package 'drc')

## Semi-parametric model
deguelin.mrr1 <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = LL.2())
deguelin.mrr1
plot(deguelin.mrr1)

## The same
gmFct <- getMeanFunctions(fname = "LL.2")
deguelin.mrr1b <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = gmFct)
deguelin.mrr1b
plot(deguelin.mrr1b)

## The same again
deguelin.mrr1c <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = list(LL2.2()))
deguelin.mrr1c
plot(deguelin.mrr1c)


deguelin.mrr1d <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = W1.2())
deguelin.mrr1d
plot(deguelin.mrr1d)

## The same
gmFct <- getMeanFunctions(fname = "W1.2")
deguelin.mrr1e <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = gmFct)
deguelin.mrr1e
plot(deguelin.mrr1e)


### Parametric models
#deguelin.mrr2 <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
#model = "parametric", fct = list(LL.2(), W1.2(), W2.2()))
#deguelin.mrr2
#plot(deguelin.mrr2)

### The same parametric models
#deguelin.mrr2b <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
#model = "parametric", fct = list(W2.2(), LL.2(), W1.2()))
#deguelin.mrr2b
#plot(deguelin.mrr2b)

## Non-parametric approach -- currently not available
#deguelin.mrr3 <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
#model = "non-parametric")
#deguelin.mrr3
#plot(deguelin.mrr3)


## Semi-parametric model with reference level 0.3
deguelin.mrr4 <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = LL.2(), reference = 0.3)
deguelin.mrr4
plot(deguelin.mrr4)


## Semi-parametric models
deguelin.mrr5 <- pnsdrm(deguelin$dose, deguelin$r, deguelin$n, type = "binomial",
model = "semi-parametric", fct = list(LL.2(), W1.2(), W2.2()))
deguelin.mrr5
plot(deguelin.mrr5)


## Analysing ryegrass (in the package 'drc')

ryegrass.mrr1 <- pnsdrm(ryegrass$conc, ryegrass$rootl, type = "continuous",
model = "semi-parametric", fct = LL.5())
ryegrass.mrr1
plot(ryegrass.mrr1)
plot(ryegrass.mrr1, log = "x")

ryegrass.mrr2 <- pnsdrm(ryegrass$conc, ryegrass$rootl, type = "continuous",
model = "semi-parametric", fct = list(LL.3(), LL.4(), LL.5()))
ryegrass.mrr2
plot(ryegrass.mrr2)

#ryegrass.mrr3 <- pnsdrm(ryegrass$conc, ryegrass$rootl, type = "continuous",
#model = "parametric", fct = list(LL.3(), LL.4(), LL.5()))
#ryegrass.mrr3
#plot(ryegrass.mrr3)

ryegrass.mrr4 <- pnsdrm(ryegrass$conc, ryegrass$rootl, type = "continuous",
model = "semi-parametric", fct = list(L.4(), LL.4(), W1.4(), W2.4()))
ryegrass.mrr4
plot(ryegrass.mrr4)


## Analysing lettuce (in the package 'drc')

lettuce.mrr1 <- pnsdrm(lettuce$conc, lettuce$weight, type = "continuous",
model = "semi-parametric", fct = LL.3())
lettuce.mrr1
plot(lettuce.mrr1)

lettuce.mrr2 <- pnsdrm(lettuce$conc, lettuce$weight, type = "continuous",
model = "semi-parametric", fct = BC.4())
lettuce.mrr2
plot(lettuce.mrr2)

#lettuce.mrr3 <- pnsdrm(lettuce$conc, lettuce$weight, type = "continuous",
#model = "semi-parametric", fct = LL.3(), robust = TRUE)

#lettuce.mrr3
#plot(lettuce.mrr3)

}
\keyword{models}
\keyword{nonlinear}
