# Hat matrix for a linear regression model at design points or subset thereof
# Observe that zero-weight obs drop out from design matrix X due to the
# internals of lm. Therefore NAs are put into design points of weight zero.
# This is perhaps not the most intuitive function.
# lm (class lm, object): lm fit object
# rows (integer, vector): rows (design points) at which hat matirx is evaluated
hat.design <- function (lm, rows=1:length(lm$residuals)) {

  QR <- lm$qr

  weights <- weights(lm)
  if ( !is.null(weights) ) {
    zero.weights <- (weights==0)
    K <- diag(!zero.weights)[!zero.weights,]
    W.sqrt <- sqrt(weights[!zero.weights])
    # Here trickery to fix rows for zero weights
    map <- 1:length(weights) - cumsum (zero.weights)
    map[zero.weights] <- NA
    rows <- map[rows]
    rows <- rows[!is.na(rows)]
  }

  if ( QR["rank"]==1 ) {

    n <- length (lm$residuals)
    if ( is.null(weights) ) {
      matrix (data=1/n, nrow=length(rows), ncol=n)
    } else {
      sum.w <- sum(weights)
      matrix (data=weights/sum.w, nrow=length(rows), ncol=n, byrow = TRUE)
    }

  } else if ( QR["rank"]==2 ) {

    Q <- qr.Q(QR)
    if ( is.null(weights) ) {
      Q[rows,] %*% t(Q)
    } else {
      t(
        apply (
          cbind( Q[rows,] , 1.0/(W.sqrt[rows]) )   ,
          1   ,
          function (x)  {x[1:2] * x[3] }
        )
      )  %*% t(Q) %*% diag(W.sqrt) %*% K
    }

  } else {

    stop ("hat.design: rank of QR decomposition not 1 or 2")

  }

}

# Hat matrix for a linear regression model at points with given covariate
# lm (class lm, object): lm fit object
# xs (numeric, vector): covariate values at which hat matrix is evaluated
hat.point <- function (lm, xs) {

  QR <- lm$qr

  weights <- weights(lm)
  if ( !is.null(weights) ) {
    zero.weights <- (weights==0)
    K <- diag(!zero.weights)[!zero.weights,]
    W.sqrt <- diag(sqrt(weights[!zero.weights]))
  }

  if ( QR["rank"]==1 ) {

    n <- length (lm$residuals)
    if ( is.null(weights) ) {
      matrix (data=1/n, nrow=length(xs), ncol=n)
    } else {
      sum.w <- sum(weights)
      matrix (data=weights/sum.w, nrow=length(xs), ncol=n, byrow = TRUE)
    }

  } else if ( QR["rank"]==2 ) {
  
    Q.t <- t(qr.Q(QR))
    R.inv <- solve(qr.R(QR))
    if ( is.null(weights) ) {
      cbind(1,xs) %*% R.inv %*% Q.t
    } else {
      cbind(1,xs) %*% R.inv %*% Q.t %*% W.sqrt %*% K
    }

  } else {

    stop ("hat.point: rank of QR decomposition not 1 or 2")

  }

}
