# Gaussian weight function
# dist (numeric, vector): distance
# bandwidth (numeric, atomic): bandwidth
gaussian.weight <- function(
  dist,
  bandwidth=1.0
) {

  if ( bandwidth<=0.0 ) {
    stop ("gaussian.weight: illegal bandwidth, positive value required (bandwidth>0)")
  }
  exp(-(dist/bandwidth)^2)

}
