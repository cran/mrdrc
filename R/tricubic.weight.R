# Tricubic weight function
# dist (numeric, vector): distance
# bandwidth (numeric, atomic): bandwidth
tricubic.weight <- function(
  dist,
  bandwidth=1.0
) {
  if ( bandwidth<=0.0 ) {
    stop ("tricubic.weight: illegal bandwidth, positive value required (bandwidth>0)")
  }
  (1-abs(dist/bandwidth)^3)^3 * (-bandwidth<dist & dist<bandwidth)
}
