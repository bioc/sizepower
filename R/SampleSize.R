##Sample Size Calculation for Completely Randomized Treatment-Control Designs 
sampleSize.randomized <- function(ER0, G0, power, absMu1, sigmad)
{
    if(ER0 < 0) { stop("'ER0' must be > 0!") }
    if(G0 < 0) { stop("'G0' must be > 0!") }
    if(power > 1 | power < 0)
    { stop("'power' must be in the interval [0, 1]!") }
    if(absMu1 < 0) { stop("'absMu1' must be > 0!") }
    if(sigmad < 0) { stop("'sigmad' must be > 0!") }
 
    ## type I error
    alpha0 <- ER0 / G0
    if(alpha0 > 1 | abs(alpha0 - 1) < 1.0e-10)
    { stop("'ER0' must be < 'G0'!") }
 
    a <- 1 - alpha0 / 2
    denom <- absMu1 / sigmad ## statistical difference
    za <- stats::qnorm(a)
    zb <- stats::qnorm(power)
    numer <- za + zb
 
    n <- (numer / denom)^2
    n <- max(1, ceiling(n))
 
    return(list(n=n, d=denom))
}

##sampleSize.randomized(ER0=1, G0=2000, power=0.9, absMu1=1, sigmad=0.566)

##Sample Size Calculation for Matched-Pairs Designs 
sampleSize.matched <- function(ER0, G0, power, absMu1, sigmad)
{
  if(ER0 < 0) { stop("'ER0' must be > 0!") }
  if(G0 < 0) { stop("'G0' must be > 0!") }
  if(power > 1 | power < 0)
  { stop("'power' must be in the interval [0, 1]!") }
  if(absMu1 < 0) { stop("'absMu1' must be > 0!") }
  if(sigmad < 0) { stop("'sigmad' must be > 0!") }

  ## type I error
  alpha0 <- ER0 / G0
  if(alpha0 > 1 | abs(alpha0 - 1) < 1.0e-10)
  { stop("'ER0' must be < 'G0'!") }

  a <- 1 - alpha0 / 2
  denom <- absMu1 / sigmad ## statistical difference
  za <- stats::qnorm(a)
  zb <- stats::qnorm(power)
  numer <- za + zb

  n <- (numer / denom)^2
  n <- max(1, ceiling(n))

  return(list(n=n, d=denom))
}

##sampleSize.matched(ER0=1, G0=2000, power=0.9, absMu1=1, sigmad=0.5)


