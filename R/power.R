## power calculation for Multiple-Treatment Designs 
## Having an Isolated Treatment Effect 
power.multi <- function(ER0, G0, numTrt, absMu1, sigma, n)
{
    if(ER0 < 0) { stop("'ER0' must be > 0!") }
    if(G0 < 0) { stop("'G0' must be > 0!") }
    if(abs(numTrt-as.integer(numTrt)) > 0 | numTrt < 2)
    { stop("'numTrt' must be integer and > 1!") }
    if(absMu1 < 0) { stop("'absMu1' must be > 0!") }
    if(sigma < 0) { stop("'sigma' must be > 0!") }
    if((abs(n - as.integer(n)) > 0) | (n < 0)) 
    { stop("'n' must be positive integer!") }
 
    ## type I error
    alpha0 <- ER0 / G0
    if(alpha0 > 1 | abs(alpha0 - 1) < 1.0e-10)
    { stop("'ER0' must be < 'G0'!") }
 
    prob <- 1 - alpha0
    df <- numTrt - 1 ## degree of freedom
    ## decision rule, if test statistic > myc, then reject H0
    myc <- stats::qchisq(prob, df)
 
    ## non-centrality parameter
    psi1 <- n * (numTrt - 1.0) / numTrt * (absMu1 / sigma)^2
 
    ## power=Pr(reject H0 | Ha true)
    res <- 1 - stats::pchisq(myc, df=df, ncp=psi1)
    return(list(power=res, psi1=psi1))
}

##power.multi(ER0=2,G0=10000, numTrt=6, absMu1=0.585, sigma=0.3,n=8)

## power calculation for Matched-Pairs Designs 
power.matched <- function(ER0, G0, absMu1, sigmad, n)
{
    if(ER0 < 0 ) { stop("'ER0' must be > 0!") }
    if(G0 < 0) { stop("'G0' must be > 0!") }
    if(absMu1 < 0) { stop("'absMu1' must be > 0!") }
    if(sigmad < 0) { stop("'sigma' must be > 0!") }
    if((abs(n - as.integer(n)) > 0) | (n < 0)) 
    { stop("'n' must be positive integer!") }
 
    ## type I error
    alpha0 <- ER0 / G0
    if(alpha0 > 1 | abs(alpha0 - 1) < 1.0e-10)
    { stop("'ER0' must be < 'G0'!") }
 
    prob <- 1 - alpha0
    numTrt <- 2 ## number of treatments
    df <- numTrt - 1 ## degree of freedom
    ## decision rule, if test statistic > myc, then reject H0
    myc <- stats::qchisq(prob, df)
 
    ## non-centrality parameter
    psi1 <- n * (absMu1 / sigmad)^2
 
    ## power=Pr(reject H0 | Ha true)
    res <- 1 - stats::pchisq(myc, df=df, ncp=psi1)
    return(list(power=res, psi1=psi1))
}

##power.matched(ER0=2,G0=5000, absMu1=1, sigmad=0.4243,n=4)

## power calculation for Completely Randomized Treatment-Control Designs 
power.randomized <- function(ER0, G0, absMu1, sigmad, n)
{
    if(ER0 < 0) { stop("'ER0' must be > 0!") }
    if(G0 < 0) { stop("'G0' must be > 0!") }
    if(absMu1 < 0) { stop("'absMu1' must be > 0!") }
    if(sigmad < 0) { stop("'sigma' must be > 0!") }
    if((abs(n - as.integer(n)) > 0) | (n < 0)) 
    { stop("'n' must be positive integer!") }
 
    ## type I error
    alpha0 <- ER0 / G0
    if(alpha0 > 1 | abs(alpha0 - 1) < 1.0e-10)
    { stop("'ER0' must be < 'G0'!") }
 
    prob <- 1 - alpha0
    numTrt <- 2 ## number of treatments
    df <- numTrt - 1 ## degree of freedom
    ## decision rule, if test statistic > myc, then reject H0
    myc <- stats::qchisq(prob, df)
 
    ## non-centrality parameter
    psi1 <- n * (absMu1 / sigmad)^2
 
    ## power=Pr(reject H0 | Ha true)
    res <- 1 - stats::pchisq(myc, df=df, ncp=psi1)
    return(list(power=res, psi1=psi1))
}

##power.randomized(ER0=2,G0=5000, absMu1=1, sigmad=0.5657,n=8)

