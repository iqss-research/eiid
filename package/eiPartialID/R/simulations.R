#' Simulate data; Example 3 from Jiang et al. 2019.
#'
#' \code{generateDataExample3} generates simulated data for example simulation 3 in Jiang et al. 2019.
#'
#' @return List object with the following attributes:
#' x (proportion of X in each geographic unit);
#' t (proportion of T in each geographic unit);
#' n (population of each geographic unit);
#' bd (the true district B);
#'
#' @importFrom stats runif
#'
#' @examples
#' generatedData <- generateDataExample3()
#'
#' @export
#'
generateDataExample3 <- function() {
  set.seed(7)
  p=1000 # number of precincts, in practice should be observed
  nn=150
  n=rep(nn,p) #precinct sizes, in practice should be observed
  x=runif(p,0,1) #proportions of blacks, in practice should be observed
  ww0=1 #white vote prob in pure white precinct, unavailable in practice
  ww1=0 #white vote prob in pure black precinct, unavailable in practice
  bb0=0  #black vote prob in pure white precinct, unavailable in practice
  bb1=0  #black vote prob in pure black precinct, unavailable in practice
  w0=log(ww0/(1-ww0))
  #in practice w0 should be estimated by regressing t<-w0+c1x+d1x^2
  w1=log(ww1/(1-ww1))-log(ww0/(1-ww0))
  b0=log(bb0/(1-bb0))
  b1=log(bb1/(1-bb1))-log(bb0/(1-bb0))
  b=rep(0.5,p)
  w=b
  s=0  #extra noise sd
  for(i in 1:p) {
    b[i]=0
    w[i]= 1-x[i]
  }
  t=(w*round(n*(1-x))+b*round(n*x))/n #comb'd vote prop, in practice is observed
  bd=sum(round(n*x)*b)/sum(round(n*x)) #true district b

  generatedData <- list()
  generatedData[["x"]] <- x
  generatedData[["t"]] <- t
  generatedData[["n"]] <- n
  generatedData[["bd"]] <- bd
  return(generatedData)
}

#' Simulate data; Example 4 from Jiang et al. 2019.
#'
#' \code{generateDataExample4} generates simulated data for example simulation 4 in Jiang et al. 2019.
#'
#' @return List object with the following attributes:
#' x (proportion of X in each geographic unit);
#' t (proportion of T in each geographic unit);
#' n (population of each geographic unit);
#' bd (the true district B);
#'
#' @importFrom stats runif rbinom rnorm
#'
#' @examples
#' generatedData <- generateDataExample4()
#'
#' @export
#'
generateDataExample4 <- function() {
  set.seed(7)
  p=1000 # number of precincts, in practice should be observed
  nn=150
  n=rep(nn,p) #precinct sizes, in practice should be observed
  x=runif(p,0,0.95) #proportions of blacks, in practice should be observed
  l=min(x)
  u=max(x)
  ww0=0.9 #white vote prob in pure white precinct, unavailable in practice
  ww1=0.9 #white vote prob in pure black precinct, unavailable in practice
  bb0=0.9 #black vote prob in pure white precinct, unavailable in practice
  bb1=0.6 #black vote prob in pure black precinct, unavailable in practice
  w0=log(ww0/(1-ww0))
  #in practice w0 should be estimated by regressing t<-w0+c1x+d1x^2
  w1=log(ww1/(1-ww1))-log(ww0/(1-ww0))
  b0=log(bb0/(1-bb0))
  b1=log(bb1/(1-bb1))-log(bb0/(1-bb0))
  b=rep(0.5,p)
  w=b
  s=0.5 #extra noise sd
  for(i in 1:p) {
    b[i]=rbinom(1,round(n[i]*x[i])+1,1/(1+exp(-b0-b1*x[i]+s*(1-x[i])*rnorm(1))))/(round(n[i]*x[i])+1)
    w[i]=rbinom(1,round(n[i]*(1-x[i]))+1,1/(1+exp(-w0-w1*x[i]+s*(1-x[i])*rnorm(1))))/(round(n[i]*(1-x[i]))+1)
  }
  t=(w*round(n*(1-x))+b*round(n*x))/n #comb'd vote prop, in practice is observed
  bd=sum(round(n*x)*b)/sum(round(n*x)) #true district b
  generatedData <- list()
  generatedData[["x"]] <- x
  generatedData[["t"]] <- t
  generatedData[["n"]] <- n
  generatedData[["bd"]] <- bd
  return(generatedData)
}

#' Simulate data; Example 5 from Jiang et al. 2019.
#'
#' \code{generateDataExample5} generates simulated data for example simulation 5 in Jiang et al. 2019.
#'
#' @return List object with the following attributes:
#' x (proportion of X in each geographic unit);
#' t (proportion of T in each geographic unit);
#' n (population of each geographic unit);
#' bd (the true district B);
#'
#' @importFrom stats runif rnorm
#'
#' @examples
#' generatedData <- generateDataExample5()
#'
#' @export
#'
generateDataExample5 <- function() {
  set.seed(7)
  p=1000 # number of precincts, in practice should be observed
  nn=150
  n=rep(nn,p) #precinct sizes, in practice should be observed
  x=runif(p,0,0.7) #proportions of blacks, in practice should be observed
  l=min(x)
  u=max(x)
  ww0=0.9 #white vote prob in pure white precinct, unavailable in practice
  ww1=0.9 #white vote prob in pure black precinct, unavailable in practice
  bb0=0.5 #black vote prob in pure white precinct, unavailable in practice
  bb1=0.5 #black vote prob in pure black precinct, unavailable in practice
  w0=log(ww0/(1-ww0))
  #in practice w0 should be estimated by regressing t<-w0+c1x+d1x^2
  w1=log(ww1/(1-ww1))-log(ww0/(1-ww0))
  b0=log(bb0/(1-bb0))
  b1=log(bb1/(1-bb1))-log(bb0/(1-bb0))
  b=rep(0.5,p)
  w=b
  s=1 #extra noise sd
  for(i in 1:p) {
    b[i]=rbinom(1,round(n[i]*x[i])+1,1/(1+exp(-b0-b1*x[i]+s*rnorm(1))))/(round(n[i]*x[i])+1)
    w[i]=rbinom(1,round(n[i]*(1-x[i]))+1,1/(1+exp(-w0-w1*x[i]+s*rnorm(1))))/(round(n[i]*(1-x[i]))+1)
  }
  t=(w*round(n*(1-x))+b*round(n*x))/n #comb'd vote prop, in practice is observed
  bd=sum(round(n*x)*b)/sum(round(n*x)) #true district b
  generatedData <- list()
  generatedData[["x"]] <- x
  generatedData[["t"]] <- t
  generatedData[["n"]] <- n
  generatedData[["bd"]] <- bd
  return(generatedData)
}


#' Run illustrative simulations.
#'
#' \code{runExampleSimulations} generates bounds for the example simulations in Jiang et al. 2019.
#'
#' @return No explicit return values. The summary of the simulations will be printed to standard out.
#'
#' @examples
#' runExampleSimulations()
#'
#' @export
#'
runExampleSimulations <- function() {
  cat("Example 3\n")
  generatedData <- generateDataExample3()
  # note that generatedData$bd is the district level DD estimate, but it is passed through as if it's the true beta_i (precinct level), which is fine here since the aggregation leaves the value unchanged
  ex3 <- generateBounds(generatedData$x, generatedData$t, generatedData$n, generatedData$bd, useXRangeOffset=FALSE, returnAdditionalStats=FALSE, printSummary=TRUE)
  cat("\nExample 4\n")
  generatedData <- generateDataExample4()
  ex4 <- generateBounds(generatedData$x, generatedData$t, generatedData$n, generatedData$bd, useXRangeOffset=FALSE, returnAdditionalStats=FALSE, printSummary=TRUE)
  cat("\nExample 5\n")
  generatedData <- generateDataExample5()
  ex5 <- generateBounds(generatedData$x, generatedData$t, generatedData$n, generatedData$bd, useXRangeOffset=FALSE, returnAdditionalStats=FALSE, printSummary=TRUE)
}

