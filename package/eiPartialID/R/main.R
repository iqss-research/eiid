#' Internal/private method. Compute bounds and summary statistics according to Jiang et al. 2019
#'
#' \code{calcSummaryOutputValues_()} is an internal/private helper method for calculating the bounds. Called by calcSummaryOutputValues().
#'
#' @param x Numeric (double-precision) vector. Contains the proportion of variable X in each precinct (or analagous geographic unit)
#' @param t Numeric (double-precision) vector. Contains the proportion of variable T in each precinct (or analagous geographic unit)
#' @param n Numeric (double-precision) vector. Contains the number of elements (people/households/etc.) in each precinct (or analagous geographic unit)
#' @param useXRangeOffset boolean If True, an offset of 0.00001 is applied to l and u to avoid division by 0 in subsequent calculations. Default TRUE.
#'
#' @return List object with the bounds and summary statistics
#'
#' @examples
#' \dontrun{
#' outputList <- calcSummaryOutputValues_(x, t, n, TRUE)
#' }
#'
calcSummaryOutputValues_ <- function(x, t, n, useXRangeOffset=TRUE) {

  if (!(useXRangeOffset)) {
    l=min(x)
    u=max(x)
  } else {
    # the offset here is to handle division by zero
    l=min(x)+0.00001  # lower end of x range where we assume linear contextual model
    u=max(x)-0.00001  # upper end of x range where we assume linear contextual model
  }

  p=length(x) # number of precincts

  nx1=sum(n*x)
  ###########################
  # DD bounds:
  bl <- numeric(p)
  bu <- numeric(p)
  for (i in 1:p) {
    lowerBound <- max( 0, (t[i] + x[i] - 1) / x[i] )
    upperBound <- min(1, t[i] / x[i] )
    if (!is.finite(lowerBound)) {
      lowerBound <- 0.0
    }
    if (!is.finite(upperBound)) {
      upperBound <- 1.0
    }
    bl[i] <- lowerBound
    bu[i] <- upperBound
  }
  ###########################
  # mean of midpoint of deterministic bound (for reference); note that this is not scaled/weighted by the population
  ddMidPointMean <- mean((bu-bl)/2 + bl)
  ###########################

  r=sum(n*x*(1-x))/sum(n*x)
  lb=0 # this is our choice of lambda (always 0 in the current version)
  h0=lb*sum(n*x*t)/sum(n*x)
  s1= sqrt(sum((n*x*((1+lb)/2-lb*x))^2)/(sum(n*x))^2)
  cl=0.90
  #this is conf level
  ex=1
  gl01=0
  gl1=c(-1/l,0,0)
  gl02=-1/(1-l)
  gl2=c(1/(1-l),1/(1-l),1/(1-l)-1)
  gu01=1/l
  gu1=c(-1/l,0,0)
  gu02=0
  gu2=c(1/(1-l),1/(1-l),1/(1-l)-1)
  ###
  gl03=0
  gl3=c(-1/u,0,0)
  gl04=-1/(1-u)
  gl4=c(1/(1-u),1/(1-u),1/(1-u)-1)
  gu03=1/u
  gu3=c(-1/u,0,0)
  gu04=0
  gu4=c(1/(1-u),1/(1-u),1/(1-u)-1)
  h=c(1-lb, sum(n*x*(1-lb*x))/sum(n*x),sum(n*x^2*(1-lb*x))/sum(n*x))
  ######

  #in practice w0 should be estimated by regressing t<-w0+c1x+d1x^2
  #d1=b1-w1 #in practice d1 should be estimated by regressing t<-w0+c1x+d1x^2
  #c1=b0-w0+w1  #in practice c1 should be estimated by regressing t<-w0+c1x+d1x^2
  reg=lm(formula = t ~ poly(x, 2,raw=TRUE)) # do quadratic regression here
  w0=coef(reg)[1]
  c1=coef(reg)[2]
  d1=coef(reg)[3]

  th=c(w0,c1,d1)
  v=sandwich::vcovHC(reg,"HC1")
  #this gives sandwich variance of STATA see
  #https://stats.stackexchange.com/questions/117052/replicating-statas-robust-option-in-r
  #information re robust sandwich variance formula:
  #https://www.stata.com/manuals/p_robust.pdf
  #
  sl1=s1+sqrt( t(h-r*gu1)%*%v%*% (h-r*gu1))
  sl2=s1+sqrt( t(h-r*gu2)%*%v%*% (h-r*gu2))
  sl3=s1+sqrt( t(h-r*gu3)%*%v%*% (h-r*gu3))
  sl4=s1+sqrt( t(h-r*gu4)%*%v%*% (h-r*gu4))
  sl=c(sl1,sl2,sl3,sl4)
  su1=s1+sqrt( t(h-r*gl1)%*%v%*% (h-r*gl1))
  su2=s1+sqrt( t(h-r*gl2)%*%v%*% (h-r*gl2))
  su3=s1+sqrt( t(h-r*gl3)%*%v%*% (h-r*gl3))
  su4=s1+sqrt( t(h-r*gl4)%*%v%*% (h-r*gl4))
  su=c(su1,su2,su3,su4)
  bl1=h0-r*gu01+t(h-r*gu1)%*%th
  bl2=h0-r*gu02+t(h-r*gu2)%*%th
  bl3=h0-r*gu03+t(h-r*gu3)%*%th
  bl4=h0-r*gu04+t(h-r*gu4)%*%th
  bbl=c(bl1,bl2,bl3,bl4)
  bu1=h0-r*gl01+t(h-r*gl1)%*%th
  bu2=h0-r*gl02+t(h-r*gl2)%*%th
  bu3=h0-r*gl03+t(h-r*gl3)%*%th
  bu4=h0-r*gl04+t(h-r*gl4)%*%th
  bbu=c(bu1,bu2,bu3,bu4)

  wuc=c(gu01+t(gu1)%*%th,gu02+t(gu2)%*%th,gu03+t(gu3)%*%th,gu04+t(gu4)%*%th)
  wu=min(gu01+t(gu1)%*%th,gu02+t(gu2)%*%th,gu03+t(gu3)%*%th,gu04+t(gu4)%*%th)
  #upperbound of w1
  wl=max(gl01+t(gl1)%*%th,gl02+t(gl2)%*%th,gl03+t(gl3)%*%th,gl04+t(gl4)%*%th)
  #lowerbound of w1

  cil=max(bbl)-ex*sl[which.max(bbl)]
  cir=min(bbu)+ex*su[which.min(bbu)]

  #these are the conservative ci at ex=1.
  ###########################

  hbdu0=min(bbu) #hat district level upperbound proposed (wtd)
  hbdl0=max(bbl) #hat district level lowererbound proposed (wtd)

  bdu=sum(n*x*bu)/sum(n*x) # district level upperbound by duncan-davis (wtd)
  bdl=sum(n*x*bl)/sum(n*x) # district level upperbound by duncan-davis (wtd)

  outputList <- list()
  outputList[["nx1"]] <- nx1 # total elements (people/households/etc.) of variable X across all geographic units
  outputList[["hbdl0"]] <- hbdl0 # CI_0 lower
  outputList[["hbdu0"]] <- hbdu0 # CI_0 upper
  outputList[["cil"]] <- cil # CI_1 lower
  outputList[["cir"]] <- cir # CI_1 upper
  outputList[["bdl"]] <- bdl # DD lower
  outputList[["bdu"]] <- bdu ## DD upper
  # for reference/analysis:
  outputList[["w0"]] <- w0
  outputList[["c1"]] <- c1
  outputList[["d1"]] <- d1
  outputList[["p"]] <- p
  outputList[["l"]] <- l
  outputList[["u"]] <- u
  outputList[["ddMidPointMean"]] <- ddMidPointMean

  return(outputList)


}

#' Internal/private method. Compute bounds and summary statistics according to Jiang et al. 2019
#'
#' \code{calcSummaryOutputValues()} is an internal/private helper method for calculating the bounds. Called by generateBounds().
#'
#' @param x Numeric (double-precision) vector. Contains the proportion of variable X in each precinct (or analagous geographic unit)
#' @param t Numeric (double-precision) vector. Contains the proportion of variable T in each precinct (or analagous geographic unit)
#' @param n Numeric (double-precision) vector. Contains the number of elements (people/households/etc.) in each precinct (or analagous geographic unit)
#' @param trueBetaB Numeric (double-precision) vector. Contains the true conditional values (beta_i) in each precinct (or analagous geographic unit).
#'  Optional. Default NULL.
#' @param useXRangeOffset boolean If True, an offset of 0.00001 is applied to l and u to avoid division by 0 in subsequent calculations. Default TRUE
#' @param returnAdditionalStats boolean If True, additional summary statistics are generated. Default FALSE.
#'
#' @return List object with the bounds and summary statistics
#'
#' @examples
#' \dontrun{
#' outputList <- calcSummaryOutputValues(x, t, n, NULL, TRUE, FALSE)
#' }
#'
calcSummaryOutputValues <- function(x, t, n, trueBetaB=NULL, useXRangeOffset=TRUE, returnAdditionalStats=FALSE) {

  outputList <- calcSummaryOutputValues_(x, t, n, useXRangeOffset)

  if (!is.null(trueBetaB)) {
    bd=sum(round(n*x)*trueBetaB)/sum(round(n*x)) # true district b
    outputList[["bd"]] <- bd
  }
  if (returnAdditionalStats) {
    # These values may be of interest for evaluation, comparisons, and future work:
    sdX <- sd(x)
    mean_nx_squared <- mean(n*x^2)
    mean_nx <- mean(n*x)
    mean_n <- mean(n)
    sd_nx <- sd(n*x)
    meanT <- mean(t)

    reg3=lm(formula = t ~ poly(x, 3,raw=TRUE))  # cubic regression
    w30=coef(reg3)[1]
    w31=coef(reg3)[2]
    w32=coef(reg3)[3]
    w33=coef(reg3)[4]
    v3=sandwich::vcovHC(reg3,"HC1")
    t3=w33/sqrt(v3[4,4])

    greg=lm(formula = t ~ poly(x, 1,raw=TRUE)) # For reference, Goodman Regression is linear regression with T,X
    ag=coef(greg)[1]
    bg=coef(greg)[2]
    gdmn=ag+bg # Goodman Regression result

    outputList[["sdX"]] <- sdX
    outputList[["mean_nx_squared"]] <- mean_nx_squared
    outputList[["mean_nx"]] <- mean_nx
    outputList[["mean_n"]] <- mean_n
    outputList[["sd_nx"]] <- sd_nx
    outputList[["meanT"]] <- meanT
    outputList[["w30"]] <- w30
    outputList[["w31"]] <- w31
    outputList[["w32"]] <- w32
    outputList[["w33"]] <- w33
    outputList[["t3"]] <- t3
    outputList[["gdmn"]] <- gdmn
  }

  return(outputList)
}

#' Compute bounds and summary statistics according to Jiang et al. 2019
#'
#' \code{generateBounds()} calculates district-level bounds. The returned object can be passed to evaluateBounds() to generate bounds across varying coverage probabilities
#' and to apply the heuristics presented in Jiang et al. 2019.
#'
#' @param x Numeric (double-precision) vector. Contains the proportion of variable X in each precinct (or analagous geographic unit)
#' @param t Numeric (double-precision) vector. Contains the proportion of variable T in each precinct (or analagous geographic unit)
#' @param n Numeric (double-precision) vector. Contains the number of elements (people/households/etc.) in each precinct (or analagous geographic unit)
#' @param trueBetaB Numeric (double-precision) vector. Contains the true conditional values (beta_i) in each precinct (or analagous geographic unit).
#'  Optional. Default NULL.
#' @param useXRangeOffset boolean If True, an offset of 0.00001 is applied to l and u to avoid division by 0 in subsequent calculations. Default TRUE
#' @param returnAdditionalStats boolean If True, additional summary statistics are generated. Default FALSE.
#' @param printSummary boolean If True, the DD bounds, l and u, CI_0, CI_1, width-ratio, and (optionally) true district B are output to standard out. Default TRUE.
#'
#' @return List object with the bounds and summary statistics:
#'
#'  nx1 Total elements (people/households/etc.) of variable X across all geographic units
#'
#'  hbdl0 CI_0 lower bound
#'
#'  hbdu0 CI_0 upper bound
#'
#'  cil CI_1 lower bound
#'
#'  cir CI_1 upper bound
#'
#'  bdl Duncan-Davis lower bound
#'
#'  bdu Duncan-Davis upper bound
#'
#'  Optional: bd True district Beta
#'
#' @examples
#' \dontrun{
#' library("MASS")
#' library("eco")
#' data("census")
#' inputDataSet <- census
#' x <- inputDataSet$X
#' t <- inputDataSet$Y
#' n <- inputDataSet$N
#' trueBetaB <- inputDataSet$W1
#' outputList <- generateBounds(x, t, n, trueBetaB=trueBetaB, useXRangeOffset=TRUE, returnAdditionalStats=FALSE, printSummary=TRUE)
#'
#' # True B: 0.674809
#' # Duncan-Davis bounds: [0.535618, 0.974010]
#' # [l,u]=[min(X_i),max(X_i)]: [0.050810, 0.939290]
#' # CI_0=[Bl_hat, Bu_hat]: [0.606101, 0.810082]
#' # CI_1: [0.572566, 0.842403]
#' # Width-ratio: |CI_0|/|DD|: 0.465295
#' }
#'
#' @export
#'
generateBounds <- function(x, t, n, trueBetaB=NULL, useXRangeOffset=TRUE, returnAdditionalStats=FALSE, printSummary=TRUE) {

  outputList <- calcSummaryOutputValues(x, t, n, trueBetaB, useXRangeOffset, returnAdditionalStats)

  if (printSummary) {
    if (!is.null(trueBetaB)) {
      cat(sprintf("True B: %f\n", outputList$bd))
    }

    cat(sprintf("Duncan-Davis bounds: [%f, %f]\n", outputList$bdl, outputList$bdu))
    cat(sprintf("[l,u]=[min(X_i),max(X_i)]: [%f, %f]\n", outputList$l, outputList$u))
    cat(sprintf("CI_0=[Bl_hat, Bu_hat]: [%f, %f]\n", outputList$hbdl0, outputList$hbdu0))
    cat(sprintf("CI_1: [%f, %f]\n", outputList$cil, outputList$cir))
    #this width ratio measures efficiency of proposed interval:
    cat(sprintf("Width-ratio: |CI_0|/|DD|: %f\n", (outputList$hbdu0-outputList$hbdl0)/(outputList$bdu-outputList$bdl)))

  }

  return(outputList)

}

#' Evaluate computed bounds, across confidence levels, applying the selection heuristic of Jiang et al. 2019
#'
#' \code{evaluateBounds()} calculates the bounds across confidence levels and generates the width-ratio relative to the deterministic DD bounds using the
#' bounds generated by generateBounds(), after applying the selection heuristic of Jiang et al. 2019. If the true district B is provided,
#' the capture of the true value is checked.
#'
#' @param outputListFromGenerateBounds List returned by generateBounds()
#'
#' @return List object with the bounds indexed across confidence levels:
#'
#'  x_for_x_in_CI_x c(0.00,0.25,0.50,0.75, 1.00, 1.25, 1.50, 1.75, 2.00), which corresponds to CI_0, CI_0.25, ..., CI_2.00 (the following vectors are parallel in indexes)
#'
#'  CI_x_lower CI_x lower bound
#'
#'  CI_x_upper CI_x upper bound
#'
#'  CI_x_isSelected If FALSE, proposed bound was not rejected by the heuristic (if TRUE, bounds are reverted to the DD bounds)
#'
#'  CI_x_widthRatio |CI_x|/|DD|
#'
#'  CI_x_nominalCoverage Nominal coverage (1-pnorm(-x_for_x_in_CI_x))
#'
#'
#'  Optional: CI_x_truthCaptured If true district Beta is provided in outputListFromGenerateBounds, then this vector contains a boolean for whether or not the
#'   true value was captured in the proposed CI_x.
#'
#' @examples
#' \dontrun{
#' library("MASS")
#' library("eco")
#' data("census")
#' inputDataSet <- census
#' x <- inputDataSet$X
#' t <- inputDataSet$Y
#' n <- inputDataSet$N
#' trueBetaB <- inputDataSet$W1
#' outputList <- generateBounds(x, t, n, trueBetaB=trueBetaB, useXRangeOffset=TRUE, returnAdditionalStats=FALSE, printSummary=TRUE)
#' summaryOutputList <- evaluateBounds(outputList)
#'
#' # $x$ & Nominal coverage (\Phi(x)) & True B in CI_x & Width-ratio: |Proposed width|/|DD| & Reverted to DD & Proposed Lower & Proposed Upper \\
#' # 0.00 & 0.5000 & TRUE & 0.4653 & FALSE & 0.6061 & 0.8101\\
#' # 0.25 & 0.5987 & TRUE & 0.5028 & FALSE & 0.5977 & 0.8182\\
#' # 0.50 & 0.6915 & TRUE & 0.5404 & FALSE & 0.5893 & 0.8262\\
#' # 0.75 & 0.7734 & TRUE & 0.5780 & FALSE & 0.5809 & 0.8343\\
#' # 1.00 & 0.8413 & TRUE & 0.6155 & FALSE & 0.5726 & 0.8424\\
#' # 1.25 & 0.8944 & TRUE & 0.6531 & FALSE & 0.5642 & 0.8505\\
#' # 1.50 & 0.9332 & TRUE & 0.6906 & FALSE & 0.5558 & 0.8586\\
#' # 1.75 & 0.9599 & TRUE & 0.7282 & FALSE & 0.5474 & 0.8666\\
#' # 2.00 & 0.9772 & TRUE & 0.7657 & FALSE & 0.5390 & 0.8747\\
#'
#' # For example, CI_0.5 (0.5893336 0.8262426) corresponds to c(summaryOutputList$CI_x_lower[3], summaryOutputList$CI_x_upper[3])
#' }
#'
#' @export
#'
evaluateBounds <- function(outputListFromGenerateBounds) {

  if (!is.null(outputListFromGenerateBounds$bd)) {
    b <- outputListFromGenerateBounds$bd
  } else {
    b <- NULL
  }
  bl1 <- outputListFromGenerateBounds$hbdl0
  bu1 <- outputListFromGenerateBounds$hbdu0
  cl <- outputListFromGenerateBounds$cil
  cu <- outputListFromGenerateBounds$cir
  dl <- outputListFromGenerateBounds$bdl
  du <- outputListFromGenerateBounds$bdu

  p <- outputListFromGenerateBounds$p
  l <- outputListFromGenerateBounds$l
  u <- outputListFromGenerateBounds$u

  sl <- bl1-cl
  su <- cu-bu1

  ## Loop through confidence levels, generating summary statistics (capture, width-ratios, etc.)
  xs=c(0.00,0.25,0.50,0.75, 1.00, 1.25, 1.50, 1.75, 2.00)
  isSelected <- numeric(length(xs))  # 0 if rejected by the heuristic (i.e., revert to DD)
  widthRatio <- numeric(length(xs))
  conf <- numeric(length(xs))

  proposedLowerBounds <- numeric(length(xs)) # lower bound across x in C_x
  proposedUpperBounds <- numeric(length(xs)) # upper bound across x in C_x

  if (!is.null(b)) {
    isCaptured <- numeric(length(xs))
  }
  for (j in 1:(length(xs) )) {

    ex=xs[j]

    bu=bu1+ex*su
    bl=bl1-ex*sl
    # intersect with DD:
    bbl = pmax(bl,dl)
    bbu = pmin(bu,du)

    if (!is.null(b)) { # if the true district B is provided, check if it appears in the proposed bound
      isCaptured[j] <- 1*(b<  bbu )*1*(b> bbl )
    }

    # heuristic for selection: check if CI_0 intersected with DD is flipped
    sel=1-1*(pmin(bu1,du)<pmax(bl1,dl))
    ddWidths <-du-dl

    if (sel == 0) {
      widthRatio[j] <- 1.0  # if not selected, set the width-ratio to 1 (i.e., using the DD bounds)
      proposedWidths <- ddWidths
      proposedLower <- dl
      proposedUpper <- du
    } else {
      widthRatio[j] <- (bbu-bbl)/(du-dl)
      proposedWidths <- bbu-bbl
      proposedLower <- bbl
      proposedUpper <- bbu
    }

    isSelected[j] <- sel
    conf[j] <- 1-pnorm(-ex)

    proposedLowerBounds[j] <- proposedLower
    proposedUpperBounds[j] <- proposedUpper

  }

  if (!is.null(b)) {
    cat("$x$ & Nominal coverage (\\Phi(x)) & True B in CI_x & Width-ratio: |Proposed width|/|DD| & Reverted to DD & Proposed Lower & Proposed Upper \\\\ \n")
    for (printIndex in 1:length(xs)) {
      cat(sprintf("%s & %s & %s & %s & %s & %s & %s\\\\ \n", format(round(xs[printIndex], 4), nsmall=2),
                  format(round(conf[printIndex], 4), nsmall=4),
                  isCaptured[printIndex]==1,
                  format(round(widthRatio[printIndex], 4), nsmall=4),
                  isSelected[printIndex]==0,
                  format(round(proposedLowerBounds[printIndex], 4), nsmall=4),
                  format(round(proposedUpperBounds[printIndex], 4), nsmall=4)))
    }
  } else {
    cat("$x$ & Nominal coverage (\\Phi(x)) & Width-ratio: |Proposed width|/|DD| & Reverted to DD & Proposed Lower & Proposed Upper \\\\ \n")
    for (printIndex in 1:length(xs)) {
      cat(sprintf("%s & %s & %s & %s & %s & %s\\\\ \n", format(round(xs[printIndex], 4), nsmall=2),
                  format(round(conf[printIndex], 4), nsmall=4),
                  format(round(widthRatio[printIndex], 4), nsmall=4),
                  isSelected[printIndex]==0,
                  format(round(proposedLowerBounds[printIndex], 4), nsmall=4),
                  format(round(proposedUpperBounds[printIndex], 4), nsmall=4)))
    }
  }
  summaryOutputList <- list()
  summaryOutputList[["x_for_x_in_CI_x"]] <- xs
  summaryOutputList[["CI_x_lower"]] <- proposedLowerBounds
  summaryOutputList[["CI_x_upper"]] <- proposedUpperBounds
  summaryOutputList[["CI_x_isSelected"]] <- isSelected == 1
  summaryOutputList[["CI_x_widthRatio"]] <- widthRatio
  summaryOutputList[["CI_x_nominalCoverage"]] <- conf
  if (!is.null(b)) {
    summaryOutputList[["CI_x_truthCaptured"]] <- isCaptured == 1
  }
  return(summaryOutputList)
}


#' Compute and evaluate bounds according to Jiang et al. 2019, illustrating usage.
#'
#' \code{bounds()} calculates district-level bounds across varying coverage probabilities, after applying the heuristics presented in Jiang et al. 2019. This is a
#' simple wrapper around calling generateBounds() followed by evaluateBounds(). Here, the returned object only contains the CI_0.5 bounds.
#'
#' @param x Numeric (double-precision) vector. Contains the proportion of variable X in each precinct (or analagous geographic unit)
#' @param t Numeric (double-precision) vector. Contains the proportion of variable T in each precinct (or analagous geographic unit)
#' @param n Numeric (double-precision) vector. Contains the number of elements (people/households/etc.) in each precinct (or analagous geographic unit)
#' @param trueBetaB Numeric (double-precision) vector. Contains the true conditional values (beta_i) in each precinct (or analagous geographic unit).
#'  Optional. Default NULL.
#'
#' @return List object with the CI_0.5 bounds:
#'
#'  CI_0.5_lower CI_0.5 lower bound
#'
#'  CI_0.5_upper CI_0.5 upper bound
#'
#'  CI_0.5_isSelected If FALSE, proposed bound was not rejected by the heuristic (if TRUE, bounds are reverted to the DD bounds)
#'
#'  CI_0.5_widthRatio |CI_x|/|DD|
#'
#'  CI_0.5_nominalCoverage Nominal coverage (1-pnorm(-0.5))
#'
#'
#'  Optional: CI_0.5_truthCaptured If true district Beta is provided as an argument to bounds(), then this variable contains a boolean for whether or not the
#'   true value was captured in the proposed CI_0.5.
#'
#' @examples
#' \dontrun{
#' library("MASS")
#' library("eco")
#' data("census")
#' inputDataSet <- census
#' x <- inputDataSet$X
#' t <- inputDataSet$Y
#' n <- inputDataSet$N
#' trueBetaB <- inputDataSet$W1
#' outputList <- bounds(x, t, n, trueBetaB=trueBetaB)
#' print(outputList)
#'
#' # > print(outputList)
#' # $CI_0.5_lower
#' # [1] 0.5893336
#' #
#' # $CI_0.5_upper
#' # [1] 0.8262426
#' #
#' # $CI_0.5_isSelected
#' # [1] TRUE
#' #
#' # $CI_0.5_widthRatio
#' # [1] 0.5404046
#' #
#' # $CI_0.5_nominalCoverage
#' # [1] 0.6914625
#' #
#' # $CI_0.5_truthCaptured
#' # [1] TRUE
#' }
#'
#' @export
#'
bounds <- function(x, t, n, trueBetaB=NULL) {
  outputList <- generateBounds(x, t, n, trueBetaB, useXRangeOffset=TRUE, returnAdditionalStats=FALSE, printSummary=FALSE)
  fullSummaryOutputList <- evaluateBounds(outputList)

  summaryOutputList <- list()
  summaryOutputList[["CI_0.5_lower"]] <- fullSummaryOutputList$CI_x_lower[3]
  summaryOutputList[["CI_0.5_upper"]] <- fullSummaryOutputList$CI_x_upper[3]
  summaryOutputList[["CI_0.5_isSelected"]] <- fullSummaryOutputList$CI_x_isSelected[3]
  summaryOutputList[["CI_0.5_widthRatio"]] <- fullSummaryOutputList$CI_x_widthRatio[3]
  summaryOutputList[["CI_0.5_nominalCoverage"]] <- fullSummaryOutputList$CI_x_nominalCoverage[3]
  if (!is.null(trueBetaB)) {
    summaryOutputList[["CI_0.5_truthCaptured"]] <- fullSummaryOutputList$CI_x_truthCaptured[3]
  }
  return(summaryOutputList)
}
