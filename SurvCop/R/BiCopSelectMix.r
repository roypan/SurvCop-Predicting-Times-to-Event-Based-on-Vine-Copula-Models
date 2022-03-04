#' Bivariate Copula Selection with Mixed continuous and discrete variables
#'
#' @description
#' Selection and Maximum Likelihood Estimation of Bivariate Copula Families with
#' Continuous-Discrete Mixed Variables.
#'
#' @param u1Plus uscore1+ or cdf at right limit.
#' @param u1Minus uscore1- or cdf at left limit.
#' @param u2Plus uscore2+ or cdf at right limit.
#' @param u2Minus uscore2- or cdf at left limit.
#' @param isDisc1 Boolean. TRUE if the first variable is discrete.
#' @param isDisc2 Boolean. TRUE if the second variable is discrete.
#' @param familyset Vector of bivariate copula families to select from.
#' The vector has to include at least one bivariate copula
#' family that allows for positive and one that allows for negative dependence.
#' If \code{familyset = NA} (default), selection among all possible families is
#' performed. If a vector of negative numbers is provided, selection among all
#' but \code{abs(familyset)} families is performed. Coding of bivariate copula
#' families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) \cr
#' @param selectioncrit Character indicating the criterion for bivariate copula
#' selection. Possible choices: \code{selectioncrit = "AIC"} (default),
#' \code{"BIC"}, or \code{"logLik"}.
#' @param max.df Numeric; upper bound for the estimation of the degrees of
#' freedom parameter of the t-copula (default: \code{max.df = 30}).
#' @param max.BB List; upper bounds for the estimation of the two parameters
#' (in absolute values) of the BB1, BB6, BB7 and BB8 copulas \cr (default:
#' \code{max.BB = list(BB1=c(5,6),BB6=c(6,6),BB7=c(5,6),BB8=c(6,1))}).
#' @param debug Debug flag; default is FALSE
#'
#' @details 
#' See \link[VineCopula]{BiCopSelect}.
#' Part of the code is adapted from \link[VineCopula]{BiCopSelect}.
#' Need to be careful what rotation means if bivariate copula is not
#' permutation symmetric. \cr
#' \eqn{(U_1,U_2)} positively dependent implies the copulas of
#' \eqn{(1-U_1,U_2)} and \eqn{(U_1,1-U_2)} have negative dependence. 
#' Refer to these as 1-reflected, 2-reflected;
#' and \eqn{(1-U_1,1-U_2)} as reflected or (12)-reflected.\cr
#' Counterclockwise rotation \eqn{(U_1,U_2)} in the unit square about the centre (1/2,1/2) leads to the following.\cr
#' Rotation of \eqn{\pi/2}:   \eqn{(U_1,U_2)} to \eqn{(1-U_2,U_1)}.\cr
#' Rotation of \eqn{2\pi/2}:  \eqn{(U_1,U_2)} to \eqn{(1-U_1,1-U_2)}.\cr
#' Rotation of \eqn{3\pi/2}:  \eqn{(U_1,U_2)} to \eqn{(U_2,1-U_1)}.\cr
#' Check if rotation in VineCopula is clockwise or counterclockwise and update here.
#'
#' @return A \code{BiCop} object.
#'
BiCopSelectMix <- function (u1Plus, u1Minus, u2Plus, u2Minus,
    isDisc1, isDisc2,
    familyset,
    selectioncrit = "AIC",
    max.df = 30,
    max.BB = list( BB1 = c(5,6), BB6 = c(6,6), BB7 = c(5,6), BB8 = c(6,1)),
    debug = FALSE) {
  # If both u1 and u2 are continuous, use VineCopula::BiCopSelect.
  if (!isDisc1 && !isDisc2) {
    return(VineCopula::BiCopSelect(u1Plus, u2Plus, familyset, selectioncrit, rotations = F))
  }
  
  ## (Internal function) Given the bivariate copula family,
  ## fit the parameters using MLE.
  ##
  ## @param family Copula family: an integer.
  ## @param type Variable types: \code{"DiscDisc"} for both discrete variables;
  ## \code{"DiscCont"} for discrete \code{u1} and continuous \code{u2};
  ## \code{"ContDisc"} for continuous \code{u1} and discrete \code{u2}.
  ##
  ## @return
  ## \item{nllk}{Negative log-likelihood.}
  ## \item{par1}{The estimated first parameter.}
  ## \item{par2}{The estimated second parameter.
  ## Zero for one-parameter copula families.}
  BiCopMLE <- function (family, type) {
    stopifnot(type %in% c("DiscDisc", "DiscCont", "ContDisc"))
    
    NegLogLik <- function (par) {
      # In case of error, return 1e7.
      tryCatch({
        if (family %in% allfams[onepar]) {
          par1 <- par
          par2 <- 0
        } else if (family %in% allfams[twopar]) {
          par1 <- par[1]
          par2 <- par[2]
        } else {
          stop("Copula family not supported.")
        }
        
        # See "Dependence Modeling with Copulas" pp. 115-116 for the definition
        # of \tilde{c}.
        if (type == "DiscDisc") {
          # C(F_1^+, F_2^+)
          cdfPlusPlus <- BiCopCDF(
            u1Plus,
            u2Plus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          # C(F_1^+, F_2^-)
          cdfPlusMinus <- BiCopCDF(
            u1Plus,
            u2Minus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          # C(F_1^-, F_2^+)
          cdfMinusPlus <- BiCopCDF(
            u1Minus,
            u2Plus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          # C(F_1^-, F_2^-)
          cdfMinusMinus <- BiCopCDF(
            u1Minus,
            u2Minus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          
          cTilde <-
            (cdfPlusPlus - cdfPlusMinus - cdfMinusPlus + cdfMinusMinus) /
            (u1Plus - u1Minus) / (u2Plus - u2Minus)
        } else if (type == "DiscCont") {
          # C_{1|2}(F_1^+ | F_2)
          cdfPlus <- BiCopHfunc2(
            u1Plus,
            u2Plus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          # C_{1|2}(F_1^- | F_2)
          cdfMinus <- BiCopHfunc2(
            u1Minus,
            u2Plus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          
          cTilde <- (cdfPlus - cdfMinus) / (u1Plus - u1Minus)
        } else if (type == "ContDisc") {
          # C_{2|1}(F_2^+ | F_1)
          cdfPlus <- BiCopHfunc1(
            u1Plus,
            u2Plus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          # C_{2|1}(F_2^- | F_1)
          cdfMinus <- BiCopHfunc1(
            u1Plus,
            u2Minus,
            family = rep(family, n),
            par = rep(par1, n),
            par2 = rep(par2, n)
          )
          
          cTilde <- (cdfPlus - cdfMinus) / (u2Plus - u2Minus)
        }
        
        cTilde[cTilde < 1e-6] <- 1e-6
        return(-sum(log(cTilde)))
      }, error = function(err) {
        return(1e7)
      })
    }
    
    # # If Kendall's tau is greater than 0, 90/270 rotated copula will throw
    # # an error. In that case, return nllk = Inf.
    # # The method "itau" won't work for two parameter copula families, it will
    # # throw a warning; the warning is suppressed.
    # initialEst <- tryCatch(
    #   suppressWarnings(BiCopEst(u1Plus, u2Plus, family, method = "itau")),
    #   error = function(e) return(NA)
    # )
    # if (all(is.na(initialEst)))
    #   return(list(family = family, nllk = Inf, par1 = NaN, par2 = NaN))
    
    # upper bound and lower bound for the copula family
    bounds <- BiCopParamBounds(family, max.df, max.BB)
    nllk <- 0
    par1 <- 0
    par2 <- 0
    
    # use a (arbitrary) convex combination of the upper bound and lower bound as inital par.
    # cannot use the average because Frank copula do not allow parameter to be zero.
    initial_par <- bounds$low * 0.4 + bounds$up * 0.6
    
    if (family %in% allfams[onepar]) {
      optimObj <- optim(
        par = initial_par,
        #par = initialEst$par,
        fn = NegLogLik,
        method = "L-BFGS-B",
        lower = bounds$low,
        upper = bounds$up
      )
      nllk <- optimObj$value
      par1 <- optimObj$par
      par2 <- 0
    } else if (family %in% allfams[twopar]) {
      optimObj <- optim(
        par = initial_par,
        #par = c(initialEst$par, initialEst$par2),
        fn = NegLogLik,
        method = "L-BFGS-B",
        lower = bounds$low,
        upper = bounds$up
      )
      nllk <- optimObj$value
      par1 <- optimObj$par[1]
      par2 <- optimObj$par[2]
    }
    
    return(list(
      family = family,
      nllk = nllk,
      par1 = par1,
      par2 = par2
    ))
  }
  
  
  ## (Internal function) Interate through the list of \code{familyset} and find
  ## the copula family that minimizes the \code{selectioncrit}
  ##
  ## @return A \code{BiCop} object.
  BiCopSelectDisc <- function () {
    optiout <- list()
    nllks  <- rep(Inf, length(familyset))
    AICs <- rep(Inf, length(familyset))
    BICs <- rep(Inf, length(familyset))
    
    for (i in seq_along(familyset)) {
      optiout[[i]] <- BiCopMLE(familyset[i], type)
      nllks[i] <- optiout[[i]]$nllk
      npars <- ifelse(familyset[i] %in% allfams[onepar], 1, 2)
      if (familyset[i] == 0)
        npars <- 0
      AICs[i] <- 2 * nllks[i] + 2 * npars
      BICs[i] <- 2 * nllks[i] + log(n) * npars
      
      # if (debug)
      #   cat("==========\nCopula family:", familyset[i], "\nnllk:", nllks[i],
      #       "\npar1:", optiout[[i]]$par1,
      #       "\npar2:", optiout[[i]]$par2, "\n==========\n\n")
    }
    
    ## select the best fitting model
    sel <- switch(
      selectioncrit,
      "logLik" = which.min(nllks),
      "AIC"    = which.min(AICs),
      "BIC"    = which.min(BICs)
    )
    obj <- BiCop(optiout[[sel]]$family,
                 optiout[[sel]]$par1,
                 optiout[[sel]]$par2,
                 check.pars = FALSE)
    obj$logLik <- -nllks[sel]
    obj$AIC <- AICs[sel]
    obj$BIC <- BICs[sel]
    obj$AIClist <- cbind(familyset, AICs)
    return(obj)
  }
  
  ## =========================
  ## Main function
  
  
  stopifnot(length(u1Plus) == length(u2Plus))
  n <- length(u1Plus)
  
  if (isDisc1)
    stopifnot(all(u1Plus - u1Minus > 0))
  if (isDisc2)
    stopifnot(all(u2Plus - u2Minus > 0))
  
  ## sets of families
  ## Use allfams[onepar] and allfams[twopar]
  allfams <- c(0:10,
               13,
               14,
               16:20,
               23,
               24,
               26:30,
               33,
               34,
               36:40,
               104,
               114,
               124,
               134,
               204,
               214,
               224,
               234)
  tawns <- which(allfams > 100)
  onepar <-
    setdiff(which(allfams %% 10 %in% c(1, 3, 4, 5, 6)), tawns)
  twopar <- seq_along(allfams)[-c(1, onepar)]
  
  
  # u1[u1 < 1e-6] <- 1e-6
  # u2[u2 < 1e-6] <- 1e-6
  # u1[u1 > 1 - 1e-6] <- 1 - 1e-6
  # u2[u2 > 1 - 1e-6] <- 1 - 1e-6
  
  # use BiCopSelectDisc.
  if (isDisc1 && isDisc2)
    type <- "DiscDisc"
  else if (isDisc1 && !isDisc2)
    type <- "DiscCont"
  else if (!isDisc1 && isDisc2)
    type <- "ContDisc"
  else
    stop("Invalid input: isDisc1 and/or isDisc2.")
  
  return(BiCopSelectDisc())
  
}

## Internal function
## Return the upper and lower bounds for a copula family.
BiCopParamBounds <- function(family, max.df, max.BB) {
  if (family == 0) {
    low <- up <- 0
  } else if (family == 1) {
    low <- -0.9999
    up <- 0.9999
  } else if (family == 2) {
    low <- c(-0.9999, 2.0001)
    up <- c(0.9999, max.df)
  } else if (family %in% c(3, 13)) {
    low <- 1e-04
    up <- 100
  } else if (family %in% c(4, 14)) {
    low <- 1.0001
    up <- 50
  } else if (family %in% c(5)) {
    low <- -100
    up <- 100
  } else if (family %in% c(6, 16)) {
    low <- 1.0001
    up <- 50
  } else if (family %in% c(23, 33)) {
    up <- -1e-04
    low <- -100
  } else if (family %in% c(24, 34)) {
    up <- -1.0001
    low <- -100
  } else if (family %in% c(26, 36)) {
    up <- -1.0001
    low <- -50
  } else if (family == 7 || family == 17) {
    low <- c(0.001, 1.001)
    up <- pmin(max.BB$BB1, c(7, 7))
  } else if (family == 8 || family == 18) {
    low <- c(1.001, 1.001)
    up <- pmin(max.BB$BB6, c(6, 8))
  } else if (family == 9 | family == 19) {
    low <- c(1.001, 0.001)
    up <- pmin(max.BB$BB7, c(6, 75))
  } else if (family == 10 | family == 20) {
    low <- c(1.001, 0.001)
    up <- pmin(max.BB$BB8, c(8, 1))
  } else if (family == 27 | family == 37) {
    up <- c(-0.001, -1.001)
    low <- -pmin(max.BB$BB1, c(7, 7))
  } else if (family == 28 | family == 38) {
    up <- c(-1.001, -1.001)
    low <- -pmin(max.BB$BB6, c(6, 8))
  } else if (family == 29 | family == 39) {
    up <- c(-1.001, -0.001)
    low <- -pmin(max.BB$BB7, c(6, 75))
  } else if (family == 30 | family == 40) {
    up <- c(-1.001, -0.001)
    low <- -pmin(max.BB$BB8, c(8, 1))
  } else {
    stop("Copula family not supported.")
  }
  
  return(list(up = up, low = low))
}
