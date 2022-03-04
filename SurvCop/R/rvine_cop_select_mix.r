#' R-vine copula selection for mixed continuous/discrete variables
#'
#' @description
#' Bivariate copula selection for mixed continuous/discrete variables
#'
#' @param data An \eqn{n}-by-\eqn{d} data matrix with uniform margins.
#' The last column corresponds the response variable
#' and the first \eqn{d-1} columns correspond to the explanatory variables.
#' @param is_disc A boolean vector of length \eqn{d} indicating
#' whether each variable is discrete.
#' @param familyset A vector of bivariate copula families to select from;
#' this uses the family codes in VineCopula.
#' @param vine_array A \eqn{d}-by-\eqn{d} matrix representing the vine array.
#' The last column corresponds to the response variable.
#' @param selectioncrit The criterion for pair-copula selection.
#' Possible choices: selectioncrit = "AIC" (default), "BIC", or "logLik".
#' @param trunc_level The level of truncation.
#' @param cores integer; if \code{cores > 1}, estimation will be parallelized
#' within each tree (using \code{\link[foreach]{foreach}}).
#' @param debug Debug flag; default is FALSE
#'
#' @return A \code{VineCopula::RVineMatrix} object
#'
#' @details
#' Adapted from \code{VineCopula::RVineCopSelect}
#' 
#' @import igraph
rvine_cop_select_mix <- function(data, is_disc,
           familyset = NA,
           vine_array,
           selectioncrit = "AIC",
           trunc_level = NA,
           cores = 1,
           debug = FALSE) {
  vine_array <- ToLowerTri(vine_array)
  
  ## sanity checks
  if (!(selectioncrit %in% c("AIC", "BIC", "logLik")))
    stop("Selection criterion not implemented.")
  
  d <- n <- ncol(data)
  stopifnot(length(is_disc) == d)
  
  N <- nrow(data)
  ## set variable names and trunc_level
  varnames <- colnames(data)
  if (is.na(trunc_level))
    trunc_level <- d
  
  ## adjust familyset
  types <- familyset
  if (trunc_level == 0)
    types <- 0
  
  ## reorder matrix to natural order
  M <- vine_array
  Mold <- M
  o <- diag(M)
  M <- reorderRVineMatrix(M)
  data <- data[, o[length(o):1]]
  is_disc <- is_disc[o[length(o):1]]
  
  # For discrete variable j, calculate F^-_j
  dataMinus <- data
  listSupport <- NULL
  if (any(is_disc)) {
    obj <- PlusToMinus(dataMinus[, is_disc])
    dataMinus[, is_disc] <- obj$dataMinus
    listSupport <- obj$listSupport
  }
  
  ## create matrices required for selection of h-functions
  MaxMat <- createMaxMat(M)
  CondDistr <- neededCondDistr(M)
  
  ## create objects for results
  Types   <- matrix(0, d, d)
  Params  <- matrix(0, d, d)
  Params2 <- matrix(0, d, d)
  LogLik  <- matrix(0, d, d)
  Ses     <- matrix(0, d, d)
  Se2s    <- matrix(0, d, d)
  emptaus <- matrix(0, d, d)
  pvals   <- matrix(0, d, d)
  nobs    <- matrix(0, d, d)
  V <- list()
  V$directPlus <- array(0, dim = c(d, N))
  V$directMinus <- array(0, dim = c(d, N))
  V$indirectPlus <- array(0, dim = c(d, N))
  V$indirectMinus <- array(0, dim = c(d, N))
  V$isDiscDirect <- rep(F, d)
  V$isDiscIndirect <- rep(F, d)
  
  V$directPlus <- t(data[, d:1])
  V$directMinus <- t(dataMinus[, d:1])
  V$isDiscDirect <- is_disc[d:1]
  
  ## register parallel backend
  if (cores != 1 | is.na(cores)) {
    if (is.na(cores))
      cores <- max(1, detectCores() - 1)
    if (cores > 1) {
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      on.exit(try(stopCluster(), silent = TRUE))
      on.exit(try(closeAllConnections(), silent = TRUE)
              , add = TRUE)
    }
  }
  
  ## loop over all trees and pair-copulas
  warn <- NULL
  for (k in d:2) {
    doEst <- function(i) {
      if (k > i) {
        ## get pseudo-observaions
        m <- MaxMat[k, i]
        
        ## zr1 are the pseudo-observation of M[i, i] given M[(k+1):d, i]
        ## zr2 are the pseudo-observation of M[k, i] given M[(k+1):d, i]
        
        # zr1 <- V$direct[i, ]
        zr1Plus <- V$directPlus[i, ]
        zr1Minus <- V$directMinus[i, ]
        zr1IsDisc <- V$isDiscDirect[i]
        
        # zr2 <- if (m == M[k, i]) {
        #   V$direct[(d - m + 1), ]
        # } else {
        #   V$indirect[(d - m + 1), ]
        # }
        
        if (m == M[k, i]) {
          zr2Plus <- V$directPlus[(d - m + 1), ]
          zr2Minus <- V$directMinus[(d - m + 1), ]
          zr2IsDisc <- V$isDiscDirect[(d - m + 1)]
        } else {
          zr2Plus <- V$indirectPlus[(d - m + 1), ]
          zr2Minus <- V$indirectMinus[(d - m + 1), ]
          zr2IsDisc <- V$isDiscIndirect[(d - m + 1)]
        }
        
        # select pair-copula
        corr <- cor(zr1Plus, zr2Plus)
        
        # TODO: This is a quick fix, see email on Feb 4, 2020, 6:32 PM
        # > Perhaps one thing we can do is to replace these 90 and 270 rotations
        # > by a Gaussian copula with the same value of Kendall's tau.
        if (corr < 0)
          familyset <- 1
        
        if (trunc_level <= (d - k))
          familyset <- 0
        
        cfit <- suppressWarnings(
          BiCopSelectMix(
            zr2Plus,
            zr2Minus,
            zr1Plus,
            zr1Minus,
            zr2IsDisc,
            zr1IsDisc,
            familyset,
            selectioncrit,
            debug = debug
          )
        )
        
        ## transform data to pseudo-oberstavions in next tree
        directPlus <-
          directMinus <- indirectPlus <- indirectMinus <- NULL
        isDiscDirect <- isDiscIndirect <- NULL
        if (CondDistr$direct[k - 1, i]) {
          # direct <- suppressWarnings(BiCopHfunc1(zr2,
          #                                        zr1,
          #                                        cfit,
          #                                        check.pars = FALSE))
          if (!zr2IsDisc) {
            directPlus <- suppressWarnings(BiCopHfunc1(zr2Plus,
                                                       zr1Plus,
                                                       cfit,
                                                       check.pars = FALSE))
            directMinus <- suppressWarnings(BiCopHfunc1(zr2Plus,
                                                        zr1Minus,
                                                        cfit,
                                                        check.pars = FALSE))
          } else {
            directPlus = (
              BiCopCDF(zr2Plus, zr1Plus, cfit, check.pars = FALSE) -
                BiCopCDF(zr2Minus, zr1Plus, cfit, check.pars = FALSE)
            ) /
              (zr2Plus - zr2Minus)
            directMinus = (
              BiCopCDF(zr2Plus, zr1Minus, cfit, check.pars = FALSE) -
                BiCopCDF(zr2Minus, zr1Minus, cfit, check.pars = FALSE)
            ) /
              (zr2Plus - zr2Minus)
          }
          isDiscDirect <- zr1IsDisc
        }
        if (CondDistr$indirect[k - 1, i]) {
          # indirect <- suppressWarnings(BiCopHfunc2(zr2,
          #                                          zr1,
          #                                          cfit,
          #                                          check.pars = FALSE))
          if (!zr1IsDisc) {
            indirectPlus <- suppressWarnings(BiCopHfunc2(zr2Plus,
                                                         zr1Plus,
                                                         cfit,
                                                         check.pars = FALSE))
            indirectMinus <-
              suppressWarnings(BiCopHfunc2(zr2Minus,
                                           zr1Plus,
                                           cfit,
                                           check.pars = FALSE))
          } else {
            indirectPlus <-
              (
                BiCopCDF(zr2Plus, zr1Plus, cfit, check.pars = FALSE) -
                  BiCopCDF(zr2Plus, zr1Minus, cfit, check.pars = FALSE)
              ) /
              (zr1Plus - zr1Minus)
            indirectMinus <-
              (
                BiCopCDF(zr2Minus, zr1Plus, cfit, check.pars = FALSE) -
                  BiCopCDF(zr2Minus, zr1Minus, cfit, check.pars = FALSE)
              ) /
              (zr1Plus - zr1Minus)
          }
          isDiscIndirect <- zr2IsDisc
        }
        
        ## return results
        return(
          list(
            directPlus = directPlus,
            directMinus = directMinus,
            isDiscDirect = isDiscDirect,
            indirectPlus = indirectPlus,
            indirectMinus = indirectMinus,
            isDiscIndirect = isDiscIndirect,
            cfit = cfit,
            warn = warn
          )
        )
        
      } else {
        return(list(cfit = BiCop(0, 0)))
      }
    }
    
    ## run pair-copula selection for tree k
    res.k <- if (cores > 1) {
      foreach(
        i = 1:(k - 1),
        .packages = c("VineCopula"),
        .export = "familyset"
      ) %dopar% doEst(i)
    } else {
      lapply(1:(k - 1), doEst)
    }
    
    for (i in 1:(k - 1)) {
      ## store info about selected pair-copula in matrices
      Types[k, i]   <- res.k[[i]]$cfit$family
      Params[k, i]  <- res.k[[i]]$cfit$par
      Params2[k, i] <- res.k[[i]]$cfit$par2
      if (debug) {
        print(res.k[[i]]$cfit)
        cat("k=", k, "i=", i, "nllk=", res.k[[i]]$cfit$logLik, "\n")
      }
      LogLik[k, i]  <- res.k[[i]]$cfit$logLik
      
      if (!is.null(res.k[[i]]$directPlus)) {
        V$directPlus[i, ] <- res.k[[i]]$directPlus
        V$directMinus[i, ] <- res.k[[i]]$directMinus
        V$isDiscDirect[i] <- res.k[[i]]$isDiscDirect
      }
      
      if (!is.null(res.k[[i]]$indirectPlus)) {
        V$indirectPlus[i, ] <- res.k[[i]]$indirectPlus
        V$indirectMinus[i, ] <- res.k[[i]]$indirectMinus
        V$isDiscIndirect[i] <- res.k[[i]]$isDiscIndirect
      }
    } # end i = 1:(d-1)
    
    V$directPlus <- clip_to_unit(V$directPlus)
    V$directMinus <- clip_to_unit(V$directMinus)
    V$indirectPlus <- clip_to_unit(V$indirectPlus)
    V$indirectMinus <- clip_to_unit(V$indirectMinus)
    
  } # end k = d:2
  
  if (!is.null(warn))
    warning(" In ", args$call[1], ": ", warn, call. = FALSE)
  
  ## store results in RVineMatrix object
  .RVM <- RVineMatrix(
    Mold,
    family = Types,
    par = Params,
    par2 = Params2,
    names = varnames
  )
  
  .RVM$isDiscPermutated <- is_disc
  .RVM$listSupportPermutated <- listSupport
  
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
  
  revo <- sapply(1:d, function(i)
    which(o[length(o):1] == i))
  like <- suppressWarnings(RVineLogLik(data[, revo], .RVM))
  .RVM$pair.logLik <- LogLik
  npar <- sum(.RVM$family %in% allfams[onepar], na.rm = TRUE) +
    2 * sum(.RVM$family %in% allfams[twopar], na.rm = TRUE)
  npar_pair <- (.RVM$family %in% allfams[onepar]) +
    2 * (.RVM$family %in% allfams[twopar])
  .RVM$pair.AIC <- -2 * .RVM$pair.logLik + 2 * npar_pair
  .RVM$pair.BIC <- -2 * .RVM$pair.logLik + log(N) * npar_pair
  .RVM$logLik <- sum(.RVM$pair.logLik)
  .RVM$AIC    <- sum(.RVM$pair.AIC)
  .RVM$BIC    <- sum(.RVM$pair.BIC)
  
  .RVM$uscorePlus <- data
  .RVM$uscoreMinus <- dataMinus
  
  ## free memory and return final object
  rm(list = ls())
  .RVM
}

################################
# Code from VineCopula package #
################################

reorderRVineMatrix <- function(Matrix) {
  oldOrder <- diag(Matrix)
  
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)
  
  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- nrow(Matrix) - i + 1
  }
  
  return(Matrix)
}


createMaxMat <- function(Matrix) {
  if (dim(Matrix)[1] != dim(Matrix)[2])
    stop("Structure matrix has to be quadratic.")
  
  MaxMat <- reorderRVineMatrix(Matrix)
  
  n <- nrow(MaxMat)
  
  for (j in 1:(n - 1)) {
    for (i in (n - 1):j) {
      MaxMat[i, j] <- max(MaxMat[i:(i + 1), j])
    }
  }
  
  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0
  
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]
  
  for (i in 1:n) {
    MaxMat[tMaxMat == i] <- oldSort[i]
  }
  
  return(MaxMat)
}

neededCondDistr <- function(Vine) {
  if (dim(Vine)[1] != dim(Vine)[2])
    stop("Structure matrix has to be quadratic.")
  
  Vine <- reorderRVineMatrix(Vine)
  
  MaxMat <- createMaxMat(Vine)
  
  d <- nrow(Vine)
  
  M <- list()
  M$direct <- matrix(FALSE, d, d)
  M$indirect <- matrix(FALSE, d, d)
  
  M$direct[2:d, 1] <- TRUE
  
  for (i in 2:(d - 1)) {
    v <- d - i + 1
    
    bw <- as.matrix(MaxMat[i:d, 1:(i - 1)]) == v
    
    direct <- Vine[i:d, 1:(i - 1)] == v
    
    M$indirect[i:d, i] <- apply(as.matrix(bw & (!direct)), 1, any)
    
    M$direct[i:d, i] <- TRUE
    
    M$direct[i, i] <-
      any(as.matrix(bw)[1, ] & as.matrix(direct)[1, ])
  }
  
  return(M)
}

ToLowerTri <- function(x) {
  ## only change matrix if not already lower triagonal
  if (all(x[lower.tri(x)] == 0)) {
    x[nrow(x):1, ncol(x):1]
  } else {
    x
  }
}
