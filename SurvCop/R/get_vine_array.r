#' Get a vine array based on different optimization criteria
#' 
#' @description
#' Find a vine structure from a correlation matrix based on different optimization criteria
#'
#' @param rmat Correlation matrix.
#' @param n Sample size.
#' @param method Currently support "enumx", "MSTx", "MSTleaf", "forwardD", "forwardR". Default: "MSTx". See Details for explanations.
#' @param cfi_bd Lower bound of CFI (comparative fit index) to stop, default 0.95, set to 1.01 for no truncation.
#' @param debug Debug flag.
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{vine_array} - vine array to truncation level
#'  \item \code{trunc_level} - truncation level to reach cfi_bd
#'  \item \code{cfi} - CFI value
#'  \item \code{cfi_vector} - vector of CFI values
#'  \item \code{treeweight} - vector of tree weights by tree, level
#'  \item \code{fitval} - cumsum(treeweight)/\eqn{sum_{1:(d-1)} treeweight}
#'  \item \code{pcmat} - matrix of partial correlations by tree
#' }
#'

#' @details 
#' (a) enumx: for number of covariates <=7; enumeration as in Section 6.13 of Joe (2014).\cr
#' (b) MSTx: sequential minimal spanning trees by vine level; see Section 6.17 of Joe (2014) and Chang and Joe (2019).\cr
#' (c) MSTleaf: sequential minimal spanning trees with response variable y forced to be a leaf on each vine tree.\cr
#' (d) forwardD: Forward selection D-vine as in Kraus and Czado (2017) when restricted to gaussian pair-copulas.
#' That is, covariates are added sequentially, with forward selection based on partial correlations.\cr
#' (e) forwardR: Forward selection as in forwarD, but then the covariates are linked in an R-vine.
#' Because the aim is to have strong correlations and partial correlations in the low level trees, the heuristic is similar to MSTx.
#' See Section 5.4 of Herrmann (2019).
#' 
#' @references
#' 1. Kraus D. and Czado C., 2017. D-vine copula based quantile regression.
#' Comput. Statist. Data Anal. 110, 1-18. \cr
#' 2. Chang B and Joe H 2019. Prediction based on conditional distributions.
#' Comput. Statist. Data Anal. 139, 45-63. \cr
#' 3. Herrmann, J S (2019).
#' Regular vine copula based quantile regression.
#' MSc Thesis, Technische Universitaet Muenchen. \cr
#' 4. Joe H (2014). Dependence Models with Copulas. Chapman-Hall/CRC
#'
#' @examples
#' ##See the donkey data vigenette for a comparison of 4 of the methods.
#'
#' @export
#'
get_vine_array <- function(rmat, n, method="MSTx", cfi_bd=0.95, debug=FALSE) {
    stopifnot(nrow(rmat) == ncol(rmat))
    d <- nrow(rmat)
    rmat_x <- rmat[-d, -d]
    
    if (method == "MSTx") {
      vine_array_x <- gaussvine.mst(rmat_x, n, 1.01, debug, debug)$VineA
      vine_array <- helper_append_y(rmat, vine_array_x)
      
    } else if (method == "enumx") {
      for (trunc_level_x in 1:(d - 2)) {
        # Find the best `trunc_level_x`-truncated vine array `vine_array_x` for the covariates,
        # then append y to it. If the resulting vine_array has better CFI than `cfi_bd`, then
        # stop, otherwise increment `trunc_level_x`.
        bestvines <-
          gausstrvine_modified(trunc_level_x, rmat_x, iprint = debug)
        vine_array_x <-
          vnum2array(d - 1, bestvines$bnum[trunc_level_x])
        vine_array <- helper_append_y(rmat, vine_array_x)
        current_cfi <-
          get_cfi_from_vine_array(rmat, vine_array, n, trunc_level_x)
        if (current_cfi > cfi_bd)
          break
      }
      
    } else if (method == "MSTleaf") {
      output <-
        yleaf.mst(
          rmat = rmat,
          n = n,
          CFIbd = cfi_bd,
          iprint = debug,
          iprint2 = debug
        )
      return(parse_output(output))
      
    } else if (method == "forwardD") {
      output <- forwsel2Darray(rmat, iprint = debug)
      vine_array <- output$A
      
    } else if (method == "forwardR") {
      diagm <- forwsel2order(rmat, iprint = debug)
      vine_array <- xorder2varray(rmat, diagm$A, iprint = debug)
      
    } else {
      stop(paste0("Method", method, "is not supported.\n"))
    }
    
    output <- trvine.CFI(rmat, vine_array, n, cfi_bd, debug)
    parse_output(output)
  }


# this function is used in get_vine_array()
parse_output <- function(output) {
  # parse the output of trvine.CFI to be compatible with the output format of get_vine_array.
  list(
    vine_array = output$VineA,
    trunc_level = output$ntrunc,
    cfi = output$CFIv[output$ntrunc],
    cfi_vector = output$CFIv,
    treeweight = output$treeweight,
    fitval = output$fitval,
    pcmat = output$pcmat
  )
}

# This function is used by get_vine_array for options MSTx and enumx.
#' Appended vine array with response variable y
#' 
#' @description
#' Given a vine array for covariates, this helper function adds the response variable so that max partial correlations are achieved subject to the proximity condition for vines.
#'
#' @param rmat Correlation matrix of covariates and response; vander Waarden/ polychoric/ polyserial correlations.
#' @param vine_array_x Vine array of the explanatory variables.
#' @param trunc_level Truncation level
#'
#' @return Full vine array (one more row/column compared with input); the response variable is in the last column.
#'
helper_append_y <- function(rmat, vine_array_x, trunc_level = nrow(rmat)-1) {
    d <- nrow(rmat)
    last_col_varray <- rep(0, d)
    last_col_rmat <- rmat[1:(d - 1), d]
    
    # level 1 and d
    last_col_varray[1] <- which.max(last_col_rmat ^ 2)
    last_col_varray[d] <- d
    
    # other levels
    for (k in 2:(d - 2)) {
      if (k > trunc_level)
        # reached the truncation level
        break
      
      candidates <- rep(0, d - k)
      pcor_candidates <- rep(0, d - k)
      
      # a matrix with all conditioning/conditioned set
      condition.set <-
        rbind(vine_array_x[1:(k - 1), k:(d - 1)], diag(vine_array_x)[k:(d - 1)])
      
      for (i in 1:(d - k)) {
        current_set <- condition.set[, i]
        conditioning <- last_col_varray[1:(k - 1)]
        
        if (all(conditioning %in% current_set)) {
          candidates[i] <- setdiff(current_set, conditioning)
          pcor_candidates[i] <-
            partcor(rmat, conditioning, d, candidates[i])
        }
      }
      
      # at least one variable available
      stopifnot(sum(pcor_candidates ^ 2) > 0)
      
      last_col_varray[k] <-
        candidates[which.max(pcor_candidates ^ 2)]
    }
    
    if (trunc_level == d - 1)
      last_col_varray[d - 1] <-
      setdiff(1:(d - 1), last_col_varray[1:(d - 2)])
    
    vine_array <- rbind(vine_array_x, rep(0, d - 1))
    vine_array <- cbind(vine_array, last_col_varray)
    colnames(vine_array) <- NULL
    
    return(vine_array)
  }
