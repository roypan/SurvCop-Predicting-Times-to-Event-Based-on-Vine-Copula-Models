#' negative log-likelihood of vine copula model for survival response
#' 
#' @description
#' negative log-likelihood of a vine copula model for a survival response (the last column of the data matrix) and mixed discrete/continuous explanatory variables (the other columns of the data matrix)
#' 
#' @param paramvec vector of length (d-1) if all the copula families associated with the survival response variable are one-parameter families; vector of length (d-1)+ otherwise
#' @param param1_complete d-1 by d-1 matrix for the first parameter of the copula models for the predictors
#' @param param2_complete d-1 by d-1 matrix for the second parameter of the copula models for the predictors
#' @param fam d by d matrix for the copula families coded by the VineCopula package
#' @param A d by d matrix for the copula structure
#' @param data_cdf_plus n by (d-1) matrix for CDF from the right
#' @param data_cdf_minus n by (d-1) matrix for CDF from the left
#' @param censor_status vector of length n representing if the response variable is censored or not for each observation
#' @param is_disc vector of length d representing if each variable is discrete or not
#'
#' @return the negative log-likelihood of the specified vine copula model for a survival response
#'
#' @import VineCopula
#' @export

survival_nllk2 = function(paramvec, param1_complete, param2_complete, fam, A, data_cdf_plus, data_cdf_minus, censor_status, is_disc) {
    d = ncol(A)
    n = nrow(data_cdf_plus)
    
    # permute all arguments based on the order of A
    A_order = diag(A)
    A = varrayperm(A, order(A_order))
    
    data_cdf_plus = data_cdf_plus[, A_order]
    data_cdf_minus = data_cdf_minus[, A_order]
    is_disc = is_disc[A_order]
    
    paramvec = parse_parameter2(paramvec, fam[1:(d-1), d], d)
    param1_vec = paramvec[1:(d-1)]
    param2_vec = paramvec[d:length(paramvec)]
    
    parammat1 = matrix(0, nrow = d, ncol = d)
    parammat1[1:(d-1), 1:(d-1)] = param1_complete
    parammat1[1:(d-1), d] = param1_vec
    
    parammat2 = matrix(0, nrow = d, ncol = d)
    parammat2[1:(d-1), 1:(d-1)] = param2_complete
    parammat2[1:(d-1), d] = param2_vec
    
    out= .Fortran("survivalnllk", as.integer(n), as.integer(d), as.double(parammat1), as.double(parammat2), as.integer(fam), as.integer(A), as.double(data_cdf_plus), as.double(data_cdf_minus), as.logical(censor_status), as.logical(is_disc), nllk=as.double(0))
    
    out$nllk
}