# vine functions from CopulaModel

# varray2M(A,iprint = F,str = "") : maximum array computed from vine array A,  
#     used in rvine_trunc_cond_quantile.R
#  (not exported)
# vnum2array(d,bnum = 0,iprint = F) : used in get_vine_array.R
#  (not exported)
# vstepb, d2b : used in vnum2array (not exported)

# varrayperm(A,iperm) : used in copreg.R, get_cfi.R
# Dvinearray(d,iNO = F) : used in varrayopt.R, rvine_trunc_cond_quantile.R
# cor2pcor.rvine(rr,A) : used in get_cfi.R

# varray2NO : not needed
# vbin2array : not needed
# varraycheck : not needed
# Cvinearray  : not needed

#============================================================

# A = dxd vine array of R-vine in standard order with 1:d on diagonal
# iprint = print flag for intermediate calculations
# str = string to describe vine for printing if iprint=T
# Output: list with two components
#  M = dxd array with m_{kj}= max a_{k1},..,a_{kj}
#    actually M could be put in lower triangle of A
#  icomp = dxd indicator array on whether back step [k,j] is needed
#    icomp[k-1,m_{kj}=1 if  a_{kj}<m_{kj} for k>=2
varray2M = function(A,iprint = F,str = "")
{
  d = ncol(A)
  d1 = d - 1
  M = A
  icomp = matrix(0,d,d)
  for (k in 2:d1)
  {
    for (j in (k + 1):d)
      M[k,j] = max(M[k - 1,j],A[k,j])
  }
  if (iprint) {
    cat("\n",str,"\n"); print(A); print(M)
  }
  for (k in 2:d1)
  {
    for (j in (k + 1):d)
    {
      if (A[k,j] < M[k,j])
        icomp[k - 1,M[k,j]] = 1
    }
  }
  if (iprint)
    print(icomp)
  list(mxarray = M,icomp = icomp)
}

#C= matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
#D= matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 3,2,1,4,0, 4,3,2,1,5), 5,5)
#B0=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 1,3,2,4,5), 5,5)
#B1=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,2,3,4,0, 2,1,3,4,5), 5,5)
#B2=matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 1,2,3,4,0, 1,2,3,4,5), 5,5)
#B3=matrix(c(1,0,0,0,0, 1,2,0,0,0, 1,2,3,0,0, 1,3,2,4,0, 2,1,3,4,5), 5,5)
#C6= matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 1,2,3,4,0,0, 1,2,3,4,5,0,
#   1,2,3,4,5,6),6,6)
#D6= matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 2,1,3,0,0,0, 3,2,1,4,0,0, 4,3,2,1,5,0,
#   5,4,3,2,1,6),6,6)

#d=5
#A=matrix(0,d,d)
#diag(A)=1:5
#A[1,2:5]=c(1,2,2,4)
#A[2,3:5]=c(1,1,2)
#A[3,4:5]=c(3,1)
#A[4,5]=3
#varray2M(A1,iprint=T,"A")
#varray2M(C ,iprint=T,"C ")
#varray2M(D ,iprint=T,"D ")
#varray2M(B0,iprint=T,"B0")
#varray2M(B1,iprint=T,"B1")
#varray2M(B2,iprint=T,"B2")
#varray2M(B3,iprint=T,"B3")
#varray2M(C6,iprint=T,"C6")
#varray2M(D6,iprint=T,"D6")

#============================================================


# Create dxd vine array based on a single number representation of the
#   binary matrix representation
# The number bnum is converted to the binary array
# d = integer>=3; if(d==3) bnum=0,
# bnum = integer between 0 and 2^((d-2)*(d-3)/2)-1 ;
#   2^((d-2)*(d-3)/2) is the number of d-dimensional vine arrays
#     in natural order.
#  bnum=0 is the D-vine, bnum=2^((d-2)*(d-3)/2)-1 is the C-vine
# iprint = print flag for intermediate steps
# To get bnum from a binary matrix representation bmat,
#   bvec as a binary vector starts increments from the last flexible position,
#   which is bmat[d-2,d];
# the (d-2)*(d-3)/2 positions in bmat: [2,4], [3,4],[3,5],...,[2,d],...,[d-2,d]
# Output: dxd vine array in natural order form
vnum2array = function(d,bnum = 0,iprint = F)
{
  dcase = (d - 2) * (d - 3) / 2
  dpow = 2 ^ dcase
  A = matrix(0,d,d); diag(A) = 1:d; A[1,3] = 1
  for (i in 2:d)
    A[i - 1,i] = i - 1
  if (d == 3)
    return(A)
  
  bnum = floor(bnum)
  if (bnum < 0 | bnum >= dpow)
    bnum = 0
  bvec = d2b(dcase,bnum)
  b = matrix(NA,d,d)
  b[1,] = 1
  diag(b) = 1
  for (i in 3:d)
    b[i - 1,i] = 1
  #for(i in 4:d) b[2:(i-2),i]=rbinom(i-3,1,0.5)
  ii = 0
  for (i in 4:d) {
    b[2:(i - 2),i] = bvec[(ii + 1):(ii + i - 3)]; ii = ii + i - 3
  }
  if (iprint)
    print(b)
  
  # column 4
  if (b[2,4] == 1)
    # C-vine for first 4 columns
  {
    A[1,4] = 1; A[2,4] = 2
  }
  else {
    A[1,4] = 2; A[2,4] = 1
  }  # D-vine for first 4 columns
  # columns 5 and higher
  if (d >= 5)
  {
    for (i in 5:d)
    {
      b0 = b[2:(i - 2),i]  # length i-3
      A[,i] = vstepb(b0,A,i)
    }
  }
  A
}

# This function is used for columns 5 or higher in vnum2array
# b0 = d-vector with length i-3,
# A =  vine array
# i = column from 4 to ncol(A)
#   A has dimension at least ixi
# Output: ith column of A based on the binary representation b0 for column i
vstepb = function(b0,A,i)
{
  itaken = rep(0,i)
  itaken[i] = 1
  itaken[i - 1] = 1
  b = c(1,b0,1,1)
  #ac=i-1  # active column
  ac = i - 2  # active column
  A1 = A
  A1[i,i] = i
  A1[i - 1,i] = i - 1
  #for(k in (i-1):1) # older version with ac=i-1
  for (k in (i - 2):1)
  {
    if (b[k] == 1)
    {
      tem = A1[ac,ac]; itaken[tem] = 1
      A1[k,i] = tem
      if (k > 1)
      {
        ac = max((1:i)[itaken == 0])
      }
    }
    else
    {
      tem = A1[k - 1,ac]; A1[k,i] = tem; itaken[tem] = 1
    }
  }
  #print(A1)
  A1[,i]
}


# decimal to binary vector
# d = dimension
# ii = integer ii  in 0 to 2^d-1
# Output: binary d-vector jj corresponding to integer ii  
d2b=function(d,ii)
{ s=ii; jj=rep(0,d)
  for(i in d:1)
  { jj[i]=s%%2; s=floor(s/2) }
  jj
}


#============================================================

#' Vine array with permuted indices
#' 
#' @description
#' Create vine array with permuted indices, useful because some
#' algorithms for vines assume ordered 1:d on the diagonal.
#' 
#' @param A dxd vine array
#' @param perm permutation of 1:d
#' @return vine array after permutation of indices
#' 
#' @examples
#' D5= matrix(c(1,0,0,0,0, 1,2,0,0,0, 2,1,3,0,0, 3,2,1,4,0, 4,3,2,1,5), 5,5)
#' D5perm= matrix(c(5,0,0,0,0, 5,4,0,0,0, 4,5,3,0,0, 3,4,5,2,0, 2,3,4,5,1), 5,5)
#' varrayperm(D5,perm=c(5,4,3,2,1))
#' 
#' @export
#'
varrayperm = function(A, perm)
{
  A2 = A
  d = ncol(A)
  for (i in 1:d)
  {
    for (j in i:d)
      A2[i,j] = perm[A[i,j]]
  }
  A2
}

#' Vine array for D-vine
#'
#' @description
#' Vine array for D-vine in standard form or normal order (latter is only useful for enumerating vines)
#'
#' @param d dimension
#' @param iNO FALSE for standard form (default) or TRUE for natural order
#'
#' @return 
#' vine array for D-vine 
#' 
#' @examples
#' Dvinearray(6) 
#' Dvinearray(5) # same as 5x5 submatrix of above
#' 
#' @export
#'
Dvinearray = function(d, iNO=FALSE)
{
  D = matrix(0,d,d)
  diag(D) = 1:d
  if (iNO) {
    D = vnum2array(d,bnum = 0)
  }
  else {
    for (j in 2:d)
      D[1:(j - 1),j] = (j - 1):1
  }
  D
}


#============================================================

#' correlation matrix to partial correlation vine
#' 
#' @description
#' Convert a correlation matrix to a partial correlation representation based on a given vin array
#' 
#' @param rr dxd positive definite correlation matrix 
#' @param A dxd vine array with 1:d on diagonal (only upper triangle is used)
#' 
#' @return 
#' List with (i) $pctree, where partial correlations for tree \eqn{\ell} are in row \eqn{\ell}, and
#'  (ii) $pcmat, \eqn{\rho_{jk;S}} in position (j,k) where S is conditioning vector of edge [j,k].
#'
#' @examples
#' rr=toeplitz(c(1,.5,.25,.125,.0625,.05))
#' D6=Dvinearray(6)
#' cor2pcor.rvine(rr,D6)  # many 0s
#' A6=matrix(c(1,0,0,0,0,0, 1,2,0,0,0,0, 1,2,3,0,0,0, 2,1,3,4,0,0, 
#'          1,3,2,4,5,0, 2,4,1,3,5,6), 6,6)  
#' cor2pcor.rvine(rr,A6)
#' 
#' @export
#' 
cor2pcor.rvine=function(rr, A)
{ 
  d=nrow(rr)
  if(d<=2) return(rr)
  pp=matrix(0,d,d); pcmat=matrix(0,d,d); diag(pcmat)=1
  for(j in 2:d) pp[1,j]=rr[A[1,j],j]
  # tree 2
  for(j in 3:d) 
  { a1=A[1,j]; a2=A[2,j] 
  pp[2,j]=(rr[j,a2]-rr[j,a1]*rr[a1,a2])/sqrt((1-rr[j,a1]^2)*(1-rr[a1,a2]^2))
  }
  # remaining trees
  if(d>3)
  { for(ell in 3:(d-1))
  { for(j in (ell+1):d)
  { given=A[1:(ell-1),j]
  pp[ell,j]=partcor(rr,given,A[ell,j],j)  # assuming A[j,j]=j
  }
  }
  }
  for(ell in 1:(d-1))
  { for(j in (ell+1):d) 
  { alj=A[ell,j]; pcmat[alj,j]=pp[ell,j]; pcmat[j,alj]=pcmat[alj,j] }
  }
  list(pctree=pp, pcmat=pcmat)
}

