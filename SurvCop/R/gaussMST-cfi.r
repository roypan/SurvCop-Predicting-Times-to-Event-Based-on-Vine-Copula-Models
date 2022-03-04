#' Sequential MST based on partial correlations 
#'
#' @description 
#' Sequential MST based on partial correlations for best truncated Gaussian vine
#' This depends on the R package igraph
#' for the minimum spanning tree algorithm and operations on graphs.
#'
#' @param rmat correlation matrix
#' @param n sample size
#' @param CFIbd lower bound of CFI (comparative fit index) to stop, default 0.95, set as 1.01 for no truncation
#' @param iprint print flag for intermediate results on truncated vine
#' @param iprint2 print flag for intermediate results on graph objects 
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{VineA} -  vine array of dimension dxd
#'  \item \code{pcmat} -  matrix of partial correlations by tree
#'  \item \code{treeweight} -  vector of length d-1 with
#'                 sum_edge log(1-rho[edge]^2) for trees 1,...d-1
#'  \item \code{edgemat} -  matrix with columns node1 node2 level vector-of-conditioning
#'  \item \code{fitval} - cumsum(treeweight)/\eqn{sum_{1:(d-1)}} treeweight
#'  \item \code{CFIv} - vector of CFI values
#'  \item \code{ntrunc} - truncation level to reach CFIbd
#' }
#'
#' @details
#' CFI = 1- numerator/denominator, \cr
#' numerator = \eqn{\max(0,D_t-\nu_t)},\cr
#' denominator = \eqn{\max(0,D_t-\nu_t,D_0-\nu_0)},\cr
#' \eqn{D_0=-n\log(\det(R))},\cr
#' \eqn{D_t = n[-L_t(V) - \log(\det(R))]} is decreasing as t increases,\cr
#' \eqn{L_t(V) = -\sum_{1:t} \sum_e -\log(1-r_e^2)} (sum pcor up to tree t),\cr
#' \eqn{\nu_t = (d-t)(d-t-1)/2; \nu_0=d(d-1)/2},\cr
#' \eqn{L_t(V)} is increasing in t.
#'
#' @references
#' For use of CFI in the context of partial correlation vines, please see: \cr
#' Brechmann E C and Joe H (2015), Truncation of vine copulas using fit indices.
#' Journal of Multivariate Analysis, 138, 19-33.
#'
#' @examples
#' n=400
#' rmat=matrix(c(
#'   1.000,  0.172, -0.062,  0.385,  0.499, 0.267,  0.578,
#'   0.172,  1.000, -0.047,  0.493,  0.602, 0.396,  0.600,
#'  -0.062, -0.062,  1.000, -0.120, -0.092, 0.070, -0.070,
#'   0.385,  0.493, -0.120,  1.000,  0.607, 0.557,  0.742,
#'   0.499,  0.602, -0.092,  0.607,  1.000, 0.672,  0.883,
#'   0.267,  0.396,  0.070,  0.557,  0.672, 1.000,  0.685,
#'   0.578,  0.600, -0.070,  0.742,  0.883, 0.685,  1.000),7,7)
#' vinestr095 = gaussvine.mst(rmat, n, CFIbd=0.95, iprint=TRUE)
#' print(vinestr095$VineA)
#' print(vinestr095$CFIv)
#' vinestr099 = gaussvine.mst(rmat, n, CFIbd=0.99, iprint=TRUE)
#' print(vinestr099$VineA)
#' print(vinestr099$CFIv)
#' vinestr100 = gaussvine.mst(rmat, n, CFIbd=1.01, iprint=FALSE)
#' print(vinestr100$VineA)
#' print(vinestr100$CFIv)
#' 
#' @export
#'
gaussvine.mst = function(rmat, n, CFIbd=0.95, iprint=FALSE, iprint2=FALSE)
{
  rmat = as.matrix(rmat)
  d = dim(rmat)[1]
  colnames(rmat) = rownames(rmat) = paste("V",1:d,sep="")
  treeweight = rep(NA,d-1)

  # set up objects for return
  dd=(d*(d-1))/2
  node1=rep(0,dd); node2=rep(0,dd); lev=rep(0,dd);
  condmat=matrix(0,dd,d-2)
  icomv=rep(0,d); idv=rep(0,d); 
  vord=rep(0,d-1)
  treeweight=rep(0,d-1)
  pcv=rep(NA,dd); 

  CFIv=rep(0,d-1)
  logdet = log(det(rmat))
  nu0=dd
  D0=-n*logdet
  sumtreewt=0
  
  # tree 1
  # initialize the first graph
  # mode=upper to match C code
  g = graph.adjacency(rmat, mode="upper",weighted=TRUE,diag=FALSE)
  E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep = ",")
  for (ii in 1:ecount(g))
  { nodes=get.edges(g,ii)
    node1[ii]=nodes[1]; node2[ii]=nodes[2]
  }
  if(iprint2) print(g)
  # R indexes from 1 ; C indexes from 0
  vord[1]=1
  for(i in 2:(d-1)) { vord[i]=vord[i-1]+(d-i+1); }
  mst = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
  treeweight[1] = sum(log(1-E(mst)$weight^2))
  sumtreewt = sumtreewt+treeweight[1]
  Dt=n*(sumtreewt-logdet); nut=((d-1)*(d-2))/2
  numer=max(0,Dt-nut); denom=max(numer,D0-nu0)
  itree=1
  CFIv[itree]=1-numer/denom
  if(iprint2) cat("\n***tree1\n")
  if(iprint2) print(mst)
  # add info to objects: node1 node2 lev pcv
  for(i in 1:(d-1))
  { nodes=get.edges(mst,i)
    v1=nodes[1]; v2=nodes[2]
    eid=vord[v1]+(v2-v1-1)
    node1[eid]=v1; node2[eid]=v2;
    lev[eid]=1; pcv[eid]=rmat[v1,v2]
  }
  if(iprint2) { tmat=cbind(node1,node2,lev,pcv,condmat); print(tmat) }
  
  # loop through remaining trees, later check if itree=2 can be included 
  if(CFIv[1]<CFIbd)
  { for(itree in 2:(d-1))
    { nres = ecount(mst)
      if(iprint2) cat("\n***setting up tree", itree,"\n")
      if (iprint2) cat("\n***ecount=", nres,"\n")
      g = graph.full(nres) # all possible edges
      V(g)$name = E(mst)$name 
      if(itree>2) V(g)$iname = E(mst)$iname 
      ii=0
      for (i in 1:ecount(g))
      # g is full graph so it goes thru all pairs and check proximity
      { if(itree==2)
        { id = ends(g,i,names=F) # internal indices for edge
          temp = get.edges(mst,id) # edges of prev mst tree
          avec=temp[1,]; bvec=temp[2,]  # should be increasing
        }
        else
        { jj = ends(g,i,names=F) # internal indices for edge, 2-vector
          # different from itree=2 in next 4 lines
          id1=idv[jj[1]]; id2=idv[jj[2]]
          avec=c(node1[id1],node2[id1],condmat[id1,1:(itree-2)])
          bvec=c(node1[id2],node2[id2],condmat[id2,1:(itree-2)])
          avec=sort(avec); bvec=sort(bvec)
        }
        chk=proxcheck(avec,bvec)
        if(iprint2) cat(chk$iok,", ", avec,":",bvec,"\n")
        ok = chk$iok
        if (ok)
        { icomv=chk$icomv[1:(itree-1)]
          icond1=chk$icond1; icond2=chk$icond2;
          # new edge
          if(itree==2)
          { # new edge
            E(g)[i]$name = paste(paste(c(icond1,icond2),collapse = ","),icomv,sep="|")
            tem=(rmat[icond1,icond2]-rmat[icomv,icond1]*rmat[icomv,icond2])/
                 sqrt((1.-rmat[icomv,icond1]*rmat[icomv,icond1])*
                      (1.-rmat[icomv,icond2]*rmat[icomv,icond2])); 
          }
          else
          { E(g)[i]$name = paste(paste(c(icond1,icond2),collapse=","),paste(icomv,collapse = ","),sep="|")
            tem = partcor(rmat,icomv,icond1,icond2)
          }
          E(g)[i]$weight = tem
          ii=ii+1
          eid=vord[icond1]+(icond2-icond1-1)
          pcv[eid]=tem;
          lev[eid]=-itree; # later change to positive if in mst
          condmat[eid,1:(itree-1)]=icomv;
          if(iprint2) cat(ii,tem,eid,"\n")
          E(g)[i]$iname=eid  
        }
        E(g)[i]$todel = !ok
      }
      g = delete.edges(g, E(g)[E(g)$todel])
      if(iprint2) { print(g); }
      # setup for next mst
      mst = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))
      treeweight[itree] = sum(log(1-E(mst)$weight^2))
      sumtreewt = sumtreewt+treeweight[itree]
      Dt=n*(sumtreewt-logdet); nut=((d-itree)*(d-itree-1))/2
      numer=max(0,Dt-nut); denom=max(numer,D0-nu0)
      CFIv[itree]=1-numer/denom
      if(iprint) cat("CFI", itree, Dt,nut, numer, denom, 1-numer/denom,"\n")
      # set selected edges to have lev[eid]=itree
      for(i in 1:(d-itree))
      { eid=E(mst)[i]$iname
        v1=node1[eid]; v2=node2[eid]
        if(iprint2) cat(i,eid,v1,v2,"\n")
        lev[eid]=itree; 
        idv[i]=eid
      }
      if(iprint2) { tmat=cbind(node1,node2,lev,pcv,condmat); print(tmat) }
      if(CFIv[itree]>=CFIbd) break;
    }
  }

  # get truncated vine array and corresponding matrix of partial correlations
  edgemat=cbind(node1,node2,lev,condmat)
  if(itree>=(d-1)) out=edges2array(d,edgemat,pcv)
  else out=edges2array.trunc(d,itree,edgemat,pcv,iprint=iprint)
  #print(out$VineA)
  # treeweight vector
  #logdet = log(det(rmat))
  #cat(logdet,sum(treeweight),"\n") # same
  fitval=cumsum(treeweight)
  fitval=fitval/logdet
  #cat(fitval,"\n")
  if(iprint) cat(CFIv,"\n")
  list(edgemat=edgemat,pcvec=pcv,VineA=out$VineA,pcmat=out$pcmat,
    treeweight=treeweight,fitval=fitval, CFIv=CFIv,ntrunc=itree)
}

# This function is used by gaussvine.mst and yleaf.mst
#' Proximity condition check
#'
#' @description
#' Check of proximity condition for ordered vectors avec, bvec of length m
#'
#' @param avec (with m distinct values) sorted in increasing order
#' @param bvec (with m distinct values) sorted in increasing order
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{iok} - equals 1 if proximity condition is satisfied, 0 if not
#'  \item \code{icomv} - vector of common variables (length m-1 if iok=1)
#'  \item \code{icond1} - smaller index of two conditioned variables if iok=1
#'  \item \code{icond2} - larger index of two conditioned variables if iok=1
#' }
#'
#' @examples
#' proxcheck(c(1,4,5),c(2,4,6))
#' proxcheck(c(1,4,5),c(1,4,6))
#'
#' @export
#'
proxcheck=function(avec,bvec)
{ m=length(avec)
  mx=max(c(avec,bvec))
  avec=c(avec,mx+1)
  bvec=c(bvec,mx+2)
  icond1=0; icond2=0
  idifv=rep(0,2*m+1)
  i1=1; i2=1; ncom=0; ndif=0;
  icomv=rep(0,m);
  a=avec[i1]; b=bvec[i2];
  while(i1<=m & i2<=m)
  { 
    if(a==b) { ncom=ncom+1; icomv[ncom]=a; i1=i1+1; i2=i2+1; a=avec[i1]; b=bvec[i2]; }
    else if(a<b) { ndif=ndif+1; idifv[ndif]=a; i1=i1+1; a=avec[i1]; }
    else { ndif=ndif+1; idifv[ndif]=b; i2=i2+1; b=bvec[i2]; }
  }
  if(i2>m & i1<=m)
  { for(i in i1:m) { ndif=ndif+1; idifv[ndif]=avec[i]; } }
  if(i1>m & i2<=m)
  { for(i in i2:m) { ndif=ndif+1; idifv[ndif]=bvec[i]; } }
  if(ndif==2) { icond1=idifv[1]; icond2=idifv[2]; iok=1; }
  else { iok=0; }
  list(iok=iok,icomv=icomv,icond1=icond1,icond2=icond2)
}

# This function is used by gaussvine.mst().
# edgemat=cbind(node1,node2,lev,condmat)  : d*(d-1)/2 x (3+(d-2))
# pcv is corresponding vector of partial correlations 
#' Edge matrix to a vine array
#' 
#' @description
#' Convert a vine stored as an edge matrix (columns node1, node2, level, conditioning variables) 
#' to a vine array
#' 
#' @param d dimension or number of variables 
#' @param edgemat edge matrix of dimension choose(d,2) x (d+1)
#'      with columns node1, node2, vine tree level, given[1], ..., given[d-2]
#'    (it is not checked if this matches a proper vine array).
#'    Each row represents an edge of the vine.
#' @param pcv corresponding vector of partial correlations for edges of vine 
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{VineA} -  vine array of dimension dxd
#'  \item \code{pcmat} -  matrix of partial correlations by tree
#' }
#' 
#' @details 
#' WARNING: It has not been checked what happens with incorrect input.
#' 
#' @examples
#' edgemat6=matrix(c(
#'    1,2, 1, 0,0,0,0,
#'    1,3, 1, 0,0,0,0,
#'    2,4, 1, 0,0,0,0,
#'    1,5, 1, 0,0,0,0,
#'    2,6, 1, 0,0,0,0,
#'    2,3, 2, 1,0,0,0,
#'    1,4, 2, 2,0,0,0,
#'    3,5, 2, 1,0,0,0,
#'    4,6, 2, 2,0,0,0,
#'    3,4, 3, 2,1,0,0,
#'    2,5, 3, 1,3,0,0,
#'    1,6, 3, 2,4,0,0,
#'    4,5, 4, 1,3,2,0,
#'    3,6, 4, 2,4,1,0,
#'    5,6, 5, 2,4,1,3), 15,7, byrow=TRUE)
#' colnames(edgemat6)=c("node1","node2","level","giv1","giv2","giv3","giv4")
#' pcvec=c(.6,.7,.6,.7,.5, .4,.4,.3,.3, .2,.2,.1, .1,.1, .1) 
#' edges2array(d=6,edgemat6,pcvec)
#' 
#' @export
#'
edges2array=function(d,edgemat,pcv)
{ A=matrix(0,d,d); pc=matrix(NA,d,d)
  mat=edgemat # this matrix will have rows deleted as columns of A are filled
  for(j in d:2)
  { k=which(mat[,3]==(j-1))  # should only be one
    #cat("col=", j, "k=", k,"\n")
    #w=mat[k,1]; A[j,j]=w; A[j-1,j]=mat[k,2]
    w=mat[k,2]; A[j,j]=w; A[j-1,j]=mat[k,1] # larger number at last col
    pc[j-1,j]=pcv[k]
    if(j==2) { A[1,1]=A[1,2] } 
    else
    { for(ii in (j-2):1)  # trees j-2 to 1
      { itree=which(mat[,3]==ii)
        #cat(ii,"\n"); print(itree)
        for(jj in itree)
        { cs=mat[jj,]; pctem=pcv[jj]
          if(cs[1]==w) { A[ii,j]=cs[2]; break }
          else if(cs[2]==w) { A[ii,j]=cs[1]; break }
        }
        pc[ii,j]=pctem
        mat[jj,]=0
      }
    }
    mat[k,]=0
    #cat("column j\n"); print(A[,j])
  }
  return(list(VineA=A,pcmat=pc))
}

# This function is used by edges2array.trunc
# find leaf variables from tree with conditioned and conditioning variables
# return leaf variables in a vector leafv
findleaf=function(conditioned,conditioning,iprint=F)
{ given=unique(c(conditioning))
  uniq=sort(unique(c(conditioned)))
  tem=table(c(conditioned))
  i0=(tem==1)
  leafs=setdiff(uniq[i0],given)
  if(iprint) cat("findleaf: ", leafs,"\n")
  leafs  
}

# This function is used by gaussvine.mst()
# d = number of varibles
# mtree = truncation level 
# edgemat=cbind(node1,node2,lev,condmat)  : d*(d-1)/2 x (3+(d-2))
# pcv is corresponding vector of partial correlations 
#' Edge matrix of a truncated vine to a vine array 
#' 
#' @description
#' Convert a truncated vine stored as an edge matrix (columns node1, node2, level, conditioning variables) 
#' to a truncated vine array
#' 
#' @param d dimension or number of variables 
#' @param mtree truncation level 
#' @param edgemat edge matrix of dimension [(d-1)+...+(d-mtree)] x (mtree+3)
#'      with columns node1, node2, vine tree level, given[1], ..., given[mtree]
#'    (it is not checked if this matches a proper vine array).
#'    Each row represents an edge of the vine.
#' @param pcv corresponding vector of partial correlations for edges of truncated vine 
#' @param iprint print flag for intermediate results, default FALSE
#'
#' @return Returns a list with the following named components:
#' \enumerate{
#'  \item \code{VineA} -  truncated vine array of dimension dxd, with 0s for trees at level greater than mtree
#'  \item \code{pcmat} -  matrix of partial correlations by tree
#' }
#' 
#' @details 
#' WARNING: It has not been checked what happens with incorrect input.
#' 
#' @examples
#' edgemat6=matrix(c(
#'    1,2, 1, 0,0,0,0,
#'    1,3, 1, 0,0,0,0,
#'    2,4, 1, 0,0,0,0,
#'    1,5, 1, 0,0,0,0,
#'    2,6, 1, 0,0,0,0,
#'    2,3, 2, 1,0,0,0,
#'    1,4, 2, 2,0,0,0,
#'    3,5, 2, 1,0,0,0,
#'    4,6, 2, 2,0,0,0,
#'    3,4, 3, 2,1,0,0,
#'    2,5, 3, 1,3,0,0,
#'    1,6, 3, 2,4,0,0,
#'    4,5, 4, 1,3,2,0,
#'    3,6, 4, 2,4,1,0,
#'    5,6, 5, 2,4,1,3), 15,7, byrow=TRUE)
#' colnames(edgemat6)=c("node1","node2","level","giv1","giv2","giv3","giv4")
#' edgemat.trunc2=edgemat6[1:9,1:5]
#' pcvec2=c(.6,.7,.6,.7,.5, .4,.4,.3,.3) 
#' edges2array.trunc(d=6,mtree=2,edgemat.trunc2,pcvec2)
#' edgemat.trunc3=edgemat6[1:12,1:6]
#' pcvec3=c(.6,.7,.6,.7,.5, .4,.4,.3,.3, .2,.2,.1) 
#' edges2array.trunc(d=6,mtree=3,edgemat.trunc3,pcvec3)
#' iperm=c(2,4,12,10,9,3,8,5,11,1,6,7)
#' edges2array.trunc(d=6,mtree=3,edgemat.trunc3[iperm,],pcvec3[iperm])
#' # should be the same, because row permutation does not affect vine
#'
#' @export
#'
edges2array.trunc=function(d,mtree,edgemat,pcv,iprint=FALSE)
{ A=matrix(0,d,d); pc=matrix(NA,d,d)
  mat=edgemat # this matrix will have rows deleted as columns of A are filled
  nc=ncol(mat)
  # complete columns d,...,m+2 of vine array (maybe through iterations)
  ntot=0; lastcol=d
  while(ntot<d-mtree-1) 
  { ind=which(mat[,3]==mtree)
    conditioned=mat[ind,1:2]
    conditioning=mat[ind,4:nc]
    leafs=findleaf(conditioned,conditioning,iprint)
    # dimension is between 1 and d-mtree-1
    # reverse order so that leaf with largest index goes to last column
    leafs=rev(leafs)
    if(iprint) cat("leafs=", leafs,"\n")
    nleaf=length(leafs)
    ntot=ntot+nleaf
    for(i in 1:nleaf)
    { icol=lastcol-i+1
      if(icol==(mtree+1)) { break }
      w=leafs[i]
      A[icol,icol]=w
      for(ii in (mtree:1))  # trees m to 1
      { itree=which(mat[,3]==ii)
        #if(iprint) { cat(ii,"\n"); print(itree) }
        for(jj in itree)
        { cs=mat[jj,]
          if(cs[1]==w) { A[ii,icol]=cs[2]; break }
          else if(cs[2]==w) { A[ii,icol]=cs[1]; break }
        }
        pc[ii,icol]=pcv[jj]
        mat[jj,]=0
      }
    }
    lastcol=icol-1
  }
  # columns mtree+1 to 2 as before
  for(j in (mtree+1):2)
  { # below Ok for columns 1,2,...,mtree
    k=which(mat[,3]==(j-1))  # should only be one
    #cat("col=", j, "k=", k,"\n")
    w=mat[k,2]; A[j,j]=w; A[j-1,j]=mat[k,1]
    pc[j-1,j]=pcv[k]
    if(j==2) { A[1,1]=A[1,2] } 
    else
    { for(ii in (j-2):1)  # trees j-2 to 1
      { itree=which(mat[,3]==ii)
        for(jj in itree)
        { cs=mat[jj,]
          if(cs[1]==w) { A[ii,j]=cs[2]; break }
          else if(cs[2]==w) { A[ii,j]=cs[1]; break }
        }
        pc[ii,j]=pcv[jj]
        mat[jj,]=0
      }
    }
    mat[k,]=0
  }
  return(list(VineA=A,pcmat=pc))
}
