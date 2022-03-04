BivCopCDF=function(u,v,fam,par1,par2=0)
{  #if(!is.loaded("pfrank"))  dyn.load("./mylib.so") # inside an R package this line can be commented
  out= .Fortran("bivcopcdf",
      as.double(u), as.double(v), as.integer(fam),
      as.double(par1),as.double(par2),
      cdf=as.double(0))
  out$cdf
}

BivCopPDF=function(u,v,fam,par1,par2=0)
{ #if(!is.loaded("pfrank"))  dyn.load("./mylib.so") # inside an R package this line can be commented
  out= .Fortran("bivcoppdf",
      as.double(u), as.double(v), as.integer(fam),
      as.double(par1),as.double(par2),
      pdf=as.double(0))
  out$pdf
}

BivCopCCDF=function(v,u,fam,par1,par2=0)
{ #if(!is.loaded("pfrank"))  dyn.load("./mylib.so") # inside an R package this line can be commented
  out= .Fortran("bivcopccdf",
      as.double(v), as.double(u), as.integer(fam),
      as.double(par1),as.double(par2),
      ccdf=as.double(0))
  out$ccdf
}
