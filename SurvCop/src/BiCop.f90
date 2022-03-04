! sample f90 main program calling functions that use C functions
!program BiCop
!  implicit none
!  double precision u,v,cpar,pdf,cdf
!  double precision BiCopCDF
!  double precision BiCopPDF
!  
!  read *, u,v,cpar
!  do while(u>0.d0 .and. u<1.d0) 
!    cdf=BiCopCDF(u,v,5,cpar,0.d0)
!    pdf=BiCopPDF(u,v,5,cpar,0.d0)
!    print '(5f9.6)', u,v,cpar,cdf,pdf
!    read *, u,v,cpar
!    end do
!  stop
!  end

subroutine bivcopcdf(u,v,fam,par1,par2,cdf)
  implicit none
  integer fam,infin(2), nu
  double precision u,v,par1,par2,cdf,lower(2),upper(2)
  double precision cpar,cpar1,tem1,tem2,tem,tem3,del,thl,sm
  double precision qnorms, mvbvn, qt, mvbvt
  select case (fam)
  case(1) ! gaussian
    cpar=par1
    tem1 = qnorms(u)
    tem2 = qnorms(v)
    upper(1) = tem1; upper(2) = tem2
    infin(1) = 0; infin(2) = 0 
    cdf = mvbvn(lower, upper, infin, cpar)
  case(2) ! t
    cpar=par1; cpar1=par2
    if (cpar1 < 1.d0) then
      cpar1 = 1.d0
    end if
    nu = nint(cpar1); cpar1 = nint(cpar1)
    tem1 = qt(u, cpar1)
    tem2 = qt(v, cpar1)
    upper(1) = tem1; upper(2) = tem2
    infin(1) = 0; infin(2) = 0 
    cdf = mvbvt(nu, lower, upper, infin, cpar)
  case(3) ! clayton
    cpar=par1
    tem1 = u**(-cpar)
    tem2 = v**(-cpar)
    cdf = (tem1 + tem2 - 1.d0)**(-1.d0/cpar)
  case(13) ! reflected clayton
    cpar=par1
    tem1 = (1.d0 - u)**(-cpar)
    tem2 = (1.d0 - v)**(-cpar)
    cdf = u + v - 1.d0 + (tem1 + tem2 - 1.d0)**(-1.d0/cpar)
  case(4) ! gumbel
    cpar=par1
    tem1 = ((-log(u))**cpar)
    tem2 = ((-log(v))**cpar)
    tem = (tem1 + tem2)**(1.d0/cpar)
    cdf = exp(-tem)
  case(14) ! reflected gumbel
    cpar=par1
    tem1 = ((-log(1.d0-u))**cpar)
    tem2 = ((-log(1.d0-v))**cpar)
    tem = (tem1 + tem2)**(1.d0/cpar)
    cdf = exp(-tem) + u + v - 1.d0
  case(5) ! frank
    cpar=par1
    if (cpar==0.d0) then 
      cdf = u*v
    else
      cpar1=1.d0-exp(-cpar)
      tem1=exp(-cpar*u); tem2=exp(-cpar*v)
      tem=cpar1-(1.d0-tem1)*(1.d0-tem2)
      cdf=(-log(tem/cpar1)/cpar)
    end if
  case(6) ! joe
    cpar=par1
    tem1 = (1.d0 - u)**cpar
    tem2 = (1.d0 - v)**cpar
    tem3 = tem1 + tem2 - tem1 * tem2
    tem = tem3**(1.d0/cpar)
    cdf = 1.d0 - tem
  case(16) ! reflected joe
    cpar=par1
    tem1 = u**cpar
    tem2 = v**cpar
    tem3 = tem1 + tem2 - tem1 * tem2
    tem = tem3**(1.d0/cpar)
    cdf = u + v - tem
  case(7) ! BB1
    cpar=par1; cpar1=par2
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = (u**(-cpar) - 1.d0)
    tem2 = (v**(-cpar) - 1.d0)
    tem1 = tem1**cpar1
    tem2 = tem2**cpar1
    sm = tem1 + tem2
    tem = sm**(del)
    cdf = (1.d0 + tem)**(-thl)
  case(17) ! reflected BB1
    cpar=par1; cpar1=par2
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = ((1.d0-u)**(-cpar) - 1.d0)
    tem2 = ((1.d0-v)**(-cpar) - 1.d0)
    tem1 = tem1**cpar1
    tem2 = tem2**cpar1
    sm = tem1 + tem2
    tem = sm**(del)
    cdf = u + v - 1.d0 + (1.d0 + tem)**(-thl)
  case(10) ! BB8
    cpar=par1; cpar1=par2
    tem1 = 1.d0 - (1.d0 - cpar1 * u)**cpar
    tem2 = 1.d0 - (1.d0 - cpar1 * v)**cpar
    tem = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    sm = (1.d0 - tem * tem1 * tem2)**(1.d0/cpar)
    cdf = (1.d0 - sm)/cpar1
  case(20) ! reflected BB8
    cpar=par1; cpar1=par2
    tem1 = 1.d0 - (1.d0 - cpar1 * (1.d0-u))**cpar
    tem2 = 1.d0 - (1.d0 - cpar1 * (1.d0-v))**cpar
    tem = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    sm = (1.d0 - tem * tem1 * tem2)**(1.d0/cpar)
    cdf = u + v - 1.d0 + (1.d0 - sm)/cpar1
  case default ! not implemented 
    cdf = -1.d0
  end select
  return
end

subroutine bivcoppdf(u,v,fam,par1,par2,pdf)
  implicit none
  integer fam
  double precision u,v,par1,par2,pdf
  double precision cpar,cpar1,tem1,tem2,tem,tem3,l1,l2,sm,pi,del,thl,r11,r12
  double precision qnorms,dnorms,lgamma,qt
  select case (fam)
  case(1) ! gaussian
    cpar=par1
    tem1 = qnorms(u)
    tem2 = qnorms(v)
    tem = (tem1**2.d0 + tem2**2.d0 - 2.d0 * cpar * tem1 * tem2) / (1.d0 - cpar ** 2.d0)
    pi = 4.d0*datan(1.d0)
    tem3 = sqrt(1.d0 - cpar**2.d0) * (2.d0 * pi)
    tem3 = exp(-tem / 2.d0)/tem3
    pdf = tem3/(dnorms(tem1) * dnorms(tem2))
  case(2) ! t
    cpar=par1; cpar1=par2
    pi = 4.d0*datan(1.d0)
    tem = exp(lgamma((cpar1 + 2.d0)/2.d0) - lgamma(cpar1/2.d0))/(pi * cpar1)
    tem = tem / sqrt(1.d0 - cpar * cpar)
    tem3 = -(cpar1 + 2.d0)/2.d0
    tem1 = qt(u, cpar1)
    tem2 = qt(v, cpar1)
    call dt(tem1, cpar1, l1)
    call dt(tem2, cpar1, l2)
    r11 = 1.d0/(1.d0 - cpar * cpar)
    r12 = -cpar * r11
    sm = 1.d0 + (tem1 * tem1 * r11 + tem2 * tem2 * r11 + 2.d0 * tem1 * tem2 * r12) / cpar1
    pdf = tem * (sm**tem3)/(l1 * l2)
  case(3) ! clayton
    cpar=par1
    tem1 = u**(-cpar)
    tem2 = v**(-cpar)
    pdf = (tem1 + tem2 - 1.d0)**(-1.d0/cpar - 2.d0) * (1.d0 + cpar) * tem1 * tem2/(u * v)
  case(13) ! reflected clayton
    cpar=par1; u=1.d0-u; v=1.d0-v
    tem1 = u**(-cpar)
    tem2 = v**(-cpar)
    pdf = (tem1 + tem2 - 1.d0)**(-1.d0/cpar - 2.d0) * (1.d0 + cpar) * tem1 * tem2/(u * v)
  case(4) ! gumbel
    cpar=par1
    l1 = -log(u); l2 = -log(v)
    tem1 = l1**cpar
    tem2 = l2**cpar
    sm = tem1 + tem2
    tem = sm**(1.d0/cpar)
    tem3 = exp(-tem)
    tem3 = tem3 * tem * tem1 * tem2 * (tem + cpar - 1.d0)
    pdf = tem3 / (sm * sm * l1 * l2 * u * v)
  case(14) ! reflected gumbel
    cpar=par1; u=1-u; v=1-v
    l1 = -log(u); l2 = -log(v)
    tem1 = l1**cpar
    tem2 = l2**cpar
    sm = tem1 + tem2
    tem = sm**(1.d0/cpar)
    tem3 = exp(-tem)
    tem3 = tem3 * tem * tem1 * tem2 * (tem + cpar - 1.d0)
    pdf = tem3 / (sm * sm * l1 * l2 * u * v)
  case(5) ! frank
    cpar=par1
    cpar1=1.d0-exp(-cpar)
    tem1=exp(-cpar*u); tem2=exp(-cpar*v)
    tem3 = cpar1 * tem1 * tem2 * cpar
    tem=cpar1-(1.d0-tem1)*(1.d0-tem2)
    pdf = tem3/(tem * tem)
  case(6) ! joe
    cpar=par1
    l1 = 1.d0 - u
    l2 = 1.d0 - v
    tem1 = l1**cpar
    tem2 = l2**cpar
    sm = tem1 + tem2 - tem1 * tem2
    tem = sm**(1.d0/cpar)
    tem = tem * ((cpar - 1.d0) * tem1 * tem2 + tem1 * tem1 * tem2 + tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
    tem = tem/(sm * sm)
    pdf = tem/(l1 * l2)
  case(16) ! reflected joe
    cpar=par1
    l1 = u
    l2 = v
    tem1 = l1**cpar
    tem2 = l2**cpar
    sm = tem1 + tem2 - tem1 * tem2
    tem = sm**(1.d0/cpar)
    tem = tem * ((cpar - 1.d0) * tem1 * tem2 + tem1 * tem1 * tem2 + tem1 * tem2 * tem2 - tem1 * tem1 * tem2 * tem2)
    tem = tem/(sm * sm)
    pdf = tem/(l1 * l2)
  case(7) ! BB1
    cpar=par1; cpar1=par2
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = (u**(-cpar) - 1.d0)
    tem2 = (v**(-cpar) - 1.d0)
    l1 = tem1**cpar1
    l2 = tem2**cpar1
    sm = l1 + l2
    tem3 = sm**(del)
    tem = (1.d0 + tem3)**(-thl - 2.d0) * (cpar * (cpar1 - 1.d0) + (cpar * cpar1 + 1.d0) * tem3)
    pdf = tem * tem3 * l1 * l2 * (tem1 + 1.d0) * (tem2 + 1.d0)/sm/sm/tem1/tem2/u/v
  case(17) ! reflected BB1
    cpar=par1; cpar1=par2
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = ((1.d0-u)**(-cpar) - 1.d0)
    tem2 = ((1.d0-v)**(-cpar) - 1.d0)
    l1 = tem1**cpar1
    l2 = tem2**cpar1
    sm = l1 + l2
    tem3 = sm**(del)
    tem = (1.d0 + tem3)**(-thl - 2.d0) * (cpar * (cpar1 - 1.d0) + (cpar * cpar1 + 1.d0) * tem3)
    pdf = tem * tem3 * l1 * l2 * (tem1 + 1.d0) * (tem2 + 1.d0)/sm/sm/tem1/tem2/(1.d0-u)/(1.d0-v)
  case(10) ! BB8
    cpar=par1; cpar1=par2
    tem1 = (1.d0 - cpar1 * u)**cpar
    tem2 = (1.d0 - cpar1 * v)**cpar
    l1 = 1.d0 - tem1
    l2 = 1.d0 - tem2
    sm = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    tem = (1.d0 - sm * l1 * l2)**(1.d0/cpar - 2.d0)
    pdf = sm * cpar1 * tem * (cpar - sm * l1 * l2) * tem1 * tem2/(1.d0 - cpar1 * u)/(1.d0 - cpar1 * v)
  case(20) ! reflected BB8
    cpar=par1; cpar1=par2
    tem1 = (1.d0 - cpar1 * (1.d0-u))**cpar
    tem2 = (1.d0 - cpar1 * (1.d0-v))**cpar
    l1 = 1.d0 - tem1
    l2 = 1.d0 - tem2
    sm = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    tem = (1.d0 - sm * l1 * l2)**(1.d0/cpar - 2.d0)
    pdf = sm * cpar1 * tem * (cpar - sm * l1 * l2) * tem1 * tem2/(1.d0 - cpar1 * (1.d0-u))/(1.d0 - cpar1 * (1.d0-v))
  case default ! not implemented 
    pdf = 0.d0
  end select 
  return
end

subroutine bivcopccdf(v,u,fam,par1,par2,ccdf)
  implicit none
  integer fam, infin(2), nu
  double precision u,v,u1,v1,par1,par2,ccdf,lower(2),upper(2)
  double precision cpar,cpar1,tem,tem1,tem2,tem3,l1,l2,eta,sm,del,thl
  double precision qnorms,pnorms,qt,mvbvt,pt
  select case (fam)
  case(1) ! gaussian
    cpar=par1
    if (v <= 0.d0 .or. u <= 0.d0 .or. u >= 1.d0) then 
      ccdf = 0.d0
    else if (v == 1.d0) then
      ccdf = 1.d0
    else
      ccdf = pnorms((qnorms(v) - cpar * qnorms(u)) / sqrt(1.d0 - cpar ** 2.d0))
    end if
  case(2) ! t
    cpar=par1; cpar1=par2
    tem1 = qt(u, cpar1)
    tem2 = qt(v, cpar1)
    tem3 = cpar * tem1
    tem = (1.d0 - cpar * cpar) * (cpar1 + tem1 * tem1)/(cpar1 + 1.d0)
    sm = (tem2 - tem3)/sqrt(tem)
    ccdf = pt(sm, cpar1 + 1.d0)
  case(3) ! clayton
    cpar=par1
    tem = v**(-cpar) - 1.d0
    tem = tem * (u**cpar) + 1.d0
    ccdf = tem**(-1.d0 - 1.d0/cpar)
  case(13) ! reflected clayton
    cpar=par1; u1=1.d0 - u; v1=1.d0 - v
    tem = v1**(-cpar) - 1.d0
    tem = tem * (u1**cpar) + 1.d0
    ccdf = 1.d0 - tem**(-1.d0 - 1.d0/cpar)
  case(4) ! gumbel
    cpar=par1
    l1 = -log(u); l2 = -log(v)
    tem1 = l1**cpar
    tem2 = l2**cpar
    sm = tem1 + tem2
    tem = sm**(1.d0/cpar)
    tem3 = exp(-tem)
    tem3 = tem3 * (1.d0 + tem2/tem1)**(-1.d0 + 1.d0/cpar)
    ccdf = tem3 / u
  case(14) ! reflected gumbel
    cpar=par1
    l1 = -log(1.d0-u); l2 = -log(1.d0-v)
    tem1 = l1**cpar
    tem2 = l2**cpar
    tem = (tem1 + tem2)**(1.d0/cpar)
    tem3 = (1.d0 + tem2/tem1)**(-1.d0 + 1.d0/cpar)
    ccdf = 1.d0 - exp(-tem) * tem3/(1.d0 - u)
  case(5) ! frank
    cpar=par1
    cpar1 = 1.d0 - exp(-cpar)
    tem = 1.d0 - exp(-cpar * u)
    ccdf = (1.d0 - tem)/(cpar1/(1.d0 - exp(-cpar * v)) - tem)
  case(6) ! joe
    cpar=par1
    tem2 = (1.d0  - v)**cpar
    tem1 = (1.d0  - u)**cpar
    tem = 1.d0 + tem2/tem1 - tem2
    tem = tem**(-1.d0 + 1.d0/cpar)
    ccdf = tem * (1.d0 - tem2)
  case(16) ! reflected joe
    cpar=par1; u1=1.d0 - u; v1=1.d0 - v
    tem2 = (1.d0  - v1)**cpar
    tem1 = (1.d0  - u1)**cpar
    tem = 1.d0 + tem2/tem1 - tem2
    tem = tem**(-1.d0 + 1.d0/cpar)
    ccdf = 1.d0 - tem * (1.d0 - tem2)
  case(7) ! BB1
    cpar=par1; cpar1=par2
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = (u**(-cpar) - 1.d0)
    tem2 = (v**(-cpar) - 1.d0)
    l1 = tem1**cpar1
    l2 = tem2**cpar1
    sm = l1 + l2
    tem3 = sm**(del)
    tem = (1.d0 + tem3)**(-thl - 1.d0)
    ccdf = tem * tem3 * l1 * (tem1 + 1.d0)/sm/tem1/u
  case(17) ! reflected BB1
    cpar=par1; cpar1=par2; u1=1.d0 - u; v1=1.d0 - v
    del = 1.d0/cpar1
    thl = 1.d0/cpar
    tem1 = (u1**(-cpar) - 1.d0)
    tem2 = (v1**(-cpar) - 1.d0)
    l1 = tem1**cpar1
    l2 = tem2**cpar1
    sm = l1 + l2
    tem3 = sm**(del)
    tem = (1.d0 + tem3)**(-thl - 1.d0)
    ccdf = 1.d0 - tem * tem3 * l1 * (tem1 + 1.d0)/sm/tem1/u1
  case(10) ! BB8
    cpar=par1; cpar1=par2
    sm = (1.d0 - cpar1 * u)**cpar
    l1 = 1.d0 - sm
    l2 = 1.d0 - (1.d0 - cpar1 * v)**cpar
    tem1 = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    tem = (1.d0 - tem1 * l1 * l2)**(1.d0/cpar - 1.d0)
    tem2 = (1.d0 - cpar1 * u)/sm
    ccdf = tem1 * l2 * tem/tem2
  case(20) ! reflected BB8
    cpar=par1; cpar1=par2; u1=1.d0 - u; v1=1.d0 - v
    sm = (1.d0 - cpar1 * u1)**cpar
    l1 = 1.d0 - sm
    l2 = 1.d0 - (1.d0 - cpar1 * v1)**cpar
    tem1 = 1.d0/(1.d0 - (1.d0 - cpar1)**cpar)
    tem = (1.d0 - tem1 * l1 * l2)**(1.d0/cpar - 1.d0)
    tem2 = (1.d0 - cpar1 * u1)/sm
    ccdf = 1.d0 -  tem1 * l2 * tem/tem2
  case default ! not implemented 
    ccdf = -1.d0
  end select
  return
end
