! Frank copula cdf
! f90 standalone with main program
! gfortran -o franksub franksub.f90
! franksub 
! link to R, use version with main program commented out
! gfortran -fpic -c franksub.f90
! gfortran -shared -o mylib.so franksub.o
!  later add other .f90 code to the .so library

! comment out this to link to R
! sample main program
!program frankmain
!  implicit none
!  integer n,i
!  double precision, dimension(:), allocatable :: uvec,vvec,cparvec,cdfvec
!
!  !read *,n  ! to read from an input file
!  n=10   ! fix a value in main program
!  allocate ( uvec(n),vvec(n),cparvec(n),cdfvec(n) )  ! dynamic allocation
!  ! set values in main program
!  do i=1,n
!    uvec(i)= i/(n+1.d0)
!    vvec(i)= 0.5d0*i/n
!    cparvec(i)= i*0.4d0
!  end do
!  call pfrank(n,uvec,vvec,cparvec,cdfvec)
!  do i=1,n
!    print "(i3,4f10.6)", i,uvec(i),vvec(i),cparvec(i),cdfvec(i)
!  end do
!  deallocate (uvec,vvec,cparvec,cdfvec)
!  stop
!  end


! Frank copula cdf via subroutine for linking to R 
! n = length of vectors u and v and cpar
! uvec,vvec vectors of values in (0,1)
! cparvec = vector of copula parameters of length n 
subroutine pfrank(n,uvec,vvec,cparvec,cdf)
  implicit none 
  integer n,i
  double precision uvec(n),vvec(n),cparvec(n),cdf(n)
  double precision u,v,cpar,cpar1,tem1,tem2,tem
  do i=1,n
    cpar=cparvec(i); u=uvec(i); v=vvec(i)
    if (cpar==0.d0) then 
      cdf(i) = u*v
    else
      cpar1=1.d0-exp(-cpar)
      tem1=exp(-cpar*u); tem2=exp(-cpar*v)
      tem=cpar1-(1.d0-tem1)*(1.d0-tem2)
      cdf(i)=(-log(tem/cpar1)/cpar)
    end if
  end do
  return 
  end

subroutine dfrank(n,uvec,vvec,cparvec,pdf)
  implicit none
  integer n,i
  double precision uvec(n),vvec(n),cparvec(n),pdf(n)
  double precision u,v,cpar,cpar1,tem1,tem2,tem,tem3
  do i=1,n
    cpar=cparvec(i); u=uvec(i); v=vvec(i)
    cpar1=1.d0-exp(-cpar)
    tem1=exp(-cpar*u); tem2=exp(-cpar*v)
    tem3 = cpar1 * tem1 * tem2 * cpar
    tem=cpar1-(1.d0-tem1)*(1.d0-tem2)
    pdf(i) = tem3/(tem * tem)
  end do
  return
  end

subroutine pcondfrank(n,vvec,uvec,cparvec,ccdf)
  implicit none
  integer n,i
  double precision uvec(n),vvec(n),cparvec(n),ccdf(n)
  double precision u,v,cpar,cpar1,tem
  do i=1,n
    cpar=cparvec(i); u=uvec(i); v=vvec(i)
    if (cpar==0.d0) then
      cpar = 1.d-10
    end if
    cpar1=1.d0-exp(-cpar)
    tem = 1.d0-exp(-cpar * u)
    ccdf(i) = (1.d0 - tem)/(cpar1/(1.d0 - exp(-cpar * v)) - tem)
  end do
  return
  end
