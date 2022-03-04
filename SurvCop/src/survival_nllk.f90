!program survivalnllkmain
!  implicit none
!  integer n,i,d
!  integer A(5,5), fam(5,5), temi(25)
!  double precision temi2(20), nllk, parammat1(5,5), parammat2(5,5)
!  double precision datacdfplus(4,5), datacdfminus(4,5)
!  logical censorstatus(4), isdisc(5)
!
!  n=4
!  d=5
!  temi = (/1,0,0,0,0, 1,3,0,0,0, 3,1,2,0,0, 3,2,1,4,0, 4,3,2,1,5/)
!  A = reshape(temi, (/5,5/))
!  temi = (/0,0,0,0,0, 5,0,0,0,0, 5,5,0,0,0, 5,5,5,0,0, 5,5,5,5,0/)
!  fam = reshape(temi, (/5,5/))
!  temi2 = (/.5,.4,.7,.5, .4,.6,.9,.7, .5,.3,.2,.4, .4,.4,.7,.4, .7,.4,.5,.4/)
!  datacdfplus = reshape(temi2, (/4,5/))
!  datacdfminus = reshape(temi2, (/4,5/))
!  parammat1 = 4.d0
!  parammat2 = 0.d0
!  censorstatus = (/.true., .false., .true., .false./)
!  isdisc = (/.false., .false., .false., .false., .false./)
!  call survivalnllk(n, d, parammat1, parammat2, fam, A, datacdfplus, datacdfminus, censorstatus, isdisc, nllk)
!end

subroutine survivalnllk(n, d, parammat1, parammat2, fam, A, datacdfplus, datacdfminus, censorstatus, isdisc, nllk)
  implicit none
  integer n, d, i, j, k, l
  integer fam(d,d), A(d,d), M(d,d), Imat(d,d)
  double precision parammat1(d,d), parammat2(d,d)
  double precision datacdfplus(n,d), datacdfminus(n,d)
  double precision vplus(n,d), vminus(n,d), vpplus(n,d), vpminus(n,d), splus(n,d), sminus(n,d), wplus(n,d), wminus(n,d)
  logical censorstatus(n), isdisc(d)
  double precision nllk, tmp, tmp1, tmp2, tmp3, tmp4
  !print*, n, d
  !print*, parammat1, parammat2
  !print*, fam
  !print*, (A(i,:), i=1,5)
  !print*, datacdfplus
  !print*, datacdfminus
  !print*, censorstatus
  !print*, isdisc
  !print*, nllk
  
  vplus=0.d0; vminus=0.d0; vpplus=0.d0; vpminus=0.d0; splus=0.d0; sminus=0.d0; wplus=0.d0; wminus=0.d0
  M=0; Imat = 0
  do i = 1,d
    do j = i,d
      M(i, j) = maxval(A(1:i, j))
    end do
  end do
  !print*, M
  
  do i = 2,d
    do j = i,d
      if (A(i, j) < M(i, j)) then
        Imat(i-1, M(i, j)) = 1
      end if
    end do
  end do
  !print*, Imat
  
  nllk = 0.d0
  do j = 1,d
    splus(:, j) = datacdfplus(:, A(1, j))
    sminus(:, j) = datacdfminus(:, A(1, j))
    wplus(:, j)= datacdfplus(:, j)
    wminus(:, j) = datacdfminus(:, j)
  end do
  
  do j = 2,d-1
    do k = 1,n
      call mixedcopulallk(splus(k, j), sminus(k, j), wplus(k, j), wminus(k, j), isdisc(A(1, j)), isdisc(j), fam(1, j),&
      parammat1(1, j), parammat2(1, j), tmp)
      nllk = nllk - tmp
    end do
  end do
  !call dblepr("level1", 8, nllk, 1)
  
  do k = 1,n
    if (.not. censorstatus(k)) then
      call mixedcopulallk(splus(k,d), sminus(k,d), wplus(k,d), wminus(k,d), isdisc(A(1,d)), isdisc(d), fam(1,d),&
      parammat1(1,d), parammat2(1,d), tmp)
      nllk = nllk - tmp
    end if
  end do
  !call dblepr("level1", 8, nllk, 1)
  
  do l = 2,d-1
    where (splus >= .999d0)
      splus = .999d0
    end where
    where (sminus >= .999d0)
      sminus = .998d0
    end where
    where (wplus >= .999d0)
      wplus = .999d0
    end where
    where (wminus >= .999d0)
      wminus = .998d0
    end where
    where (splus <= .001d0)
      splus = .002d0
    end where
    where (sminus <= .001d0)
      sminus = .001d0
    end where
    where (wplus <= .001d0)
      wplus = .002d0
    end where
    where (wminus <= .001d0)
      wminus = .001d0
    end where
    
    do j = 1,d
      if (Imat(l-1, j) == 1) then
        if (isdisc(j)) then
          do k = 1,n
            call bivcopcdf(splus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
            call bivcopcdf(splus(k,j), wminus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
            vpplus(k,j) = (tmp1 - tmp2) / (wplus(k,j) - wminus(k,j))
            call bivcopcdf(sminus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp3)
            call bivcopcdf(sminus(k,j), wminus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp4)
            vpminus(k,j) = (tmp3 - tmp4) / (wplus(k,j) - wminus(k,j))
          end do
        else
          do k = 1,n
            call bivcopccdf(splus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
            vpplus(k,j) = tmp1
            call bivcopccdf(sminus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
            vpminus(k,j) = tmp2
          end do
        end if
      end if
      if (isdisc(A(l-1, j))) then
        do k = 1,n
          call bivcopcdf(splus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
          call bivcopcdf(sminus(k,j), wplus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
          vplus(k,j) = (tmp1 - tmp2) / (splus(k,j) - sminus(k,j))
          call bivcopcdf(splus(k,j), wminus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp3)
          call bivcopcdf(sminus(k,j), wminus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp4)
          vminus(k,j) = (tmp3 - tmp4) / (splus(k,j) - sminus(k,j))
        end do
      else
        do k = 1,n
          call bivcopccdf(wplus(k,j), splus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
          vplus(k,j) = tmp1
          call bivcopccdf(wminus(k,j), splus(k,j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
          vminus(k,j) = tmp2
        end do
      end if
    end do
    
    do j = l+1,d
      if (A(l,j) == M(l,j)) then
        splus(:, j) = vplus(:, M(l,j))
        sminus(:, j) = vminus(:, M(l,j))
      else
        splus(:, j) = vpplus(:, M(l,j))
        sminus(:, j) = vpminus(:, M(l,j))
      end if
      wplus(:, j) = vplus(:, j)
      wminus(:, j) = vminus(:, j)
    end do
    
    if (l < (d - 1)) then
      do j = l+1,d-1
        do k = 1,n
          call mixedcopulallk(splus(k,j), sminus(k,j), wplus(k,j), wminus(k,j), isdisc(A(l,j)), isdisc(j), fam(l,j),&
          parammat1(l,j), parammat2(l,j), tmp)
          nllk = nllk - tmp
        end do
      end do
      !call dblepr("levell", 8, nllk, 1)
      do k = 1,n
        if (.not. censorstatus(k)) then
          call mixedcopulallk(splus(k,d), sminus(k,d), wplus(k,d), wminus(k,d), isdisc(A(l,d)), isdisc(d), fam(1,d),&
          parammat1(l,d), parammat2(l,d), tmp)
          nllk = nllk - tmp
        end if
      end do
      !call dblepr("levell", 8, nllk, 1)
    else
      do k = 1,n
        if (.not. censorstatus(k)) then
          call mixedcopulallk(splus(k,d), sminus(k,d), wplus(k,d), wminus(k,d), isdisc(A(d-1,d)), isdisc(d), fam(d-1,d),&
          parammat1(d-1,d), parammat2(d-1,d), tmp)
          nllk = nllk - tmp
        end if
      end do
      !call dblepr("levell", 8, nllk, 1)
      if (isdisc(A(d-1, d))) then
        do k = 1,n
          if (censorstatus(k)) then
            call bivcopcdf(splus(k,d), wplus(k,d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp1)
            call bivcopcdf(sminus(k,d), wplus(k,d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp2)
            nllk = nllk - log(1 - (tmp1 - tmp2) / (splus(k,d) - sminus(k,d)))
          end if
        end do
        !call dblepr("levell", 8, nllk, 1)
      else
        do k = 1,n
          if (censorstatus(k)) then
            call bivcopccdf(wplus(k,d), splus(k,d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp)
            nllk = nllk - log(1 - tmp)
          end if
        end do
        !call dblepr("levell", 8, nllk, 1)
      end if
    end if
  end do
  return
end