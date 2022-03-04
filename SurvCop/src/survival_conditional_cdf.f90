subroutine survivalconditionalcdf(n, d, u, parammat1, parammat2, fam, A, datacdfplus, datacdfminus, isdisc, output)
  implicit none
  integer n, d, i, j, l
  integer fam(d,d), A(d,d), M(d,d), Imat(d,d)
  double precision parammat1(d,d), parammat2(d,d)
  double precision u(n), datacdfplus(n,d-1), datacdfminus(n,d-1), datacdfplusi(d), datacdfminusi(d), output(n)
  double precision vplus(d), vminus(d), vpplus(d), vpminus(d), splus(d), sminus(d), wplus(d), wminus(d)
  logical isdisc(d)
  double precision ucond, tmp, tmp1, tmp2, tmp3, tmp4
  !print*, n, d
  !print*, parammat1, parammat2
  !print*, fam
  !print*, (A(i,:), i=1,5)
  !print*, datacdfplus
  !print*, datacdfminus
  !print*, isdisc
  
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
  
  do i = 1,n
    datacdfplusi(1:(d-1)) = datacdfplus(i, :)
    datacdfplusi(d) = u(i)
    datacdfminusi(1:(d-1)) = datacdfminus(i, :)
    datacdfminusi(d) = u(i)
    do j = 1,d
      splus(j) = datacdfplusi(A(1, j))
      sminus(j) = datacdfminusi(A(1, j))
      wplus(j)= datacdfplusi(j)
      wminus(j) = datacdfminusi(j)
    end do
    
    if (d == 2) then
      if (isdisc(A(d-1, d))) then
        call bivcopcdf(splus(d), wplus(d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp1)
        call bivcopcdf(sminus(d), wplus(d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp2)
        ucond = (tmp1 - tmp2) / (splus(d) - sminus(d))
      else
        call bivcopccdf(wplus(d), splus(d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp)
        ucond = tmp
      end if
    else
      do l = 2,d-1
        where (splus >= 1.d0)
          splus = .999d0
        end where
        where (sminus >= .999d0)
          sminus = .998d0
        end where
        where (wplus >= 1.d0)
          wplus = .999d0
        end where
        where (wminus >= .999d0)
          wminus = .998d0
        end where
        where (splus <= .001d0)
          splus = .002d0
        end where
        where (sminus <= 0.d0)
          sminus = .001d0
        end where
        where (wplus <= .001d0)
          wplus = .002d0
        end where
        where (wminus <= 0.d0)
          wminus = .001d0
        end where
        
        do j = 1,d
          if (Imat(l-1, j) == 1) then
            if (isdisc(j)) then
              call bivcopcdf(splus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
              call bivcopcdf(splus(j), wminus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
              vpplus(j) = (tmp1 - tmp2) / (wplus(j) - wminus(j))
              call bivcopcdf(sminus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp3)
              call bivcopcdf(sminus(j), wminus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp4)
              vpminus(j) = (tmp3 - tmp4) / (wplus(j) - wminus(j))
            else
              call bivcopccdf(splus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
              vpplus(j) = tmp1
              call bivcopccdf(sminus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
              vpminus(j) = tmp2
            end if
          end if
          if (isdisc(A(l-1, j))) then
            call bivcopcdf(splus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
            call bivcopcdf(sminus(j), wplus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
            vplus(j) = (tmp1 - tmp2) / (splus(j) - sminus(j))
            call bivcopcdf(splus(j), wminus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp3)
            call bivcopcdf(sminus(j), wminus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp4)
            vminus(j) = (tmp3 - tmp4) / (splus(j) - sminus(j))
          else
            call bivcopccdf(wplus(j), splus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp1)
            vplus(j) = tmp1
            call bivcopccdf(wminus(j), splus(j), fam(l-1,j), parammat1(l-1,j), parammat2(l-1,j), tmp2)
            vminus(j) = tmp2
          end if
        end do
        
        do j = l+1,d
          if (A(l,j) == M(l,j)) then
            splus(j) = vplus(M(l,j))
            sminus(j) = vminus(M(l,j))
          else
            splus(j) = vpplus(M(l,j))
            sminus(j) = vpminus(M(l,j))
          end if
          wplus(j) = vplus(j)
          wminus(j) = vminus(j)
        end do
        
        if (l == (d - 1)) then
          if (isdisc(A(d-1, d))) then
            call bivcopcdf(splus(d), wplus(d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp1)
            call bivcopcdf(sminus(d), wplus(d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp2)
            ucond = (tmp1 - tmp2) / (splus(d) - sminus(d))
          else
            call bivcopccdf(wplus(d), splus(d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp)
            ucond = tmp
          end if
        end if
      end do
    end if
    output(i) = ucond
  end do
  return
end