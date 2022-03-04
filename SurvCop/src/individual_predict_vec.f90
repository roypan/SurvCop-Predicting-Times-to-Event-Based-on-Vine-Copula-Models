subroutine individualpredictvec(n, d, parammat1, parammat2, fam, A, M, Imat, uvec, plus, minus, isdisc, p, output)
  implicit none
  integer n, d, i, j, k, l
  integer fam(d,d), A(d,d), M(d,d), Imat(d,d)
  double precision parammat1(d,d), parammat2(d,d)
  double precision plus(n,d-1), minus(n,d-1), datacdfplus(n,d), datacdfminus(n,d), output(n), uvec(n), ucond(n)
  double precision vplus(n,d), vminus(n,d), vpplus(n,d), vpminus(n,d), splus(n,d), sminus(n,d), wplus(n,d), wminus(n,d)
  logical isdisc(d)
  double precision tmp, tmp1, tmp2, tmp3, tmp4, p
  
  vplus=0.d0; vminus=0.d0; vpplus=0.d0; vpminus=0.d0; splus=0.d0; sminus=0.d0; wplus=0.d0; wminus=0.d0
  
  datacdfplus(:, 1:(d-1)) = plus; datacdfplus(:, d) = uvec
  datacdfminus(:, 1:(d-1)) = minus; datacdfminus(:, d) = uvec
  
  do j = 1,d
    splus(:,j) = datacdfplus(:,A(1, j))
    sminus(:,j) = datacdfminus(:,A(1, j))
    wplus(:,j)= datacdfplus(:,j)
    wminus(:,j) = datacdfminus(:,j)
  end do
  
  if (d == 2) then
    do k = 1,n
      if (isdisc(A(d-1, d))) then
        call bivcopcdf(splus(k,d), wplus(k,d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp1)
        call bivcopcdf(sminus(k,d), wplus(k,d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp2)
        ucond(k) = (tmp1 - tmp2) / (splus(k,d) - sminus(k,d))
      else
        call bivcopccdf(wplus(k,d), splus(k,d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp)
        ucond(k) = tmp
      end if
    end do
  else
    do l = 2,d-1
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
          splus(:,j) = vplus(:,M(l,j))
          sminus(:,j) = vminus(:,M(l,j))
        else
          splus(:,j) = vpplus(:,M(l,j))
          sminus(:,j) = vpminus(:,M(l,j))
        end if
        wplus(:,j) = vplus(:,j)
        wminus(:,j) = vminus(:,j)
      end do
          
      if (l == (d - 1)) then
        if (isdisc(A(d-1, d))) then
          do k = 1,n
            call bivcopcdf(splus(k,d), wplus(k,d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp1)
            call bivcopcdf(sminus(k,d), wplus(k,d), fam(d-1,d), parammat1(d-1,d), parammat2(d-1,d), tmp2)
            ucond(k) = (tmp1 - tmp2) / (splus(k,d) - sminus(k,d))
          end do
        else
          do k =1,n
            call bivcopccdf(wplus(k,d), splus(k,d), fam(d-1, d), parammat1(d-1, d), parammat2(d-1, d), tmp)
            ucond(k) = tmp
          end do
        end if
      end if
    end do
  end if
  output = ucond - p
  return
end