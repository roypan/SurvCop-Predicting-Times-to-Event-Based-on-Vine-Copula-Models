subroutine mixedcopulallk(uplus, uminus, vplus, vminus, isdiscu, isdiscv, fam, param1, param2, llk)
  implicit none
  integer fam
  double precision uplus, uminus, vplus, vminus, llk, numerator, tmp, tmp1, tmp2, tmp3, tmp4
  double precision param1, param2
  logical isdiscu, isdiscv

  if (uplus >= .999d0) then
    uplus = .999d0
  end if
  if (uminus >= .999d0) then
    uminus = .998d0
  end if
  if (vplus >= .999d0) then
    vplus = .999d0
  end if
  if (vminus >= .999d0) then
    vminus = .998d0
  end if
  if (uplus <= .001d0) then
    uplus = .002d0
  end if
  if (uminus <= .001d0) then
    uminus = .001d0
  end if
  if (vplus <= .001d0) then
    vplus = .002d0
  end if
  if (vminus <= .001d0) then
    vminus = .001d0
  end if
  
  if ((.not. isdiscu) .and. (.not. isdiscv)) then
    call bivcoppdf(uplus, vplus, fam, param1, param2, tmp)
    llk = log(tmp)
  else if (isdiscu .and. (.not. isdiscv)) then
    call bivcopccdf(uplus, vplus, fam, param1, param2, tmp2)
    call bivcopccdf(uminus, vplus, fam, param1, param2, tmp1)
    numerator = tmp2 - tmp1
    llk = log(numerator / (uplus - uminus))
  else if ((.not. isdiscu) .and. isdiscv) then
    call bivcopccdf(vplus, uplus, fam, param1, param2, tmp2)
    call bivcopccdf(vminus, uplus, fam, param1, param2, tmp1)
    numerator = tmp2 - tmp1
    !call dblepr("numerator", 8, numerator,1)
    !call intpr("family",4,fam,1)
    !call dblepr("vplus", 8, vplus,1)
    !call dblepr("vplus", 8, vminus,1)
    llk = log(numerator / (vplus - vminus))
  else if (isdiscu .and. isdiscv) then
    call bivcopcdf(uplus, vplus, fam, param1, param2, tmp4) 
    call bivcopcdf(uminus, vplus, fam, param1, param2, tmp3)
    call bivcopcdf(uplus, vminus, fam, param1, param2, tmp2)
    call bivcopcdf(uminus, vminus, fam, param1, param2, tmp1)
    numerator = tmp4 - tmp3 - tmp2 + tmp1
    llk = log(numerator / (uplus - uminus) / (vplus - vminus))
  end if
  return
end