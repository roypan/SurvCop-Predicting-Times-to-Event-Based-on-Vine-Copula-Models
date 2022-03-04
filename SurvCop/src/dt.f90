subroutine dt(x,nu, density)
  implicit none
  double precision x, nu, lgamma, pi, tmp, density
  pi = 4.d0 * datan(1.d0)
  tmp = (nu + 1.d0) / 2.d0
  density = exp(lgamma(tmp)) / sqrt(nu * pi) / exp(lgamma(nu / 2.d0)) * (1.d0 + x**2.d0 / nu) ** (-tmp)
  return
end