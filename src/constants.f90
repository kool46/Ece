module constants
  use kinds, only : i8, r8
  implicit none
  private i8, r8
  real(kind=r8)   , parameter :: zero = 0._r8, one = 1._r8
  real(kind=r8)   , parameter :: two = 2._r8, half = 0.5_r8
  complex(kind=r8), parameter :: im_unit = (zero, one)
  real(kind=r8)   , parameter :: pi     = 3.1415926535897932384626433832795_r8
  real(kind=r8)   , parameter :: two_pi = two*pi , half_pi = half*pi
  real(kind=r8)   , parameter :: pi_sqrd= pi*pi
  real(kind=r8)   , parameter :: pi_sqrt= 1.7724538509055160272981674833411_r8
  real(kind=r8)   , parameter :: deg2rad= pi/180._r8, rad2deg = 180._r8/pi
  real(kind=r8)   , parameter :: huge_r8= huge(one), tiny_r8=tiny(one)
  real(kind=r8)   , parameter :: eps_r8=epsilon(one)
end module constants
