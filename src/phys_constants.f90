module phys_constants
  use kinds, only : r8
  use constants, only : pi, two_pi
  implicit none
  private r8, pi, two_pi
  real(kind=r8), parameter :: vc_si = 2.99792458e8_r8 , vc = vc_si*1e2_r8
  real(kind=r8), parameter :: qe_si = 1.602176487e-19_r8, qe = qe_si*1e1_r8*vc_si
  real(kind=r8), parameter :: me_si = 9.10938215e-31_r8 , me = me_si*1e3_r8
  real(kind=r8), parameter :: me_kev = me_si*vc_si*vc_si/qe_si*1e-3_r8
  real(kind=r8), parameter :: mu0 = 4.e-7_r8*pi
  real(kind=r8), parameter :: cxg19 = 1e-5_r8/pi*qe*qe/me
  real(kind=r8), parameter :: cai_mev = 1e15_r8/(two_pi*vc_si**2)
  real(kind=r8), parameter :: cai_pws = cai_mev*qe_si*1e18_r8
  real(kind=r8), parameter :: k0_per_ghz = 1e9_r8*two_pi/vc
  real(kind=r8), parameter :: bres_per_ghz = 1e5_r8*two_pi*me*vc/qe
  real(kind=r8), parameter :: kboltz_si = 1.3806504e-23_r8, kboltz_evk = kboltz_si/qe_si
end module phys_constants
