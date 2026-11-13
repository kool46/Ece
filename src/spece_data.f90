module spece_data
	use kinds, only : r8
	implicit none
	integer, parameter :: nmx=20000,               &
	                      jmx=12,                  &
	                      kmx=20,                  &
	                      nnw = 500, nnh = 500,    &
	                      nbb = 1000,              &
	                      npmx = 250,              &
	                      ncmx = 10,               &
	                      nfmax=2000,              &
	                      nhlmx=300,               &
	                      naxmx=3000,              &
	                      nscmx=10,                &
	                      ncuts = 100,             &
	                      n_psi = 100,             &
	                      npts=500
	real(kind=r8), parameter :: tmax=5.0_r8
end module spece_data
