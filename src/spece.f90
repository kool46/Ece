#ifndef REVISION
#define REVISION 'unknown'
#endif
program spece
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!           SPECE: Raytracing code for Tokamak ECE modelling           C
!           ------------------------------------------------           C
!                    D.Farina, L.Figini, P.Platania                    C
!           IFP-CNR Milano, v. Cozzi 53, 20125 Milano, Italy           C
!                     EURATOM-ENEA-CNR Association                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  use kinds, only : r8
  use constants, only : pi
  use phys_constants, only : qe, me, vc
  use spece_data, only : nmx
  implicit none
  integer nstep,nom,iox,nray,ktx,ierr,iom,j,k,kktx,narg,iargc
  integer istop,istep,i,iprtr,iprfr,it
  integer lid,lcase
  real(kind=r8) st0,st,dst,ak0,fghz0,fghzf,dfghz,fghz,bres,sox,ddrmx
  character(len=255) idstring,prefix,sarg
  logical prefnotset,alloutput,goutput

  common/ss/st
  common/init/nstep
  common/dsds/dst
  common/ak0/ak0
  common/fghz0f/fghz0,fghzf
  common/fghz/fghz
  common/nom/nom
  common/parbres/bres
  common/mode/sox
  common/imode/iox
  common/ierr/ierr
  common/ddrmx/ddrmx
  common/nrayktx/nray,ktx
  common/spectinfo/idstring,prefix
  common/iprint/iprtr,iprfr
  common/out/alloutput,goutput
  
  interface
    subroutine read_param
    end subroutine read_param
    subroutine read_data
      use kinds, only: r8
    end subroutine read_data
    subroutine prfile
    end subroutine prfile
    subroutine print_current
    end subroutine print_current
    subroutine pweight
    end subroutine pweight
    subroutine datinit
    end subroutine datinit
    subroutine icrtrn(j,k,st0,ierr)
      use kinds, only : r8
      integer, intent(in) :: j, k
      real(kind=r8), intent(out) :: st0
      integer, intent(out) :: ierr
    end subroutine icrtrn
!    subroutine pol_limit
!    end subroutine pol_limit
    subroutine rkint4
    end subroutine rkint4
    subroutine after_onestep(istep,istop)
      integer, intent(in)  :: istep
      integer, intent(out) :: istop
    end subroutine after_onestep
    subroutine rayenddata
    end subroutine rayenddata
    subroutine raysections
    end subroutine raysections
    subroutine taonaxis(j,k,iox)
      integer, intent(in) :: j, k, iox
    end subroutine taonaxis
    subroutine sub_trad(j,k)
      integer, intent(in) :: j, k
    end subroutine sub_trad
    subroutine trad_back
    end subroutine trad_back
    subroutine printtraj
    end subroutine printtraj
    subroutine  trad_tot(j,k)
      integer, intent(in) :: j, k
    end subroutine  trad_tot
    subroutine trad_sum
    end subroutine trad_sum
    subroutine taaxis_sum(goutput)
      logical, intent(in) :: goutput
    end subroutine taaxis_sum
  end interface

  prefnotset = .true.
  alloutput = .true.
  goutput = .false.
  narg=iargc()
  if (narg>0) then
    do i=1,narg
      call getarg(i,sarg)
      if (sarg(1:1)=='-') then
        do j=2,len_trim(sarg)
          select case(sarg(j:j))
          case('h','H')
            write(*,*) 'Usage: spece [options] [target dir]'
            write(*,*) 'Options: -h    show this help'
            write(*,*) '         -g    enables output of G function on axis'
            write(*,*) '         -s    final spectra output only'
            write(*,*) '         -v    show version number'
            write(*,*) 'If "target dir" is present, output files will be located there'
            write(*,*) 'Input files needed:'
            write(*,*) '         1) spece.data: input file names and code parameters'
            write(*,*) '         2) magnetic equilibrium description in G EQDSK format'
            write(*,*) '         3) electron density profile vs normalized psi'
            write(*,*) '                        density in 10^19 m^-3'
            write(*,*) '         4) electron temperature profile vs normalized psi'
            write(*,*) '                        temperature in keV'
            write(*,*) '         5) electron distribution function parameters'
            write(*,*) '         Name of files 2)-5) must be given in spece.data'
            stop
          case('g','G')
            write(*,*) 'G function output enabled'
            goutput = .true.
          case('s','S')
            alloutput = .false.
            write(*,*) 'Final spectra output only'
          case('v','V')
            write(*,*) '  SVN Rev. '//REVISION
            stop
          end select
        end do
      else if (prefnotset) then
        if (scan(sarg,'/',back=.true.)==len_trim(sarg)) then
          write(prefix,'(a)') trim(sarg)
        else
          write(prefix,'(a)') trim(sarg)//'/'
        end if
        prefnotset=.false.
        write(*,*) 'Output files in '//trim(prefix)
      else
        write(*,*) 'Unexpected argument: '//trim(sarg)
        write(*,*) 'Type spece -h for usage'
        stop
      end if
    end do
  end if
  if (prefnotset) write(prefix,'(a)') ''
  call read_param
! read data plus initialization
  call read_data
  if (nstep>nmx) nstep=nmx
  open(15,file=trim(prefix)//'spectra'//trim(idstring)//'.txt',action='write')
  if (alloutput) then
    open(11,file=trim(prefix)//'O'//trim(idstring)//'.txt',action='write')
    open(12,file=trim(prefix)//'X'//trim(idstring)//'.txt',action='write')
    open(17,file=trim(prefix)//'err'//trim(idstring)//'.txt',action='write')
    if (iprfr>0) then
      open(21,file=trim(prefix)//'engyO'//trim(idstring)//'.txt',action='write')
      open(22,file=trim(prefix)//'engyX'//trim(idstring)//'.txt',action='write')
      if (iprtr>0) then
        open( 9,file=trim(prefix)//'trajO'//trim(idstring)//'.txt',action='write')
        open(10,file=trim(prefix)//'trajX'//trim(idstring)//'.txt',action='write')
      end if
    end if
    open(23,file=trim(prefix)//'resolO'//trim(idstring)//'.txt',action='write')
    open(24,file=trim(prefix)//'resolX'//trim(idstring)//'.txt',action='write')
    open(25,file=trim(prefix)//'centray'//trim(idstring)//'.txt',action='write')
    if (goutput) open(26,file=trim(prefix)//'Gaxis'//trim(idstring)//'.txt',action='write')
  end if
  call prfile
  if (alloutput) call print_current
  call pweight
! loop on frequency: begin
  dfghz=0.0d0
  if(nom>1) dfghz=(fghzf-fghz0)/dble(nom-1)
  do iom=1,nom
    fghz=fghz0+dfghz*(iom-1)
    ak0=2.0d9*pi*fghz/vc
    bres=2.0d+5*pi*fghz*me*vc/qe
!   loop on rays: begin
    do j=1,nray
      kktx=1
      if(j>1) kktx=ktx
      do k=1,kktx
!       loop on OM-XM: begin
        do iox=1,2
          if(iox==1) then
            sox=-1.0d0
          else
            sox=1.0d0
          end if
          ierr=0
          call datinit
          call icrtrn(j,k,st0,ierr)
          if (ierr==0) then
!            if(iox==1) call pol_limit
!           ray integration: begin
!            st0=st0-dst
            st0=-dst
            do i=2,nstep
              st=st0+i*dst
              istep=i
!             one step integration:
              call rkint4
!             calculations after one step:
              call after_onestep(istep,istop)
              if(istop==1) exit
            end do
!            write(33,*) ' '
!           ray integration: end
          end if
          call rayenddata
          call raysections
          call taonaxis(j,k,iox)
          call sub_trad(j,k)
          call trad_back
          if(alloutput.and.((j==1.and.k==1).or.nom==1)) then
            write(17,98) iox,fghz,st,ddrmx
            if (iprfr>0.and.mod(iom-1,iprfr)==0) call printtraj
          end if
        end do
!       loop on OM-XM: end
        call  trad_tot(j,k)
      end do
    end do
!   loop on rays: end
    call trad_sum
    if (alloutput) call taaxis_sum(goutput)
  end do
! loop on frequency: end
  close(15)
  if (alloutput) then
    close(11)
    close(12)
    close(17)
    if (iprfr>0) then
      close(21)
      close(22)
      if (iprtr>0) then
        close(9)
        close(10)
      end if
    end if
    close(23)
    close(24)
    close(25)
    if (goutput) close(26)
  end if
  98 format(i5,20(1x,e16.8e3))
end program spece

subroutine los_inout(r_min,r_max,z_min,z_max,x0,y0,z0,nx,ny,nz, &
                     tau_in,tau_out,ierr)
! finds entrance/exit for a given valid R,z region and a line of sight
! line described parametrically as (x,y,z) = (x0,y0,z0) - tau (nx,ny,nz)
  use kinds, only : r8
  use constants, only : zero,eps_r8, huge_r8, deg2rad
  use intervals, only : interval, intersect, assignment(=)
  implicit none
! arguments
  real(kind=r8), intent(in)  :: r_min, r_max, z_min, z_max
  real(kind=r8), intent(in)  :: x0, y0, z0, nx, ny, nz
  real(kind=r8), intent(out) :: tau_in, tau_out
  integer      , intent(out) :: ierr
! local variables
  real(kind=r8) :: n_norm2, r02, nxy2, xydotn
  real(kind=r8) :: r_staz, tau_staz, delta_r_dn, delta_r_up
  type(interval) :: tau_z, tau_r_dn(2), tau_r_up, tau(2), tau_temp, tau_ok

  ierr = 0
  r02= x0**2+y0**2
  nxy2 = nx**2 + ny**2 ! squared horizontal component
  n_norm2 = nxy2 + nz**2
  if (n_norm2>zero) then
    xydotn = -x0*nx - y0*ny
! coordinate of intersection with inner and outer R boundaries
    if (nxy2>zero) then
      tau_staz = -xydotn/nxy2
      r_staz = sqrt(r02+tau_staz*xydotn)
      if (r_staz<r_max) then    ! line intersects outer boundary
        delta_r_up = sqrt(xydotn**2-nxy2*(r02-r_max**2))
        tau_r_up%lo = (-xydotn-delta_r_up)/nxy2
        tau_r_up%hi = (-xydotn+delta_r_up)/nxy2
      else
        tau_r_up = -huge_r8
        ierr = ierr + 100
      end if
      tau_r_dn(1)%lo = -huge_r8
      tau_r_dn(2)%hi =  huge_r8
      if (r_staz<r_min) then    ! line intersects inner boundary
        delta_r_dn = sqrt(xydotn**2-nxy2*(r02-r_min**2))
        tau_r_dn(1)%hi = (-xydotn-delta_r_dn)/nxy2
        tau_r_dn(2)%lo = (-xydotn+delta_r_dn)/nxy2
      else
        tau_r_dn(1)%hi = -huge_r8
        tau_r_dn(2)%lo = -huge_r8
        ierr = ierr + 10
      end if
    else                            ! vertical line
      r_staz = sqrt(r02)
      if (r_staz<r_max) then
        tau_r_up = interval(-huge_r8,huge_r8)
      else
        tau_r_up = interval(-huge_r8,-huge_r8)
      end if
      tau_r_dn(1) = -huge_r8
      if (r_staz>r_min) then
        tau_r_dn(2) = interval(-huge_r8,huge_r8)
      else
        tau_r_dn(2) = interval(-huge_r8,-huge_r8)
      end if
    end if
! coordinate of intersection with upper and lower z boundaries
    if (nz>0) then                           ! 2 intersections exist
      tau_z%lo = (z_min-z0)/nz
      tau_z%hi = (z_max-z0)/nz
    else if (nz<0) then
      tau_z%lo = (z_max-z0)/nz
      tau_z%hi = (z_min-z0)/nz
    else if (z_max>z0 .and. z_min<z0) then   ! between the 2 planes
      tau_z = interval(-huge_r8, huge_r8)
    else                                     ! outside the 2 planes
      tau_z = interval(-huge_r8,-huge_r8)
    end if
! intersection of valid intervals
    tau_temp = interval(0._r8,huge_r8)
    tau_temp = intersect(tau_temp,tau_z)
    tau_temp = intersect(tau_temp,tau_r_up)
    tau(1) = intersect(tau_temp,tau_r_dn(1))
    tau(2) = intersect(tau_temp,tau_r_dn(2))
    if (tau(1)%lo >= 0._r8) then
      tau_ok = tau(1)
    else
      tau_ok = tau(2)
    end if
    if (tau_ok%lo < 0._r8) then
      tau_in=0._r8
      tau_out=0._r8
      ierr = ierr + 1
      print *, 'Line of sight does not intersect computational box'
    else
      tau_in = tau_ok%lo
      tau_out = tau_ok%hi
      ierr = 0
    end if
  else
    print *, 'Invalid direction vector'
    ierr = 222
  end if
end subroutine los_inout

subroutine goto_bound(x0, y0, z0, nx, ny, nz, st, ierr)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! GOTO_BOUND computation of ray intersection with plasma boundary and  C
!            re-setting of initial conditions                          C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : zero, tiny_r8, eps_r8
  implicit none
! arguments
  real(kind=r8), intent(inout) :: x0, y0, z0
  real(kind=r8), intent(in)    :: nx, ny, nz
  real(kind=r8), intent(out)   :: st
  integer,       intent(out)   :: ierr
! local variables
  integer       :: ier
  real(kind=r8) :: r_min, r_max, z_min, z_max
  real(kind=r8) :: tau_in, tau_out, tau_psimin, tau_bound, tau_mid
  real(kind=r8), parameter :: eps = 1e-8_r8, eta = 1e-8_r8
! common variables
  real(kind=r8) :: rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz
  common/pareq2/rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz
! external subroutines
  interface
    subroutine los_inout(r_min,r_max,z_min,z_max,x0,y0,z0,nx,ny,nz, &
                     tau_in,tau_out,ierr)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: r_min, r_max, z_min, z_max
      real(kind=r8), intent(in)  :: x0, y0, z0, nx, ny, nz
      real(kind=r8), intent(out) :: tau_in, tau_out
      integer      , intent(out) :: ierr
    end subroutine los_inout
    subroutine zeroin(a,b,eps,eta,f,x,ifail)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: a, b, eps, eta
      real(kind=r8), intent(out)   :: x
      integer      , intent(inout) :: ifail
      real(kind=r8), external :: f
    end subroutine zeroin
    function dpsi_dtau(tau)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: dpsi_dtau
      real(kind=r8), intent(in) :: tau
    end function dpsi_dtau
    function dpsi_bound(tau)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: dpsi_bound
      real(kind=r8), intent(in) :: tau
    end function dpsi_bound
  end interface
! executable statements
! grid boundaries in cm
  r_min = (rmnm+0.1_r8*dr)*100._r8
  r_max = (rmxm-0.1_r8*dr)*100._r8
  z_min = (zmnm+0.1_r8*dz)*100._r8
  z_max = (zmxm-0.1_r8*dz)*100._r8
  call los_inout(r_min,r_max,z_min,z_max,x0,y0,z0,nx,ny,nz, &
                     tau_in,tau_out,ier)
  if (mod(ier,10)/=0) then
    ! intersection with grid boundaries not found
    print*, 'intersection with grid boundaries not found'
    st = 0
    ierr = 1
  else
    tau_mid=(tau_out+tau_in)/2._r8
    if (mod(ier,100)/=0) then
      ! line does not intersect inner boundary
      print*, 'line does not intersect inner boundary'
      tau_mid=tau_mid/2._r8
!      if (dpsi_dtau(tau_mid)*dpsi_dtau(tau_in)<zero) then
!        call zeroin(tau_in,tau_mid,eps,eta,dpsi_dtau,tau_psimin,ier)
!      else if (dpsi_dtau(tau_mid)*dpsi_dtau(tau_out)<zero) then
!        call zeroin(tau_mid,tau_out,eps,eta,dpsi_dtau,tau_psimin,ier)
!      else
!        tau_psimin = tau_mid
!      end if
!    else
!      call zeroin(tau_in,tau_out,eps,eta,dpsi_dtau,tau_psimin,ier)
    end if
!    if (dpsi_bound(tau_psimin)>zero) then
    if (dpsi_bound(tau_mid)>zero) then
      ! line does not intersect plasma boundary
      print*, 'line does not intersect plasma'
      print*, tau_in, tau_out, tau_mid
      x0 = x0 - tau_mid*nx
      y0 = y0 - tau_mid*ny
      z0 = z0 - tau_mid*nz
      print*, x0,y0,z0
      st = 0
      ierr = 1
    else
      if (dpsi_bound(tau_in)<zero) then
        ! antenna inside plasma boundary
        print*, 'warning: antenna inside plasma boundary'
        tau_bound = tau_in
      else
!        call zeroin(tau_in,tau_psimin,1e-10_r8,sqrt(tiny_r8),dpsi_bound,tau_bound,ier)
        call zeroin(tau_in,tau_mid,1e-10_r8,sqrt(tiny_r8),dpsi_bound,tau_bound,ier)
      end if
      x0 = x0 - tau_bound*nx
      y0 = y0 - tau_bound*ny
      z0 = z0 - tau_bound*nz
      st = tau_bound*sqrt(nx**2+ny**2+nz**2)
      ierr = 0
    end if
  end if
end subroutine goto_bound

function dpsi_dtau(tau)
  use kinds, only : r8
  implicit none
  real(kind=r8) :: dpsi_dtau
  real(kind=r8), intent(in) :: tau
  real(kind=r8) :: x0, y0, z0, taum, nx, ny, nz, rm, zm, xm, ym
  integer, parameter :: ndim = 6
  real(kind=r8), dimension(ndim) :: ywrk0,ypwrk0
  real(kind=r8) :: dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/wrk/ywrk0,ypwrk0
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface
  
  x0 =.01_r8*ywrk0(1)
  y0 =.01_r8*ywrk0(2)
  z0 =.01_r8*ywrk0(3)
  taum =.01_r8*tau
  nx = ywrk0(4)
  ny = ywrk0(5)
  nz = ywrk0(6)
  xm = x0-nx*taum
  ym = y0-ny*taum
  rm = sqrt(xm**2+ym**2)
  zm = z0-nz*tau
  call equinum(rm,zm)
  dpsi_dtau = -dpsidr*(xm*nx+ym*ny)/rm-dpsidz*nz
end function dpsi_dtau

function dpsi_bound(tau)
  use kinds, only : r8
  implicit none
  real(kind=r8) :: dpsi_bound
  real(kind=r8), intent(in) :: tau
  real(kind=r8) :: x0, y0, z0, taum, nx, ny, nz, rm, xm, ym, zm
  integer, parameter :: ndim = 6
  real(kind=r8), dimension(ndim) :: ywrk0,ypwrk0
  real(kind=r8) :: psinv, psdbnd
  common/wrk/ywrk0,ypwrk0
  common/psival/psinv
  common/densbnd/psdbnd
  
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface
  
  x0 =.01_r8*ywrk0(1)
  y0 =.01_r8*ywrk0(2)
  z0 =.01_r8*ywrk0(3)
  taum =.01_r8*tau
  nx = ywrk0(4)
  ny = ywrk0(5)
  nz = ywrk0(6)
  xm = x0-nx*taum
  ym = y0-ny*taum
  rm = sqrt(xm**2+ym**2)
  zm = z0-nz*tau
  call equinum(rm,zm)
  dpsi_bound = psinv-psdbnd
end function dpsi_bound

subroutine after_onestep(i,istop)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! AFTER_ONESTEP(i,istop) computation absorption/emission at ith step   C
!   uses                                                               C
!   sets                                                               C
!   calls EMISSION                                                     C
!         DENSITY                                                      C
!         POL_LIMIT                                                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi
  use spece_data, only : nmx
  implicit none
  integer, intent(in)  :: i
  integer, intent(out) :: istop
  integer iox,iif,ierr,iwarm,ilarm,ipolc,iproj
  real(kind=r8) taumn,taumx,psinv,psdbnd,dd,an2s,an2,ddrmx
  real(kind=r8) rr,rrm,yg,dens,ddens,st,xg
  real(kind=r8) anpl,anpl2,anpr
  integer      , parameter :: ndim = 6
  real(kind=r8) anplv(nmx),psjki(nmx)
  real(kind=r8) tauv(nmx),alphav(nmx)
  real(kind=r8) tev(nmx),rrr(nmx),ygv(nmx)
  real(kind=r8) ywrk(ndim),ypwrk(ndim)
  real(kind=r8) zzv(nmx)
  real(kind=r8) stv(nmx),phiv(nmx),xgv(nmx),anprv(nmx),ddv(nmx),ddrmxv(nmx)

  common/atjki/tauv,alphav,tev
  common/ygv/ygv
  common/warm/iwarm,ilarm
  common/imode/iox
  common/iiv/iif
  common/psjki/psjki
  common/psival/psinv
  common/ierr/ierr
  common/taumnx/taumn,taumx
  common/rrr/rrr
  common/rr/rr
  common/densbnd/psdbnd
  common/ygyg/yg
  common/ddd/dd,an2s,an2
  common/ddrmx/ddrmx
  common/nplr/anpl,anpl2,anpr
  common/nplv/anplv
  common/ss/st
  common/wrk/ywrk,ypwrk
  common/xgxg/xg
  common/traj/stv,zzv,phiv,xgv,anprv
  common/ddv/ddv,ddrmxv

  common/dens/dens,ddens
  common/ipol/ipolc
  common/iproj/iproj

  interface
    integer function islhfs(r,z)
      real(kind=8), intent(in) :: r, z
    end function islhfs
    subroutine emission(i)
      integer, intent(in) :: i
    end subroutine emission
    subroutine pol_limit
    end subroutine pol_limit
  end interface

  taumn=1d+30
  taumx=-1d+30
  istop=0
  if(ierr>0) then
     print*,'            IERR = ', ierr
     istop=1
  end if
  rrr(i)=rr
  rrm=rr*1.d-2
  zzv(i)=ywrk(3)*1.0d-2
  stv(i)=st*1.0d-2
  phiv(i)=acos(ywrk(1)/dsqrt(ywrk(1)**2+ywrk(2)**2))
  if(ywrk(2)<0.0d0) phiv(i)=-phiv(i)
  phiv(i)=phiv(i)*180.0d0/pi
  if(iproj.eq.1) then
    if(islhfs(rrm,zzv(i))==1) then
      psjki(i)=psinv
    else
      psjki(i)=-psinv
    end if
  else
    psjki(i)=psinv
  end if
  xgv(i)=xg
  ygv(i)=yg
  anplv(i)=anpl
  anprv(i)=anpr
  ddv(i)=dd
  if (iwarm>0.and.i>1) then
    if(psinv<=4.9_r8 .and. psinv<=psdbnd) then
      call emission(i)
      iif=i
    else if(iif>1) then
      tauv(i)=tauv(i-1)
      istop=1
    end if
    if(tauv(i)<=taumn) taumn=tauv(i)
    if(tauv(i)>=taumx) taumx=tauv(i)
  end if
  if(abs(dd)>ddrmx .and. istop==0) ddrmx=abs(dd)
  ddrmxv(i)=ddrmx
  call density(psinv)
  if(ipolc==0 .and. dens>0._r8 .and. iox==1) then
    call pol_limit
    ipolc=1
  end if
end subroutine after_onestep

subroutine prfile
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PRFILE print head of files for different subroutines                 C
!   uses                                                               C
!   sets                                                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  implicit none
  integer nom,nlin,iprtr,iprfr
  logical alloutput,goutput
  common/nom/nom
  common/out/alloutput,goutput
  common/iprint/iprtr,iprfr
  nlin=nom
  write(15,*) nlin
  write(15,*) ' #fghz sumtradx sumtrady citxtot citytot'
  if (alloutput) then
    write(11,*) nlin
    write(12,*) nlin
    write(23,*) nlin
    write(24,*) nlin
    write(25,*) nlin
    nlin=2*nlin
    write(17,*) nlin
    if (iprfr>0) then
      write(21,*) ' #fghz s R N// tau alfa T G Emin Emax nharm'
      write(22,*) ' #fghz s R N// tau alfa T G Emin Emax nharm'
      if (iprtr>0) then
        write(9,*)  ' #fghz s R z phi X Y Nperp N// psi dd ddmx'
        write(10,*) ' #fghz s R z phi X Y Nperp N// psi dd ddmx'
      end if
    end if
    write(11,*) ' #fghz sf Tradc Int Rmean DRmn tau iovlp nh_mn Rres'//&
                ' R50 R2.5 R97.5 z50 z2.5 z97.5 psi50 psi2.5 psi97.5'//&
                ' T50 T2.5 T97.5 Tres2 Tres3 Tres4 tau2 tau3 tau4'
    write(12,*) ' #fghz sf Tradc Int Rmean DRmn tau iovlp nh_mn Rres'//&
                ' R50 R2.5 R97.5 z50 z2.5 z97.5 psi50 psi2.5 psi97.5'//&
                ' T50 T2.5 T97.5 Tres2 Tres3 Tres4 tau2 tau3 tau4'
    write(23,*) ' #fghz Rax50 Rax2.5 Rax97.5 Rax50mix Rax2.5mix Rax97.5mix'//&
                ' psi50 psi2.5 psi97.5 rhotor50 rhotor2.5 rhotor97.5'
    write(24,*) ' #fghz Rax50 Rax2.5 Rax97.5 Rax50mix Rax2.5mix Rax97.5mix'//&
                ' psi50 psi2.5 psi97.5 rhotor50 rhotor2.5 rhotor97.5'
    if (goutput) write(26,*) ' #fghz R GaxO GaxX GaxqO GaxqX'
    write(25,"(a)") ' #fghz Tradx11 Trady11 citx11 city11 trqom trqxm'//&
                    ' TradO TradX tauvmO tauvmX TrObk TrXbk cpol cutoff'
    write(17,*) ' #iox fghz stf ddrmx'
!    write(33,*) ' #fghz xg yg mu Npl Nprf2 Npr2 i'
  end if
end subroutine prfile

subroutine read_param
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! READ_PARAM reading of input data and parameters                      C
!   uses                                                               C
!   sets                                                               C
!   calls PROFDATA                                                     C
!         EQUIDATA                                                     C
!         CUTOFF                                                       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : tiny_r8, pi
  implicit none
  integer versign,iprtr,iprfr,ncmp0,fddtype
  character(len=255) idstring,prefix
  character(len=255) filenmeqq,filenmte,filenmne,filenmlos,filenmtail
  integer nom,nray,ktx,nstep,iwarm,ilarm,imx,ixp,iproj
  integer ilos,ierr
  integer ipsinorm,icocos,idsc,iffmt,iscatter,iang
  real(kind=r8)  sgnbphi,sgniphi
  real(kind=r8) fghz0,fghzf,x00,y00,z00,anx0c,any0c,anz0c,rmx
  real(kind=r8) w0xt,pw0xt,w0yt,pw0yt,awr,refl,polar,anglp,anglt,cvdr,polrad,torrad
  real(kind=r8) psdbnd,bfact,bvfact,sspld,dst,sspl
  real(kind=r8) nfact,tfact

  common/fghz0f/fghz0,fghzf
  common/nom/nom
  common/nrayktx/nray,ktx
  common/rhomx/rmx
  common/rp/refl,polar
  common/init/nstep
  common/warm/iwarm,ilarm
  common/idisp/imx
  common/dsds/dst
  common/sspldn/sspld
  common/ixp/ixp
  common/eqconv/ipsinorm,icocos,sgnbphi,sgniphi
  common/ssplps/sspl
  common/anic/anx0c,any0c,anz0c
  common/mirr/x00,y00,z00
  common/parbeam/w0xt,pw0xt,w0yt,pw0yt,awr
  common/densbnd/psdbnd
  common/bfact/bfact,bvfact
  common/tfact/tfact
  common/nfact/nfact
  common/spectinfo/idstring,prefix
  common/eqlfnm/filenmeqq,filenmte,filenmne
  common/vs/versign
  common/ncmp0/ncmp0
  common/fddtyp/fddtype
  common/iprint/iprtr,iprfr
  common/eqdskfmt/idsc,iffmt,iscatter
  common/iproj/iproj

  interface
    subroutine taildata(filenmtail)
      character(len=255), intent(in) :: filenmtail
    end subroutine taildata
  end interface

  open(2,file='spece.data',status= 'old')
! identificative string to put in output file names
  read(2,*) idstring
  idstring=adjustl(idstring)
!   filenmeqq   input magnetic equilibrium file name (G EQDSK format)
!   filenmne    input ne profile file name (10^19 m^-3 units, vs normalized psi)
!   filenmte    input Te profile file name (keV units, vs normalized psi)
  read(2,*) filenmeqq
  read(2,*) filenmne
  read(2,*) filenmte
!   fddtype     distribution function model 0: maxwellian, 1: LH-like edf
!               2,3: multi-maxwellian with absolute/relative parameters
!               (>0: read parameters from filenmtail)
  read(2,*) fddtype,filenmtail
!   ilos        0: perpendicular line of sight
!               other: line of sight read from filenmlos
  read(2,*) ilos,filenmlos
!   fghz        minimum/maximum wave frequency (GHz), number of freq. points
  read(2,*) fghz0,fghzf,nom
!   nray        radial number of rays
!   ktx         angular number of rays
!   rmx         normalized maximum radius of beam power
!               rmx=1: last shell of ray at radius = beam size (asymptotic)
  read(2,*) nray,ktx,rmx
!   beam parameters in local reference system
!   w0t         waist size (cm)
!   pw0         waist distance from launching point (cm)
!   awr         angle of beam ellipse (deg)
  read(2,*) w0xt,w0yt,pw0xt,pw0yt,awr
!   refl,polar  wall reflection coefficients
  read(2,*) refl,polar
!   iwarm=0     no absorption / emission
!   iwarm=1     weakly relativistic absorption
!   iwarm=2     relativistic absorption, n<1 asymptotic expansion
!   ilarm       order of larmor expansion
  read(2,*) iwarm,ilarm
!   imx         max iterations for dispersion relation convergence
!   versign     1: standard sign choice algorithm in dispersion relation
!               2: use previous solution as initial guess
  read(2,*) imx,versign
!   dst         integration step in c*t (cm)
!   nstep       maximum number of integration steps
  read(2,*) dst,nstep
!   sspld       spline parameter for density interpolation
!   psbnd       value of psi ( > 1 ) of density boundary
  read(2,*) sspld,psdbnd
!   sspl        spline parameter for psi interpolation
!   ixp         0: no X point, -1: lower, +1: upper, +2: double X point
  read(2,*) sspl,ixp,ipsinorm,icocos,sgnbphi,sgniphi
!   bfact bvfact  total and vacuum magnetic field scaling factors
!                 B' = bfact*sqrt[B^2 + (bvfact^2-1)*B_vac^2]
!   tfact       scaling factor for temperature
!   nfact       scaling factor for density
  read(2,*) bfact, bvfact, nfact, tfact
  if (bfact<=0.0d0) bfact=1.0d0
  if (bvfact<=0.0d0) bvfact=1.0d0
  if (tfact<=0.0d0) tfact=1.0d0
  if (nfact<=0.0d0) nfact=1.0d0
!   iprtr       raytracing steps between trajectory data output
!   iprfr       frequency steps between trajectory/energy data output
  read(2,*) iprtr,iprfr
!   idsc       id string in 1st EQDSK line (0: no, 1: yes) [default=1]
!   iffmt      EQDSK format (0: fixed, 1: free)            [default=0]
!   iscatter   valid psi data outside LCFS (0: yes, 1: no) [default=0]
  idsc=1
  iffmt=0
  iscatter=0
  read(2,*) idsc,iffmt,iscatter
!  iang        angles convention (0: ITER [alpha beta], 1: TCV [theta phi]) [default=1]
  iang=0
  read(2,*) iang
!  iproj       0:projection on 0<psi<+1, 1:projection on -1<psi<+1
  read(2,*) iproj
  close(unit=2)
  call taildata(filenmtail)
!   x00,y00,z00 coordinates of launching point
!   anglp,anglt considered positive clockwise (when looking towards
!               positive y and negative z from a point (x>0,0,0) )
  select case (ilos)
  case (3)
    x00=413.345d0
    z00=24.8d0
    anglp=0.0d0
    anglt=0.0d0
    write(*,*) 'JET KK3 new line of sight (JPN>80289)'
  case (-3)
    x00=413.345d0
    z00=13.3d0
    anglp=0.0d0
    anglt=0.0d0
    write(*,*) 'JET KK3 old line of sight (JPN<=80289)'
  case (-1)
    x00=413.345d0
    z00=36.288d0
    anglp=0.0d0
    anglt=0.0d0
    write(*,*) 'JET KK1 line of sight'
  case (0)
    x00=437.802d0
    z00=25.497d0
    anglp=0.0d0
    anglt=0.0d0
    write(*,*) 'JET KK5 0 deg line of sight'
  case (1)
    x00=435.45d0
    z00=30.26d0
    anglp=0.073d0
    anglt=9.46d0
    write(*,*) 'JET KK5 10 deg line of sight'
  case (2)
    x00=433.56d0
    z00=30.26d0
    anglp=0.073d0
    anglt=20.59d0
    write(*,*) 'JET KK5 20 deg line of sight'
  case (10)
    x00=240.0d0
    z00=0.0d0
    anglp=0.0d0
    anglt=0.0d0
    write(*,*) 'EAST ECE16/32 line of sight'
  case default
    open(4,file=trim(adjustl(filenmlos)),status='old',action='read')
    read(4,*) x00,z00
    read(4,*) anglp,anglt
    close(4)
    write(*,*) 'Custom line of sight'
    if(iang.eq.1) then ! TCV angles -> alpha beta
      print*,'TCV angles convention: transforming to (alpha,beta)...'
      cvdr=pi/180.0_r8
      polrad=anglp*cvdr
      torrad=anglt*cvdr
      anglp=atan2(sin(polrad)*cos(torrad),cos(polrad))/cvdr
      anglt=asin(-sin(polrad)*sin(torrad))/cvdr
      write(*,"(' theta,phi -> alpha,beta=',f8.4,1x,f8.4)") anglp,anglt
    end if
  end select
  y00=0.0d0
  cvdr=pi/180.0d0
  anglt=anglt*cvdr
  anglp=anglp*cvdr
! anglt=phi, anglp=-theta
!  anx0c=cos(anglt)*cos(anglp)
!  any0c=-sin(anglt)*cos(anglp)
!  anz0c=sin(anglp)
! anglt=beta, anglp=alpha
  anx0c=cos(anglt)*cos(anglp)
  any0c=-sin(anglt)
  anz0c=cos(anglt)*sin(anglp)
end subroutine read_param

subroutine read_data
  use kinds, only: r8
  implicit none
  character(len=255) idstring,prefix
  character(len=255) filenmeqq,filenmte,filenmne
  integer neqdsk,nprofn,nproft,nfitn
  logical alloutput,goutput
  common/out/alloutput,goutput
  common/nfile/neqdsk,nprofn,nproft,nfitn
  common/spectinfo/idstring,prefix
  common/eqlfnm/filenmeqq,filenmte,filenmne

  interface
    subroutine profdata
    end subroutine profdata
    subroutine equidata
    end subroutine equidata
    subroutine cutoff
    end subroutine cutoff
  end interface

! unit numbers for density, temperature, density fit and eqdsk files
  neqdsk=99
  nprofn=98
  nproft=97
  nfitn=96
  open(nprofn,file=trim(filenmne),status='old',action='read')
  open(nproft,file=trim(filenmte),status='old',action='read')
  if(alloutput) open(nfitn,file=trim(prefix)//'ne'//trim(idstring)//'.fit.prf',action='write')
  call profdata
  close(nprofn)
  close(nproft)
  if(alloutput) close(nfitn)
  open(neqdsk,file=trim(filenmeqq),status='old',action='read')
  call equidata
  close(neqdsk)
  if(alloutput) call cutoff
end subroutine read_data

subroutine equidata
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! EQUIDATA reads equilibrium from EQDSK file                           C
!          computes coefficients for psi and fpol spline interpolation C
!   uses                                                               C
!   sets                                                               C
!   calls DENSITY                                                      C
!         POL_LIMIT                                                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi
  use spece_data, only : nnw,nnh,nbb,nhlmx,naxmx
  use utils, only: inside
  implicit none
  integer nu
  integer, parameter :: kspl = 3
  integer, parameter :: nrest = nnw+4, nzest = nnh+4,                    &
         lwrk = 4*(nnw+nnh) + 11*(nrest+nzest) + nrest*nnh + nzest + 54, &
         liwrk = nnh + nnw + nrest + nzest + 3
  integer, parameter :: lwrkbsp = 4*(nnw+nnh), liwrkbsp = nnw+nnh
  integer, parameter :: lw10 = nnw*3 + nnh*4 + nnw*nnh, &
                        lw01 = nnw*4 + nnh*3 + nnw*nnh, &
                        lw20 = nnw*2 + nnh*4 + nnw*nnh, &
                        lw02 = nnw*4 + nnh*2 + nnw*nnh, &
                        lw11 = nnw*3 + nnh*3 + nnw*nnh, &
                        lwrkf = nnw*4 + nrest*16,       &
                        ldiwrk = nnw + nnh
  real(kind=r8)  fpol(nnw),pres(nnw),qpsi(nnw)
  real(kind=r8)  ffprim(nnw),pprim(nnw)
  real(kind=r8)  psi(nnw,nnh),rv(nnw),zv(nnh),psin(nnw,nnh),psinr(nnw)
  real(kind=r8)  fvpsi(nnw*nnh),cc(nnw*nnh),ffvpsi(nnw*nnh)
  real(kind=r8), allocatable :: rv1d(:),zv1d(:),wf1d(:)
  real(kind=r8)  tr(nrest),tz(nzest),wrk(lwrk)
  integer iwrk(liwrk)
  real(kind=r8)  cc01(lw01),cc10(lw10),cc02(lw02),cc20(lw20),cc11(lw11)
  integer iwrkd(ldiwrk)
  real(kind=r8)  tfp(nrest),cfp(nrest),wrkf(lwrkf),wf(nnw)
  real(kind=r8)  tffp(nrest),cffp(nrest),wrkff(lwrkf),wrkfd(lwrkf)
  integer iwrkf(nrest)
  real(kind=r8)  wrkbsp(lwrkbsp)
  integer iwrkbsp(liwrkbsp)
  real(kind=r8)  fpoli(nnw),ffpol(nnw),dfpolv(nnw)
  integer neqdsk,nprofn,nproft,nfitn,nr,nz,nppn,nppt,i,j,ij,nrz
  integer nur,nuz,nsrt,nszt,nsft,nsftt,ier,iopt,nsr,nsz,ixp
  real(kind=r8) drnr1,dznz1,rcen,rleft,zmid,rmaxis,zmaxis,psiax,psiedge
  real(kind=r8) btrcen,current,simag,xdum,sibry,bfact,bvfact,fpolas,psia,psiant
  real(kind=r8) fp,xb,xe,sspl,ssfp
  real(kind=r8) psia0,btaxis
  real(kind=r8) dr,dz,dpsinr,rmnm,rmxm,zmnm,zmxm
  integer nbbbs,limitr
  real(kind=r8) rbmin,zbmin,rbmax,zbmax,delpsinox,psin1 !,deltapsi
  real(kind=r8) rbbbs(nbb),zbbbs(nbb)
  integer info
  integer nhl,nax,iproj
  real(kind=r8) hlr(nhlmx),hlz(nhlmx),chl(nhlmx,4)
  real(kind=r8) hlri,hlzi,hlpsii,dzhl
  real(kind=r8) psinv,psieq(naxmx),req(naxmx)
  real(kind=r8) drax,dst,psinoptmp,psinxptmp,psinxp2tmp
  real(kind=r8) r10,z10,rxp,zxp,r1,rn,z1,rop,zop,psinxp,psinop,rxp2,zxp2,psinxp2
  character(len=50) strdum
  integer idum
  integer ipsinorm,icocos,icocosmod,idsc,iffmt,iscatter
  real(kind=r8)  sgnbphi,sgniphi
  real(kind=r8), dimension(nnw,4) :: cq

  common/pareq1/psia
  common/pareq1t/psiant,psinop
  common/pareq2/rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz
  common/baxis/btaxis
  common/bound/zbmin,zbmax,rbmin,rbmax
  common/psinr/psinr
  common/cpsin/rv,zv,psin
  common/cpsi/psi
  common/eqnn/nr,nz,nppn,nppt
  common/nfile/neqdsk,nprofn,nproft,nfitn
  common/ixp/ixp
  common/eqconv/ipsinorm,icocos,sgnbphi,sgniphi
  common/ssplps/sspl
  common/coffeqt/tr,tz,tfp,tffp
  common/coffeq1/cc,cc01,cc10
  common/coffeq2/cc20,cc02,cc11
  common/coffeqn/nsrt,nszt,nsft,nsftt
  common/cofffp/cfp,cffp
  common/fpas/fpolas
  common/bfact/bfact,bvfact
  common/hline/chl,hlr,hlz,nhl
  common/psival/psinv
  common/axis/psieq,req,nax
  common/dsds/dst
  common/eqdskfmt/idsc,iffmt,iscatter
  common/iproj/iproj
  common/coffq/cq
  
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
    subroutine points_ox(rz,zz,rf,zf,psinvf,info)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: rz, zz
      real(kind=r8), intent(out) :: rf, zf, psinvf
      integer      , intent(out) :: info
    end subroutine points_ox
    subroutine points_tgo(rz,zz,rf,zf,psin,info)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: rz,zz,psin
      real(kind=r8), intent(out) :: rf,zf
      integer      , intent(out) :: info
    end subroutine points_tgo
    subroutine points_mino(rz,zz,rf,psinf,info)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: rz, zz
      real(kind=r8), intent(out) :: rf, psinf
      integer      , intent(out) :: info
    end subroutine points_mino
    subroutine swap(a,b)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(inout) :: a, b
    end subroutine swap
    subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: nx, ny, kx, ky, mx, my, lwrk, kwrk
      integer, intent(out)   :: ier
      integer, intent(inout) :: iwrk(kwrk)
      real(kind=r8), intent(in)    :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
      real(kind=r8), intent(out)   :: z(mx*my)
      real(kind=r8), intent(inout) :: x(mx), y(my), wrk(lwrk)
    end subroutine bispev
    subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,&
         wrk,lwrk,iwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: iopt, m, k, nest, lwrk
      integer, intent(out) :: ier
      integer, intent(inout) :: n, iwrk(nest)
      real(kind=r8), intent(in)    :: xb, xe, s, x(m), y(m), w(m)
      real(kind=r8), intent(out)   :: fp, c(nest)
      real(kind=r8), intent(inout) :: t(nest), wrk(lwrk)
    end subroutine curfit
    subroutine difcs(xx,yy,iopt,cc,ier)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: iopt
      real(kind=r8), intent(in)    :: xx(:),yy(:)
      real(kind=r8), intent(inout) :: cc(:,:)
      integer      , intent(out)   :: ier
    end subroutine difcs
    subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,&
       wrk,lwrk,iwrk,kwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: nx, ny, kx, ky, nux, nuy, mx, my, lwrk, kwrk
      integer, intent(out)   :: ier
      integer, intent(inout) :: iwrk(kwrk)
      real(kind=r8), intent(in)    :: tx(nx), ty(ny), c((nx-kx-1)*(ny-ky-1)), &
                                      x(mx), y(my)
      real(kind=r8), intent(out)   :: z(mx*my)
      real(kind=r8), intent(inout) :: wrk(lwrk)
    end subroutine parder
    subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,&
       nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: iopt, mx, my, kx, ky, nxest, nyest, lwrk, kwrk
      integer, intent(out)   :: ier
      integer, intent(inout) :: nx, ny, iwrk(kwrk)
      real(kind=r8), intent(in)    :: xb, xe, yb, ye, x(mx), y(my), z(mx*my), s
      real(kind=r8), intent(out)   :: fp, c((nxest-kx-1)*(nyest-ky-1))
      real(kind=r8), intent(inout) :: tx(nxest), ty(nyest), wrk(lwrk)
    end subroutine regrid
    subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, nu, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)    :: t(n), c(n), x(m)
      real(kind=r8), intent(out)   :: y(m), wrk(n)
    end subroutine splder
    subroutine splev(t,n,c,k,x,y,m,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)  :: t(n), c(n), x(m)
      real(kind=r8), intent(out) :: y(m)
    end subroutine splev
    subroutine rhotor
      use kinds, only : r8
      use spece_data, only : nnw
      implicit none
    end subroutine rhotor
  end interface

!     read from file eqdsk
!     (see http://fusion.gat.com/efit/g_eqdsk.html)
!     spline smoothing of psi and derivatives
  if(idsc==1) then
    read (neqdsk,2000) strdum,idum,nr,nz
  else
    read (neqdsk,*) nr,nz
  end if
  if(iffmt==0) then
    read (neqdsk,2020) drnr1,dznz1,rcen,rleft,zmid
    read (neqdsk,2020) rmaxis,zmaxis,psiax,psiedge,btrcen
    read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
    read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
    read (neqdsk,2020) (fpol(i),i=1,nr)
    read (neqdsk,2020) (pres(i),i=1,nr)
    read (neqdsk,2020) (ffprim(i),i=1,nr)
    read (neqdsk,2020) (pprim(i),i=1,nr)
    read (neqdsk,2020) ((psi(i,j),i=1,nr),j=1,nz)
    read (neqdsk,2020) (qpsi(i),i=1,nr)
  else
    read (neqdsk,*) drnr1,dznz1,rcen,rleft,zmid
    read (neqdsk,*) rmaxis,zmaxis,psiax,psiedge,btrcen
    read (neqdsk,*) current,simag,xdum,rmaxis,xdum
    read (neqdsk,*) zmaxis,xdum,sibry,xdum,xdum
    read (neqdsk,*) (fpol(i),i=1,nr)
    read (neqdsk,*) (pres(i),i=1,nr)
    read (neqdsk,*) (ffprim(i),i=1,nr)
    read (neqdsk,*) (pprim(i),i=1,nr)
    read (neqdsk,*) ((psi(i,j),i=1,nr),j=1,nz)
    read (neqdsk,*) (qpsi(i),i=1,nr)
  end if
!     read plasma boundaries from eqdsk file
  read (neqdsk,*) nbbbs,limitr
  if(nbbbs>0) then
    if(iffmt==0) then
        read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
    else
        read (neqdsk,*) (rbbbs(i),zbbbs(i),i=1,nbbbs)
    end if
  end if
  2000 format (a48,3i4)
  2020 format (5e16.9)
!  remove psi normalization
  if (ipsinorm/=0) then 
    do j=1,nz
      do i=1,nr
        psi(i,j)=psi(i,j)*(psiedge-psiax)+psiax
      end do
    end do
  end if
!  compensate for different reference systems
  icocosmod=mod(icocos,10)
  if (mod(icocos,2)==0.and.icocosmod>0) then
!  icocos mod 2 = 0: toroidal angle phi CW (opposite to gray convention)
    btrcen=-btrcen
    current=-current
    do i=1,nr
      fpol(i)=-fpol(i)
    end do
  end if
!
  if (icocosmod==1 .or. icocosmod==4 .or. &
      icocosmod==5 .or. icocosmod==8) then
!  icocos mod 10 = 1,4,5,8: psi increasing with CCW Ip
!  icocos mod 10 = 2,3,6,7: psi decreasing with CCW Ip
    psiedge=-psiedge
    psiax=-psiax
    do j=1,nz
      do i=1,nr
        psi(i,j)=-psi(i,j)
      end do
    end do
  end if
! add check for Ip/psi and B0/Fpol sign consistency?
  current=sign(current,psiax-psiedge)
  btrcen=sign(btrcen,fpol(nr))

!     length in m !!!
  dr=drnr1/dble(nr-1)
  dz=dznz1/dble(nz-1)
  rv(1)=rleft
  zv(1)=zmid-dznz1/2.0d0
  dpsinr=1.0d0/dble(nr-1)
  do i=1,nr
    psinr(i)=(i-1)*dpsinr
    rv(i)=rv(1)+(i-1)*dr
  end do
  do j=1,nz
    zv(j)=zv(1)+(j-1)*dz
  end do
  rmnm=rv(1)
  rmxm=rv(nr)
  zmnm=zv(1)
  zmxm=zv(nz)
!  psi function
  psia0=psiedge-psiax
!  icocos=0: adapt psi to force Ip sign, otherwise maintain psi
  if (icocosmod/=0) sgniphi=sign(1.0_r8,-psia0)
  current=sign(current,sgniphi)
  psia=-sgniphi*abs(psia0)*bfact
!  icocos>10: input psi is in Wb
!  icocos<10: input psi is in Wb/rad (gray convention)
  if (icocos>=10) psia=psia/(2.0d0*pi)
  do j=1,nz
    do i=1,nr
      psin(i,j)=(psi(i,j)-psiax)/psia0
      psi(i,j)=psin(i,j)*psia
    enddo
  enddo

  if (iscatter.eq.1) then
! valid psi(R,z) data only inside plasma boundary
! determine number of valid grid points
    nrz=0
    do j=1,nz
      do i=1,nr
        if (nbbbs.gt.0) then
          if(.not.inside(rbbbs,zbbbs,nbbbs,rv(i),zv(j))) cycle
        else
          if(psin(i,j).le.0.0d0) cycle
        end if
        nrz=nrz+1
      end do
    end do
! store valid data
    allocate(rv1d(nrz),zv1d(nrz),wf1d(nrz))
    ij=0
    do j=1,nz
      do i=1,nr
        if (nbbbs.gt.0) then
          if(.not.inside(rbbbs,zbbbs,nbbbs,rv(i),zv(j))) cycle
        else
          if(psin(i,j).le.0.0d0) cycle
        end if
        ij=ij+1
        rv1d(ij)=rv(i)
        zv1d(ij)=zv(j)
        fvpsi(ij)=psin(i,j)
        wf1d(ij)=1.0d0
      end do
    end do
!     fit as a scattered set of points
!     use reduced number of knots to limit memory comsumption ?
    nsr=nr/4+4
    nsz=nz/4+4
    sspl=sspl*nrz*1.0d-5
    call scatterspl(rv1d,zv1d,fvpsi,wf1d,nrz,kspl,sspl,               &
                    rmnm,rmxm,zmnm,zmxm,tr,nsr,tz,nsz,cc,ier)
!     if ier=-1 data are fitted using sspl=0    
    if(ier.eq.-1) then
      sspl=0.0d0
      nsr=nr/4+4
      nsz=nz/4+4
      call scatterspl(rv1d,zv1d,fvpsi,wf1d,nrz,kspl,sspl,             &
                      rmnm,rmxm,zmnm,zmxm,tr,nsr,tz,nsz,cc,ier)
    end if
    deallocate(rv1d,zv1d,wf1d)
  else
! valid psi(R,z) data on the whole grid
    do j=1,nz
      do i=1,nr
        ij=nz*(i-1)+j
        fvpsi(ij)=psin(i,j)
      enddo
    enddo
!     spline interpolation of psi(i,j) and derivatives
    sspl=sspl*nr*nz*1.0d-5
    iopt=0
    call regrid(iopt,nr,rv,nz,zv,fvpsi,rmnm,rmxm,zmnm,zmxm,           &
                  kspl,kspl,sspl,nrest,nzest,nsr,tr,nsz,tz,cc,fp,     &
                  wrk,lwrk,iwrk,liwrk,ier)
! if ier=-1 data are fitted using sspl=0    
    if(ier.eq.-1) then
      sspl=0.0d0
      call regrid(iopt,nr,rv,nz,zv,fvpsi,rmnm,rmxm,zmnm,zmxm,         &
                  kspl,kspl,sspl,nrest,nzest,nsr,tr,nsz,tz,cc,fp,     &
                  wrk,lwrk,iwrk,liwrk,ier)
    end if
  end if
  nsrt=nsr
  nszt=nsz
  if (sspl>0.0d0) then
    call bispev(tr,nsr,tz,nsz,cc,kspl,kspl,rv,nr,zv,nz,ffvpsi,      &
        wrkbsp,lwrkbsp,iwrkbsp,liwrkbsp,ier)
    do j=1,nz
      do i=1,nr
        ij=nz*(i-1)+j
        psin(i,j)=ffvpsi(ij)
        psi(i,j)=psin(i,j)*psia
      enddo
    enddo
  end if
  nur=1
  nuz=0
  call parder(tr,nsr,tz,nsz,cc,kspl,kspl,nur,nuz,rv,nr,zv,nz,       &
   ffvpsi,cc10,lw10,iwrkd,ldiwrk,ier)
  nur=0
  nuz=1
  call parder(tr,nsr,tz,nsz,cc,kspl,kspl,nur,nuz,rv,nr,zv,nz,       &
   ffvpsi,cc01,lw01,iwrkd,ldiwrk,ier)
  nur=2
  nuz=0
  call parder(tr,nsr,tz,nsz,cc,kspl,kspl,nur,nuz,rv,nr,zv,nz,       &
   ffvpsi,cc20,lw20,iwrkd,ldiwrk,ier)
  nur=0
  nuz=2
  call parder(tr,nsr,tz,nsz,cc,kspl,kspl,nur,nuz,rv,nr,zv,nz,       &
   ffvpsi,cc02,lw02,iwrkd,ldiwrk,ier)
  nur=1
  nuz=1
  call parder(tr,nsr,tz,nsz,cc,kspl,kspl,nur,nuz,rv,nr,zv,nz,       &
   ffvpsi,cc11,lw11,iwrkd,ldiwrk,ier)
! scaling of f_poloidal
! icocos=0: adapt fpol to force Ip sign, otherwise maintain fpol
  if (icocosmod/=0) sgnbphi=sign(1.0_r8,fpol(nr))
    btrcen=sign(btrcen,sgnbphi)
  do i=1,nr
    qpsi(i)=qpsi(i)*sqrt(1.0d0+(bvfact**2-1.d0)*(rcen*btrcen/fpol(i))**2)
    fpol(i)=sgnbphi*sign(1.d0,fpol(nr))*bfact*sign(1.d0,fpol(i))*   &
              sqrt(fpol(i)**2+(bvfact**2-1.d0)*(rcen*btrcen)**2)
    pres(i)=pres(i)*bfact*bfact
    wf(i)=1.0d0
  end do
!  wf(nr)=1.0d2
! spline interpolation of fpol(i)
  iopt=0
  xb=0.0d0
  xe=1.0d0
  ssfp=0.0d0
!  ssfp=0.01d0
  call curfit(iopt,nr,psinr,fpol,wf,xb,xe,kspl,ssfp,nrest,nsft,     &
   tfp,cfp,fp,wrkf,lwrkf,iwrkf,ier)
  call splev(tfp,nsft,cfp,3,psinr,fpoli,nr,ier)
  fpolas=fpoli(nr)
  nu=1
  call splder(tfp,nsft,cfp,3,nu,psinr,dfpolv,nr,wrkfd,ier)
  do i=1,nr
    ffpol(i)=fpol(i)*dfpolv(i)/psia
    wf(i)=1.0d0
  end do
  iopt=0
  xb=0.0d0
  xe=1.0d0
  call curfit(iopt,nr,psinr,ffpol,wf,xb,xe,kspl,ssfp,nrest,nsftt,   &
   tffp,cffp,fp,wrkff,lwrkf,iwrkf,ier)
!   if(ixp/=0) then 
!   compute max and min z of last closed surface 
  rbmin=rmaxis
  rbmax=rmaxis
  zbmin=1.0d+30
  zbmax=-1.0d+30
  if (nbbbs>1) then
    do i=1,nbbbs
      if(zbbbs(i)<zbmin) then
        zbmin=zbbbs(i)
        rbmin=rbbbs(i)
      end if
      if(zbbbs(i)>zbmax) then
        zbmax=zbbbs(i)
        rbmax=rbbbs(i)
      end if
    end do 
  else
    zbmin=(2.0d0*zmnm+zmaxis)/3.0d0
    rbmin=rmaxis
    zbmax=(2.0d0*zmxm+zmaxis)/3.0d0
    rbmax=rmaxis
  end if
  if(zbmin<=zmnm) zbmin=zmnm+dz
  if(rbmin<=rmnm) rbmin=rmnm+dr
  if(zbmax>=zmxm) zbmax=zmxm-dz
  if(rbmax>=rmxm) rbmax=rmxm-dr
!   end if
!     start with uncorrected normalized psi
  psinop=0.0d0
  psinxp=1.0d0
  psiant=1.0d0
!     search for O-point
  call points_ox(rmaxis,zmaxis,rop,zop,psinoptmp,info)
  rmaxis=rop
  zmaxis=zop
  write (*,"('O point: (',f7.3,',',f7.3,') psi=',e11.4)") rop,zop,psinoptmp
!     compute B_toroidal on axis
  btaxis=fpol(1)/rop
!     search for X-point
  if (ixp/=0) then
    psinxptmp=5.0d0
    if (ixp/=1) then
      ! search lower X
      r10=rbmin
      z10=zbmin
      call points_ox(r10,z10,rxp,zxp,psinxptmp,info)
    end if
    psinxp2tmp=5.0d0
    if (ixp/=-1) then
      ! search upper X
      r10=rop
      z10=zbmax
      call points_ox(r10,z10,rxp2,zxp2,psinxp2tmp,info)
    end if
    if (psinxptmp<5.0d0) then
      if (psinxp2tmp<5.0d0) then
        ! 2 X points found. check which is innermost.
        if (psinxp2tmp<psinxptmp) then
          call swap(psinxp2tmp,psinxptmp)
          call swap(rxp2,rxp)
          call swap(zxp2,zxp)
          ixp=1
        else
          ixp=-1
        end if
      else
        ! only lower X
        ixp=-1
      end if
    else if (psinxp2tmp<5.0d0) then
      ! only upper X
      psinxptmp=psinxp2tmp
      rxp=rxp2
      zxp=zxp2
      ixp=1
      psinxp2tmp=5.0d0
    else
      ! no X point
      ixp=0
    end if
  end if
!     psi re-normalization
  if (ixp==0) then
    psin1=1.0d0
    psinop=psinoptmp
    psiant=psin1-psinop
    r10=rmaxis
    z10=(zbmax+zmaxis)/2.0d0
    call points_tgo(r10,z10,r1,z1,psin1,info)
    rbmax=r1
    zbmax=z1
    z10=(zbmin+zmaxis)/2.0d0
    call points_tgo(r10,z10,r1,z1,psin1,info)
    zbmin=z1
    rbmin=r1
    write(*,"('no X-point  ',4f8.4)") rbmin,zbmin,rbmax,zbmax
  else
    psinop=psinoptmp
    psinxp=psinxptmp
    psiant=psinxp-psinop
    psin1=1.0d0
    r10=rmaxis
    if (ixp<0) then
      rbmin=rxp
      zbmin=zxp
      z10=(zbmax+zmaxis)/2.0d0
      call points_tgo(r10,z10,r1,z1,psin1,info)
      rbmax=r1
      zbmax=z1
    else
      rbmax=rxp
      zbmax=zxp
      z10=(zbmin+zmaxis)/2.0d0
      call points_tgo(r10,z10,r1,z1,psin1,info)
      rbmin=r1
      zbmin=z1
    end if
    write (*,"('X point: (',f7.3,',',f7.3,') psi=',e11.4)") rxp,zxp,psinxp
    if (psinxp2tmp<5.0d0) then
      psinxp2=psinxp2tmp
      write (*,"('2nd Xpt: (',f7.3,',',f7.3,') psi=',e11.4)") rxp2,zxp2,psinxp2
    end if
  end if  
!   compute LFS/HFS boundary
  nhl=int((zbmax-zbmin)/dz)+1
  if (nhl>nhlmx) then
    nhl=nhlmx
  else if (nhl<2) then
    nhl=2
  end if
  dzhl=(zbmax-zbmin)/dble(nhl-1)
  do i=1,nhl
    hlzi=zbmin+dble(i-1)*dzhl
    call points_mino(rmaxis,hlzi,hlri,hlpsii,info)
    hlr(i)=hlri
    hlz(i)=hlzi
  end do
  iopt=0
  call difcs(hlz(:nhl),hlr(:nhl),iopt,chl(:nhl,:),ier)
!   compute psi on equatorial plane
  r1=rmxm-dr
  rn=rmnm+dr
  drax=dst*1.d-2
  nax=int((r1-rn)/drax)+1
  if (nax>naxmx) then
    nax=naxmx
  else if (nax<2) then
    nax=2
  end if
  if(iproj.eq.1) then
    drax=(rn-r1)/dble(nax-1)
    do i=1,nax
      req(i)=r1+dble(i-1)*drax
      call equinum(req(i),zop)
      if (req(i)>=rop) then
        psieq(i)=psinv
      else
        psieq(i)=-psinv
      end if
    end do
  else
    drax=(rop-r1)/dble(nax-1)
    do i=1,nax
      req(i)=r1+dble(i-1)*drax
      call equinum(req(i),zop)
      psieq(i)=psinv
    end do
  end if
! 
! q interpolation
  iopt=0
  call difcs(psinr(1:nr),qpsi(1:nr),iopt,cq,ier)
! rhot interpolation
  call rhotor
end subroutine equidata

subroutine profdata
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PROFDATA reads temperature and density plofiles from file            C
!          computes coefficients for ne and Te interpolation           C
!          polynomial transition to zero density at plasma boundary    C
!          linear transition to zero temperature at plasma boundary    C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : zero,one,tiny_r8
  use cubic, only : qdcrt, cbcrt
  use spece_data, only : npmx
  implicit none
  integer, parameter :: npest = npmx+4, lwrkf = npmx*4 + npest*16
  real(kind=r8) psradn(npmx),psradt(npmx),terad(npmx),derad(npmx)
  real(kind=r8) ct(npmx,4)
  real(kind=r8) tfn(npest),cfn(npest),wrkf(lwrkf)
  integer iwrkf(npest)
  real(kind=r8) wfn(npmx)
  real(kind=r8) densi(npest),ddensi(1),d2densi(1),d3densi(1),wrkfd(npest)
  real(kind=r8) psdbnd,psnpp,sspld,tfact,nfact
  real(kind=r8) xe,xb,fp
  real(kind=r8) dnpp,ddnpp,d2dnpp,d3dnpp,dpsb
  real(kind=r8) a(0:3), psi_solr(3)
  real(kind=r8) dens_sol, eta, eps, psiin, psiout
  real(kind=r8) sol_p, sol_m, delta
  complex(kind=r8) psi_solc(3)
  integer npsi_sol
  integer nppn,nppt,neqdsk,nprofn,nproft,nfitn,i,iopt,ier,kspl,n_int,isol
  integer nr,nz,nn,nn1,nn2,nu,nsfdt
  character(len=10) strdum
  real(kind=r8) psntpp
  real(kind=r8) delta2,xnv,xxm,xxp
!    DENSITY PROFILE DEBUG
  real(kind=r8) psistep,psii,dpsimax,dpsii
  logical alloutput,goutput,sharp_dbnd

  common/out/alloutput,goutput
  common/densbnd/psdbnd
  common/denspp/psnpp
  common/cdenspp/a
  common/temppp/psntpp
  common/eqnn/nr,nz,nppn,nppt
  common/nfile/neqdsk,nprofn,nproft,nfitn
  common/crad/psradt,terad
  common/coffte/ct
  common/coffdt/tfn
  common/coffdnst/nsfdt
  common/cofffn/cfn
  common/sspldn/sspld
  common/tfact/tfact
  common/nfact/nfact
  common/dbndtype/sharp_dbnd

  interface
    subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,&
         wrk,lwrk,iwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: iopt, m, k, nest, lwrk
      integer, intent(out) :: ier
      integer, intent(inout) :: n, iwrk(nest)
      real(kind=r8), intent(in)    :: xb, xe, s, x(m), y(m), w(m)
      real(kind=r8), intent(out)   :: fp, c(nest)
      real(kind=r8), intent(inout) :: t(nest), wrk(lwrk)
    end subroutine curfit
    subroutine difcs(xx,yy,iopt,cc,ier)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: iopt
      real(kind=r8), intent(in)    :: xx(:),yy(:)
      real(kind=r8), intent(inout) :: cc(:,:)
      integer      , intent(out)   :: ier
    end subroutine difcs
    subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, nu, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)    :: t(n), c(n), x(m)
      real(kind=r8), intent(out)   :: y(m), wrk(n)
    end subroutine splder
    subroutine splev(t,n,c,k,x,y,m,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)  :: t(n), c(n), x(m)
      real(kind=r8), intent(out) :: y(m)
    end subroutine splev
    function denspoly(psi)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: denspoly
      real(kind=r8), intent(in) :: psi
    end function denspoly
  end interface

!     read plasma profiles from files
  read(nprofn,*) nppn
  read(nprofn,*) strdum
  do i=1,nppn
    read(nprofn,*) psradn(i),derad(i)
    derad(i)=derad(i)*nfact
    wfn(i)=one
  end do
  read(nproft,*) nppt
  read(nproft,*) strdum
  do i=1,nppt
    read(nproft,*) psradt(i),terad(i)
    terad(i)=terad(i)*tfact
  end do
  if(psradn(1)/=zero) psradn(1)=zero
  if(psradt(1)/=zero) psradt(1)=zero
!     spline approximation of density
  iopt=0
  xb=zero
  xe=psradn(nppn)
  kspl=3
  call curfit(iopt,nppn,psradn,derad,wfn,xb,xe,kspl,sspld,npest,    &
      nsfdt,tfn,cfn,fp,wrkf,lwrkf,iwrkf,ier)
  call splev(tfn,nsfdt,cfn,3,psradn,densi,nppn,ier)
  dnpp=densi(nppn)
  call splder(tfn,nsfdt,cfn,3,1,psradn(nppn:nppn),ddensi(1:1),1,wrkfd,ier)
  ddnpp=ddensi(1)
  call splder(tfn,nsfdt,cfn,3,2,psradn(nppn:nppn),d2densi(1:1),1,wrkfd,ier)
  d2dnpp=d2densi(1)
  call splder(tfn,nsfdt,cfn,3,3,psradn(nppn:nppn),d3densi(1:1),1,wrkfd,ier)
  d3dnpp=d3densi(1)
  psnpp=psradn(nppn)
!     density transition to 0
  if(dnpp<tiny_r8 .or. psdbnd<zero) then
    psdbnd=psnpp
    sharp_dbnd=.true.
  else
    sharp_dbnd=.false.
  end if
  if (sharp_dbnd) then
    a=(/dnpp,ddnpp,d2dnpp/2._r8,d3dnpp/6._r8/)
    call cbcrt(a, psi_solc, npsi_sol, psi_solr)
    isol=0
    do i=1,npsi_sol
      if (psi_solr(i)>zero) then
        isol = i
        exit
      end if
    end do
    if (isol>0) then
      psdbnd = psnpp + psi_solr(isol)
    else
      do i=1,npsi_sol
        if (psnpp+psi_solr(i)>psradn(nppn-1)) then
          isol = i
          exit
        end if
      end do
      if (isol>0) then
        psdbnd = psnpp + psi_solr(isol)
        nppn=nppn-1
        dnpp=densi(nppn)
        psnpp=psradn(nppn)
      else
        write(*,*) 'Bad ne spline extrapolation. Try with polynomial extension'
        sharp_dbnd=.false.
      end if
    end if
    if (sharp_dbnd) then
      dens_sol = denspoly(psdbnd)
      eta = sqrt(tiny_r8)
      eps = 1e-12_r8
      if (abs(dens_sol)>eta) then
        if (dens_sol<zero) then
          psiin = psnpp - 1e-2_r8
          psiout = psdbnd
        else if (dens_sol>zero) then
          psiin = psdbnd
          psiout = psdbnd + 1e-2_r8
        end if
        call zeroin(psiin,psiout,eps,eta,denspoly,psdbnd,ier)
      end if
    end if
  end if
  if (.not.sharp_dbnd) then
    a(0)=dnpp
    a(1)=ddnpp
    a(2)=d2dnpp
    delta2=(a(1)/a(2))**2_r8-2.0_r8*a(0)/a(2)
    xnv=psnpp-a(1)/a(2)
    if(psdbnd<=psnpp) psdbnd=psnpp*1.01_r8
    
    if(delta2 < zero) then
      if(xnv > psnpp) psdbnd=min(psdbnd,xnv)
    else
      xxm=xnv-sqrt(delta2)
      xxp=xnv+sqrt(delta2)
      if(xxm > psnpp) then
        psdbnd=min(psdbnd,xxm)
      else if (xxp > psnpp) then
        psdbnd=min(psdbnd,xxp)
      end if    
    end if
  end if
  write(*,"(' Density boundary at psi =',f11.8)") psdbnd
!  interpolated density output
  if(alloutput) then
    psistep=1.d-2
    n_int=nint((psdbnd-psnpp)/psistep)
    n_int=max(n_int,0)
    write(nfitn,*) nppn+n_int
    write(nfitn,*) ' #psin ne_fit'
    do i=1,nppn
      write(nfitn,*) psradn(i),densi(i)
    enddo
    if (n_int>0) then
      if (sharp_dbnd) then
        psistep=(psdbnd-psnpp)/n_int
        do i=1,n_int
          psii=i*psistep
          write(nfitn,*) psnpp+psii,a(0)+a(1)*psii+a(2)*psii**2+a(3)*psii**3
        enddo
      else
        dpsii=zero
        dpsimax=psdbnd-psnpp
        do i=1,n_int
          dpsii=dpsii+psistep
          write(nfitn,*) psnpp+dpsii, (a(0) + a(1)*dpsii + 0.5_r8*a(2)*dpsii**2) *                &
                                      (one - dpsii/dpsimax)**3_r8 *                               &
                                      (one + 3.0_r8*dpsii/dpsimax + 6.0_r8*(dpsii/dpsimax)**2_r8)
        end do
      end if
    end if
  end if
!     temperature check at boundary
  do while (psradt(nppt)>psdbnd.and.nppt>1)
    nppt=nppt-1
  end do
  if (terad(nppt)>zero) then
    if (psradt(nppt)<psdbnd) then
      nppt=nppt+1
      psradt(nppt)=psdbnd
    end if
    terad(nppt)=zero
  end if
  psntpp=psradt(nppt)
  iopt=0
  call difcs(psradt(:nppt-1),terad(:nppt-1),iopt,ct(:nppt-1,:),ier)
end subroutine profdata

function denspoly(psi)
  use kinds, only : r8
  implicit none
  real(kind=r8) :: denspoly
  real(kind=r8), intent(in) :: psi
  real(kind=r8) ad(0:3),psnpp,dpsib
  common/cdenspp/ad
  common/denspp/psnpp
  dpsib=psi-psnpp
  denspoly = ad(0) + dpsib*(ad(1) + dpsib*(ad(2) + dpsib*ad(3)))
end function denspoly

subroutine taildata(filenmtail)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TAILDATA                                                             C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use phys_constants, only : me_kev
  use spece_data, only : ncmx
  implicit none
  character(len=255), intent(in) :: filenmtail
  integer, parameter :: uedf=3

  integer neqdsk,nprofn,nproft,nfitn
  integer icmp,ncmp0,fddtype
  real(kind=r8) u,te,et0,ps0,psc
  real(kind=r8) tekev(ncmx),terad
  real(kind=r8) udr(ncmx),gam(ncmx)
  real(kind=r8) urel(ncmx),trel(ncmx)
  real(kind=r8) eta0(ncmx),psi0(ncmx),psic(ncmx)
  real(kind=r8) eta(ncmx),amu(ncmx)
  real(kind=r8) umax,deltau

  common/nfile/neqdsk,nprofn,nproft,nfitn
  common/etaprf/eta0,psi0,psic
  common/lhpar/umax,deltau
  common/fddpar/urel,trel
  common/tete/tekev,terad
  common/udrf/udr,gam
  common/multc/eta,amu
  common/ncmp0/ncmp0
  common/fddtyp/fddtype

  udr(1)=0.0d0
  gam(1)=1.0d0
  if (fddtype<0.or.fddtype>3) fddtype=0
  if (fddtype==0) then
    write(*,*) 'Maxwellian e.d.f.'
    ncmp0=1
  else
    open(uedf,file=trim(adjustl(filenmtail)),status='old',action='read')
    if (fddtype==1) then              ! LH tail model
      write(*,*) 'LH tail model e.d.f.'
      read(uedf,*) umax, te, et0, ps0, psc
      do icmp=2,ncmx
        tekev(icmp)=te  ! constant parameters for all components
        eta0(icmp)=et0  ! except for u_drift
        psi0(icmp)=ps0  ! (calculated in subroutine fdd_cmposition)
        psic(icmp)=psc
        amu(icmp)=me_kev/te
      enddo
      if (umax>=0.0d0) then
        deltau=dsqrt(4.0d0*te/me_kev)  ! du between components
      else
        deltau=-dsqrt(4.0d0*te/me_kev)
      endif
    else                        ! manual fdd composition
      write(*,*) 'Generic multi-Maxwellian e.d.f.'
      read(uedf,*) ncmp0           ! number of tail components
      ncmp0=ncmp0+1             ! total components (incl. bulk)
      if (ncmp0<2) then
        ncmp0=2
      elseif (ncmp0>ncmx) then
        ncmp0=ncmx
        write(*,*) 'ncmp = ',ncmp0,' > ncmx'
      endif
      do icmp=2,ncmp0           ! read tail components parameters
        read(uedf,*) u, te, eta0(icmp), psi0(icmp), psic(icmp)
        if (fddtype==2) then  ! absolute parameters
          udr(icmp)=u
          gam(icmp)=dsqrt(1.0d0+udr(icmp)**2)
          tekev(icmp)=te
          amu(icmp)=me_kev/te
        else                    ! relative parameters
          urel(icmp)=u          ! (scaled to bulk temperature)
          trel(icmp)=te         ! abs params will be calculated in
        endif                   ! subroutine fdd_composition
      enddo
    endif
    close(uedf)
  endif
end subroutine taildata

subroutine fdd_composition(psinv)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! FDD_COMPOSITION(psi)                                                 C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : huge_r8
  use phys_constants, only : me_kev
  use spece_data, only : ncmx
  implicit none
  real(kind=r8), intent(in) :: psinv
  real(kind=r8), parameter  :: etacr = 1.0d-6
  integer icmp,ncmp,ncmp0,fddtype
  real(kind=r8) tekev(ncmx),terad,eta(ncmx),amu(ncmx)
  real(kind=r8) urel(ncmx),trel(ncmx)
  real(kind=r8) u1,uthbulk,udr(ncmx),gam(ncmx)
  real(kind=r8) umax,deltau
  real(kind=r8) etasum
  common/tete/tekev,terad
  common/ndrf/ncmp
  common/ncmp0/ncmp0
  common/udrf/udr,gam
  common/lhpar/umax,deltau
  common/fddpar/urel,trel
  common/multc/eta,amu
  common/fddtyp/fddtype

  interface
    function temp(arg)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: temp
      real(kind=r8), intent(in)  :: arg
    end function temp
    function etafunc(icmp,psi)
      use kinds, only : r8
      implicit none
      real(kind=r8) ::etafunc
      integer, intent(in) :: icmp
      real(kind=r8), intent(in) :: psi
    end function etafunc
  end interface

  tekev(1)=temp(psinv)
  if (tekev(1)>0.0d0) then
    amu(1)=me_kev/tekev(1)
    if (fddtype==0) then          ! 0= maxwellian fdd: eta=0
      etasum=0.0d0
    elseif (fddtype==1) then  ! LH model: eta(2)=eta of whole tail
      etasum=etafunc(2,psinv)
    else                        ! manual mode: sum over components
      etasum=0.0d0
      do icmp=2,ncmp0
        eta(icmp)=etafunc(icmp,psinv)
        etasum=etasum+eta(icmp)
      enddo
    endif
    if (etasum>etacr) then
      if (fddtype==1) then        ! 1= LH model
        if (umax>=0.0d0) then
          u1=1.0d0/dsqrt(amu(1))
        else
          u1=-1.0d0/dsqrt(amu(1))
        endif
        ncmp=2+nint((umax-u1)/deltau)
        if (ncmp>1) then
          if (ncmp>ncmx) then
            ncmp=ncmx
            write(*,*) 'LH tail incomplete: ncmp = ',ncmp,' > ncmx'
          endif
          eta(1)=1.0d0-etasum
          etasum=etasum/dble(ncmp-1)
          do icmp=2,ncmp
            udr(icmp)=u1+(icmp-2)*deltau
            gam(icmp)=dsqrt(1.0d0+udr(icmp)**2)
            eta(icmp)=etasum
          enddo
        else
          ncmp=1
          eta(1)=1.0d0
        endif
      else                          ! 2,3= multi-maxwellian
        eta(1)=1.0d0-etasum
        ncmp=ncmp0
        if (fddtype==3) then      ! 3= parameters scaled to T_bulk
          uthbulk=1.d0/dsqrt(amu(1))              !** CHECK U_TH **
          do icmp=2,ncmp                          !** DEFINITION **
            udr(icmp)=urel(icmp)*uthbulk
            gam(icmp)=dsqrt(1.0d0+udr(icmp)**2)
            tekev(icmp)=trel(icmp)*tekev(1)
            amu(icmp)=me_kev/tekev(icmp)
          enddo
        endif
      endif
    else         ! only maxwellian bulk if deformation is negligible
      ncmp=1
      eta(1)=1.0d0
    endif
  else           ! only bulk with T=0
    amu(1)=huge_r8
    eta(1)=1.0d0
    ncmp=1
  endif
end subroutine fdd_composition

function etafunc(icmp,psi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ETAFUNC(icmp,psi)                                                    C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : ncmx
  implicit none
! arguments
  real(kind=r8) ::etafunc
  integer, intent(in) :: icmp
  real(kind=r8), intent(in) :: psi
! common variables
  real(kind=r8) eta0(ncmx),psi0(ncmx),psic(ncmx)
  common/etaprf/eta0,psi0,psic
! executable statements
  etafunc=eta0(icmp)*exp( -((psi-psi0(icmp))/psic(icmp))**2 )
end function etafunc

subroutine print_current
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PRINT_CURRENT                                                        C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8  
  use constants, only : zero, one, pi
  use phys_constants, only : vc_si, qe_si, mu0
  use spece_data, only : ncuts, n_psi, ncmx
  implicit none
  integer, parameter :: ncon = 2*ncuts + 1
  real(kind=r8) :: rcon(ncon), zcon(ncon)
  real(kind=r8) :: dpsi, psi, psi0, currdens, dlprec, dlnext, dl, rnow, znow
  real(kind=r8) :: bphi, br, bz, bpol, btot, coeff_j, curr_edf, curr_pla
  real(kind=r8) :: dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  real(kind=r8) :: btaxis, dens, ddens, jtor_edf0, jtor_edf
  real(kind=r8) :: psia
  integer :: i_psi, icmp, icon, inext, ncmp
  real(kind=r8) :: eta(ncmx), amu(ncmx)
  real(kind=r8) :: udr(ncmx), gam(ncmx)
  logical, parameter :: cont_prn = .true.

  common/ndrf/ncmp
  common/udrf/udr,gam
  common/multc/eta,amu
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/pareq1/psia
  common/baxis/btaxis
  common/dens/dens,ddens

  character(len=255) :: idstring,prefix
  common/spectinfo/idstring,prefix
  
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
    subroutine contours_psi(psin,ncuts,rcnt,zcnt,cont_prn)
      use kinds, only : r8  
      implicit none
      integer, intent(in) :: ncuts
      real(kind=r8), intent(in) :: psin
      real(kind=r8), intent(out) :: rcnt(2*ncuts+1), zcnt(2*ncuts+1)
      logical, intent(in) :: cont_prn
    end subroutine contours_psi
  end interface

  if (cont_prn) then
    open(71,file=trim(prefix)//'contours'//trim(idstring)//'.txt',action='write')
    write(71,*) n_psi
    write(71,*)' #i psi R z'
  end if
  open(72,file=trim(prefix)//'current'//trim(idstring)//'.txt',action='write')
  write(72,*) n_psi+1
  write(72,*)' #psi cfjpl j//edf I_edf I_plasma'
  coeff_j = 2.0d0*pi/btaxis
  curr_edf = zero
  jtor_edf = zero
  psi = zero
  do i_psi=0,n_psi
    jtor_edf0 = jtor_edf
    psi0 = psi
    psi = (i_psi*(0.9999d0/n_psi))**2
    dpsi = psi-psi0
    call fdd_composition(psi)
    currdens = zero
    do icmp=1,ncmp
      currdens = currdens + eta(icmp)*udr(icmp)/gam(icmp)
    end do
    call density(psi)
    currdens = -currdens*dens*1.d19*qe_si*vc_si
    if (psi <= zero) then
      coeff_j = 2.0d0*pi/btaxis
      curr_pla = zero
    else
      call contours_psi(psi,ncuts,rcon,zcon,cont_prn)
      dlprec = sqrt((rcon(ncon)-rcon(ncon-1))**2+(zcon(ncon)-zcon(ncon-1))**2)
      coeff_j = zero
      curr_pla = zero
      do icon=1,ncon-1
        rnow = rcon(icon)
        znow = zcon(icon)
        inext = icon + 1
        dlnext = sqrt((rcon(inext)-rnow)**2+(zcon(inext)-znow)**2)
        dl = (dlprec + dlnext)/2
        call equinum(rnow,znow)
        bphi = fpolv/rnow
        br = -dpsidz/rnow
        bz =  dpsidr/rnow
        btot = br**2+bz**2
        bpol = sqrt(btot)
        btot = sqrt(btot+bphi**2)
        coeff_j = coeff_j + dl/(bpol*rnow)*bphi/btot
        curr_pla = curr_pla + bpol*dl
        dlprec = dlnext
      end do
      curr_pla = sign(1.0d0,psia)*curr_pla/mu0
    end if
    jtor_edf = currdens*coeff_j
    curr_edf = curr_edf + 0.5d0*(jtor_edf0+jtor_edf)*dpsi
    write(72,'(12(1x,e12.5))') psi, coeff_j, currdens, curr_edf, curr_pla
  end do
  if (cont_prn) close(71)
  close(72)
end subroutine print_current

subroutine contours_psi(psin,ncuts,rcnt,zcnt,cont_prn)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CONTOURS_PSI                                                         C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi
  use spece_data, only : nnw, nnh
  implicit none
  integer, intent(in) :: ncuts
  real(kind=r8), intent(in) :: psin
  real(kind=r8), intent(out) :: rcnt(2*ncuts+1), zcnt(2*ncuts+1)
  logical, intent(in) :: cont_prn
  integer, parameter :: kspl=3
  integer, parameter :: mest=4, nrest=nnw+4, nzest=nnh+4
  integer, parameter :: lw10 = nnw*3 + nnh*4 + nnw*nnh, lw01 = nnw*4 + nnh*3 + nnw*nnh
  real(kind=r8) :: cc(nnw*nnh), tr(nrest), tz(nzest), tfp(nrest), tffp(nrest)
  real(kind=r8) :: cc01(lw01),cc10(lw10)
  real(kind=r8) :: czc(nrest), rzero(mest)
  real(kind=r8) :: psia,psiant,psinop
  real(kind=r8) :: rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz,zbmin,zbmax,rbmin,rbmax
  real(kind=r8) :: r0, z0, rup, zup, rlw, zlw, z_hsum, z_hdif, theta, zcut
  real(kind=r8) :: rhopol,value
  integer :: ic, iopt, info, m, ii, ier
  integer ::is(2*ncuts+1),ms(2*ncuts+1)
  integer :: nsr, nsz, nsft, nsftt

  common/pareq1/psia
  common/pareq1t/psiant,psinop
  common/pareq2/rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz
  common/bound/zbmin,zbmax,rbmin,rbmax
  common/coffeqn/nsr,nsz,nsft,nsftt
  common/coffeq1/cc,cc01,cc10
  common/coffeqt/tr,tz,tfp,tffp

  interface
    subroutine profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
      use kinds, only : r8
      implicit none
      integer,intent(in) :: iopt,nx,ny,kx,ky,nu
      integer,intent(out) :: ier
      real(kind=r8), intent(in) :: u
      real(kind=r8),intent(in) :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
      real(kind=r8),intent(out) :: cu(nu)
    end subroutine profil
    subroutine sproota(val,t,n,c,zero,mest,m,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n,mest
      integer, intent(out) :: m,ier
      real(kind=r8), intent(in) :: val
      real(kind=r8), intent(in) :: t(n),c(n)
      real(kind=r8), intent(out) :: zero(mest)
    end subroutine sproota
  end interface

  rhopol=sqrt(psin)
!  r0=(rbmax+rmaxis)/2.0d0
!  z0=(zbmax+zmaxis)/2.0d0
  r0=rhopol*(rbmax-rmaxis)+rmaxis
  z0=rhopol*(zbmax-zmaxis)+zmaxis
  call points_tgo(r0,z0,rup,zup,psin,info)
!  r0=(rmaxis+rbmin)/2.0d0
!  z0=(zmaxis+zbmin)/2.0d0
  r0=rhopol*(rbmin-rmaxis)+rmaxis
  z0=rhopol*(zbmin-zmaxis)+zmaxis
  call points_tgo(r0,z0,rlw,zlw,psin,info)  

  theta=pi/dble(ncuts)
  z_hsum=(zup+zlw)/2.0d0
  z_hdif=(zup-zlw)/2.0d0
  rcnt(1)=rlw
  zcnt(1)=zlw
  rcnt(2*ncuts+1)=rlw
  zcnt(2*ncuts+1)=zlw
  rcnt(ncuts+1)=rup
  zcnt(ncuts+1)=zup
  is(1)=1
  ms(1)=1
  is(2*ncuts+1)=1
  ms(2*ncuts+1)=1
  is(ncuts+1)=1
  ms(ncuts+1)=1

  do ic=2,ncuts
    zcut=z_hsum-z_hdif*cos(theta*(ic-1))
    iopt=1
    call profil(iopt,tr,nsr,tz,nsz,cc,kspl,kspl,zcut,nsr,czc,ier)
    if(ier>0) print*,'    profil =',ier
    value=psin*psiant+psinop
    call sproota(value,tr,nsr,czc,rzero,mest,m,ier)
    if (m>1) then
      ii=1
      do while (m-ii>1)
!        if (rzero(ii+1)<rlw .and. rzero(ii+1)<rup) then
        if ((rzero(ii+1)-rcnt(ic-1))**2 + &
            (rzero(ii+2)-rcnt(2*ncuts+2-ic+1))**2 < &
            (rzero(ii)-rcnt(ic-1))**2 + &
            (rzero(ii+1)-rcnt(2*ncuts+2-ic+1))**2 ) then
          ii=ii+1
        else
          exit
        end if
      end do
      rcnt(ic)=rzero(ii)
      zcnt(ic)=zcut
      rcnt(2*ncuts+2-ic)=rzero(ii+1)
      zcnt(2*ncuts+2-ic)=zcut
      is(ic)=ii
      ms(ic)=m
      is(2*ncuts+2-ic)=ii+1
      ms(2*ncuts+2-ic)=m
    else
      rcnt(ic)=rcnt(ic-1)
      zcnt(ic)=zcnt(ic-1)
      rcnt(2*ncuts+2-ic)=rcnt(2*ncuts+2-ic+1)
      zcnt(2*ncuts+2-ic)=zcnt(2*ncuts+2-ic+1)
      is(ic)=0
      ms(ic)=0
      is(2*ncuts+2-ic)=0
      ms(2*ncuts+2-ic)=0
    end if
  end do

  if (cont_prn) then
    do ii=1,2*ncuts+1
      write(71,'(i6,3(1x,e12.5),2i6)') ii,psin,rcnt(ii),zcnt(ii),is(ii),ms(ii)
    end do
    write(71,*) ' '
  end if
end subroutine contours_psi

subroutine datinit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DATINIT                                                              C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx, npts, tmax
  implicit none
  integer i,nstep,iwarm,ilarm
  integer iovlp,iharm,nharm,iif,icall,ipolc
  real(kind=r8) ddrmx,ttar1,ttar2
  real(kind=r8), parameter :: dt = 2.0d0*tmax/npts
  real(kind=r8) psjki(nmx)
  real(kind=r8) tauv(nmx),alphav(nmx)
  real(kind=r8) tev(nmx),tta(nmx),rrr(nmx),ygv(nmx),anplv(nmx)
  real(kind=r8) ttv(npts+1),extv(npts+1)
  real(kind=r8) zzv(nmx)
  real(kind=r8) stv(nmx),phiv(nmx),xgv(nmx),anprv(nmx)


  common/harm/iovlp,iharm
  common/ttar/ttar1,ttar2
  common/nharm/nharm
  common/ddrmx/ddrmx
  common/iiv/iif
  common/icall/icall
  common/ttex/ttv,extv
  common/psjki/psjki
  common/atjki/tauv,alphav,tev
  common/tta/tta
  common/init/nstep
  common/ygv/ygv
  common/warm/iwarm,ilarm
  common/nplv/anplv
  common/traj/stv,zzv,phiv,xgv,anprv
  common/rrr/rrr

  common/ipol/ipolc

  ttar1=0.0d0
  ttar2=0.0d0
  iovlp=0
  iharm=99
  ipolc=0
  nharm=0
  iif=1
  icall=-2
  ddrmx=0.0d0
  do i=1,nstep
    psjki(i)=0.0d0
    tauv(i)=0.0d0
    alphav(i)=0.0d0
    tev(i)=0.0d0
    tta(i)=0.0d0
    ygv(i)=0.0d0
    anplv(i)=0.0d0
    stv(i)=0.0d0
    zzv(i)=0.0d0
    phiv(i)=0.0d0
    xgv(i)=0.0d0
    anprv(i)=0.0d0
  end do
  if(iwarm>1) then
    
    do i = 1, npts+1
      ttv(i) = -tmax+dble(i-1)*dt
      extv(i)=exp(-ttv(i)*ttv(i))
    end do
  end if
end subroutine datinit

subroutine rkint4
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! RKINT4                                                               C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  integer :: ieq
  integer, parameter :: ndim = 6
  real(kind=r8)  yy(ndim)
  real(kind=r8)  fk1(ndim),fk2(ndim),fk3(ndim),fk4(ndim)
! common variables
  real(kind=r8)  ywrk(ndim),ypwrk(ndim)
  real(kind=r8)  dst,h,hh,h6
  common/dsds/dst
  common/wrk/ywrk,ypwrk
! external subroutines
  interface
    subroutine fwork(y,dery)
      use kinds, only : r8
      implicit none
      integer, parameter :: ndim = 6
      real(kind=r8), intent(in)  :: y(ndim)
      real(kind=r8), intent(out) :: dery(ndim)
    end subroutine fwork
  end interface
! executable statements
!  queste assegnazioni potrebbero essere fatte una volta sola
      h=-dst
  hh=h*0.5d0
  h6=h/6.0d0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  do ieq=1,ndim
    fk1(ieq)=ypwrk(ieq)
    yy(ieq)=ywrk(ieq)+fk1(ieq)*hh
  end do
  call fwork(yy,fk2)
  do ieq=1,ndim
    yy(ieq)=ywrk(ieq)+fk2(ieq)*hh
  end do
  call fwork(yy,fk3)
  do ieq=1,ndim
    yy(ieq)=ywrk(ieq)+fk3(ieq)*h
  end do
  call fwork(yy,fk4)
  do ieq=1,ndim
    ywrk(ieq)=ywrk(ieq)+h6*(fk1(ieq)+2.0d0*fk2(ieq)+2.0d0*fk3(ieq)+fk4(ieq))
  end do
  call fwork(ywrk,ypwrk)
end subroutine rkint4

subroutine fwork(y,dery)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! FWORK                                                                C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  integer, parameter :: ndim = 6
  real(kind=r8), intent(in)  :: y(ndim)
  real(kind=r8), intent(out) :: dery(ndim)
  real(kind=r8) anv(3),vgv(3),bv(3),derbv(3,3),derxg(3),deryg(3)
  real(kind=r8) derdxv(3),danpldxv(3),derdnv(3)
  integer ierr,iv
  real(kind=r8) xx,yy,zz,rr,anx,any,anz,an2,anpl,anpl2,dnl,anpr,anpr2,an2s
  real(kind=r8) xg,yg,yg2,duh,del,sox,vgm,denom,derdnm,derdom,dd
  real(kind=r8) dan2sdxg,dan2sdyg,dan2sdnpl

  common/ddd/dd,an2s,an2
  common/mode/sox
  common/nplr/anpl,anpl2,anpr
  common/resah/dnl
  common/bb/bv
  common/dbb/derbv
  common/xgxg/xg
  common/ygyg/yg
  common/dxgyg/derxg,deryg
  common/vgv/vgm,derdnm
  common/ierr/ierr
  common/rr/rr
  common/anv/anx,any,anz

  interface
    subroutine plas_deriv(xx,yy,zz)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: xx, yy, zz
    end subroutine plas_deriv
  end interface

  xx=y(1)
  yy=y(2)
  zz=y(3)
  anx=y(4)
  any=y(5)
  anz=y(6)
  rr=sqrt(xx**2+yy**2)
  anv(1) = anx
  anv(2) = any
  anv(3) = anz
  call plas_deriv(xx,yy,zz)
  an2  = 0.0d0
  anpl = 0.0d0
  do iv=1,3
    an2=an2+anv(iv)*anv(iv)
    anpl=anpl+anv(iv)*bv(iv)
  end do
  if(abs(anpl)>=1) then
    print*,'      |Nparallel| > 1 !'
    ierr=98
  end if
  anpl2=anpl*anpl
  dnl=1.0d0-anpl2
  anpr2=an2-anpl2
  anpr=0.0d0
  if(anpr2>0.0d0) anpr=sqrt(anpr2)
  yg2=yg**2
!  if(xg>0.0d0) then
    duh=1.0d0-xg-yg2
    del=sqrt(dnl*dnl+4.0d0*anpl*anpl*(1.0d0-xg)/yg2)
    an2s=1.0d0-xg -xg*yg2*(1.0d0+anpl2+sox*del)/duh/2.0d0
    dan2sdxg=-yg2*(1.0d0-yg2)*(1.0d0+anpl2+sox*del)/duh**2/2.0d0    &
             +sox*xg*anpl2/(del*duh)-1.0d0
    dan2sdyg=-xg*yg*(1.0d0-xg)*(1.0d0+anpl2+sox*del)/duh**2         &
             +2.0d0*sox*xg*(1.0d0-xg)*anpl2/(yg*del*duh)
    dan2sdnpl=-xg*yg2*anpl/duh                                      &
              -sox*xg*anpl*(2.0d0*(1.0d0-xg)-yg2*dnl)/(del*duh)
    derdnm=0.0d0
    do iv=1,3
      danpldxv(iv)=anv(1)*derbv(1,iv)+anv(2)*derbv(2,iv)+anv(3)*derbv(3,iv)
      derdxv(iv)=-(derxg(iv)*dan2sdxg+deryg(iv)*dan2sdyg+danpldxv(iv)*dan2sdnpl)
      derdnv(iv)=2.0d0*anv(iv)-bv(iv)*dan2sdnpl
      derdnm=derdnm+derdnv(iv)**2
    end do
    derdnm=sqrt(derdnm)
    derdom=-2.0d0*an2+2.0d0*xg*dan2sdxg+yg*dan2sdyg+anpl*dan2sdnpl
!  else
!    del=sqrt(dnl*dnl+4.0d0*anpl*anpl/yg2)
!    duh=1.0d0-yg2
!    dan2sdxg=-yg2*(1.0d0+anpl2+sox*del)/duh/2.0d0
!    derdnm=0.0d0
!    do iv=1,3
!      derdxv(iv)=-derxg(iv)*dan2sdxg
!      derdnv(iv)=2.0d0*anv(iv)
!      derdnm=derdnm+derdnv(iv)**2
!    end do
!    derdnm=sqrt(derdnm)
!    derdom=-2.0d0*an2
!  end if
  denom=-derdom
  dery(1) = derdnv(1)/denom
  dery(2) = derdnv(2)/denom
  dery(3) = derdnv(3)/denom
  dery(4) = -derdxv(1)/denom
  dery(5) = -derdxv(2)/denom
  dery(6) = -derdxv(3)/denom
!     vgv :  group velocity
!  vgm=0
!  do iv=1,3
!    vgv(iv)=-derdnv(iv)/derdom
!    vgm=vgm+vgv(iv)**2
!    vgm=vgm+derdnv(iv)**2
!  end do
!  vgm=sqrt(vgm)
!  vgm=-sqrt(vgm)/derdom
  vgm=-derdnm/derdom
!     dd :  dispersion relation
  dd=an2-an2s
!  print*,an2,an2s,duh
end subroutine fwork

subroutine plas_deriv(xx,yy,zz)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PLAS_DERIV(x,y,z)                                                    C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  real(kind=r8), intent(in) :: xx, yy, zz
  real(kind=r8) bv(3),derxg(3),deryg(3),derbv(3,3),dbtot(3)
  real(kind=r8) bvc(3),dbvcdc(3,3),dbvdc(3,3),dbv(3,3)
  integer iv,jv
  real(kind=r8) xg,yg,psinv,rr,csphi,snphi,zzm,rrm
  real(kind=r8) bphi,brr,bzz,btot,b2tot,bres
  real(kind=r8) dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv,dxgdpsi
  real(kind=r8) derxgdr,derxgdz,dphidx,dphidy,drrdx,drrdy,dfpolv

  common/parbres/bres
  common/btot/btot
  common/bb/bv
  common/dbb/derbv
  common/xgxg/xg
  common/dxgdps/dxgdpsi
  common/ygyg/yg
  common/dxgyg/derxg,deryg
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/psival/psinv
! external subroutines
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
    subroutine sub_xg
    end subroutine sub_xg
  end interface
! executable statements
  xg=0.0d0
  yg=9.9d1
  do iv=1,3
    derxg(iv)=0.0d0
    deryg(iv)=0.0d0
    bv(iv)=0.0d0
    dbtot(iv)=0.0d0
    do jv=1,3
      dbv(iv,jv)=0.0d0
      derbv(iv,jv)=0.0d0
      dbvcdc(iv,jv)=0.0d0
      dbvcdc(iv,jv)=0.0d0
      dbvdc(iv,jv)=0.0d0
    end do
  end do
! cylindrical coordinates
  rr=sqrt(xx**2+yy**2)
  csphi=xx/rr
  snphi=yy/rr
  bv(1)=snphi
  bv(2)=-csphi
! convert from cm to meters
  zzm=zz/1.0d2
  rrm=rr/1.0d2
  call equinum(rrm,zzm)
  yg=abs(fpolv)/(rrm*bres)
  bphi=fpolv/rrm
  btot=abs(bphi)
!  if(psinv<0.0d0) return
  if(psinv>4.9d0) return
!     B = f(psi)/R e_phi+  grad psi x e_phi/R
!     bvc(i) = B_i in cylindrical coordinates
!     bv(i)  = B_i in cartesian coordinates
  bphi=fpolv/rrm
  brr=-dpsidz/rrm
  bzz= dpsidr/rrm
  dfpolv=ffpv/fpolv
  bvc(1)=brr
  bvc(2)=bphi
  bvc(3)=bzz
  bv(1)=bvc(1)*csphi-bvc(2)*snphi
  bv(2)=bvc(1)*snphi+bvc(2)*csphi
  bv(3)=bvc(3)
  b2tot=bv(1)**2+bv(2)**2+bv(3)**2
  btot=sqrt(b2tot)
!     dbvcdc(iv,jv) = d Bcil(iv) / dxvcil(jv)
  dbvcdc(1,1)=-ddpsidrz/rrm-brr/rrm
  dbvcdc(1,3)=-ddpsidzz/rrm
  dbvcdc(2,1)= dfpolv*dpsidr/rrm-bphi/rrm
  dbvcdc(2,3)= dfpolv*dpsidz/rrm
  dbvcdc(3,1)= ddpsidrr/rrm-bzz/rrm
  dbvcdc(3,3)= ddpsidrz/rrm
!     dbvdc(iv,jv) = d Bcart(iv) / dxvcil(jv)
  dbvdc(1,1) = dbvcdc(1,1)*csphi-dbvcdc(2,1)*snphi
  dbvdc(1,2) = -bv(2)
  dbvdc(1,3) = dbvcdc(1,3)*csphi-dbvcdc(2,3)*snphi
  dbvdc(2,1) = dbvcdc(1,1)*snphi+dbvcdc(2,1)*csphi
  dbvdc(2,2) = bv(1)
  dbvdc(2,3) = dbvcdc(1,3)*snphi+dbvcdc(2,3)*csphi
  dbvdc(3,1) = dbvcdc(3,1)
  dbvdc(3,2) = dbvcdc(3,2)
  dbvdc(3,3) = dbvcdc(3,3)
  drrdx=csphi
  drrdy=snphi
  dphidx=-snphi/rrm
  dphidy=csphi/rrm
!     dbv(iv,jv) = d Bcart(iv) / dxvcart(jv)
  dbv(1,1)=drrdx*dbvdc(1,1)+dphidx*dbvdc(1,2)
  dbv(1,2)=drrdy*dbvdc(1,1)+dphidy*dbvdc(1,2)
  dbv(1,3)=dbvdc(1,3)
  dbv(2,1)=drrdx*dbvdc(2,1)+dphidx*dbvdc(2,2)
  dbv(2,2)=drrdy*dbvdc(2,1)+dphidy*dbvdc(2,2)
  dbv(2,3)=dbvdc(2,3)
  dbv(3,1)=drrdx*dbvdc(3,1)+dphidx*dbvdc(3,2)
  dbv(3,2)=drrdy*dbvdc(3,1)+dphidy*dbvdc(3,2)
  dbv(3,3)=dbvdc(3,3)
  dbtot(1)=(bv(1)*dbv(1,1)+bv(2)*dbv(2,1)+bv(3)*dbv(3,1))/btot
  dbtot(2)=(bv(1)*dbv(1,2)+bv(2)*dbv(2,2)+bv(3)*dbv(3,2))/btot
  dbtot(3)=(bv(1)*dbv(1,3)+bv(2)*dbv(2,3)+bv(3)*dbv(3,3))/btot
!     bv(i) = B_i / B ; derbv(i,j) = d (B_i / B) /d x,y,z
  do iv=1,3
    bv(iv)=bv(iv)/btot
  end do
  do iv=1,3
    do jv=1,3
      derbv(iv,jv)=(dbv(iv,jv)-bv(iv)*dbtot(jv))/btot
    end do
  end do
  yg=btot/Bres
  do iv=1,3
    deryg(iv)=dbtot(iv)/Bres
  end do
! derivatives of X
  call sub_xg
  derxgdr=dpsidr*dxgdpsi
  derxgdz=dpsidz*dxgdpsi
  derxg(1)=drrdx*derxgdr
  derxg(2)=drrdy*derxgdr
  derxg(3)=derxgdz
! convert spatial derivatives from dummy/m -> dummy/cm
! to be used in fwork
  do iv=1,3
    deryg(iv)=1.0d-2*deryg(iv)
    derxg(iv)=1.0d-2*derxg(iv)
    do jv=1,3
      derbv(iv,jv)=1.0d-2*derbv(iv,jv)
    end do
  end do
end subroutine plas_deriv

subroutine equinum(rpsim,zpsim)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! EQUINUM(r,z)                                                         C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nnw, nnh
  implicit none
  real(kind=r8), intent(in) :: rpsim, zpsim
  integer, parameter :: nrest = nnw+4, nzest = nnh+4, &
                       lwrk = 8, liwrk = 2, lw11 = nnw*3 + nnh*3 + nnw*nnh, &
            lw10 = nnw*3 + nnh*4 + nnw*nnh, lw01 = nnw*4 + nnh*3 + nnw*nnh, &
            lw20 = nnw*2 + nnh*4 + nnw*nnh, lw02 = nnw*4 + nnh*2 + nnw*nnh
  real(kind=r8) ccspl(nnw*nnh)
  real(kind=r8) tr(nrest),tz(nzest),wrk(lwrk)
  integer iwrk(liwrk)
  real(kind=r8) rrs(1),zzs(1),ffspl(1)
  real(kind=r8) cc01(lw01),cc10(lw10),cc02(lw02),cc20(lw20),cc11(lw11)
  real(kind=r8) tfp(nrest),cfp(nrest)
  real(kind=r8) tffp(nrest),cffp(nrest)
  integer nr,nz,nppn,nppt,nsrt,nszt,nsft,nsftt,ier
  integer nrs,nzs,nsr,nsz,nur,nuz,kkr,kkz,iwr,iwz
  real(kind=r8) psinv,psia,psiant,psinop
  real(kind=r8) dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv,fpolas
  real(kind=r8) rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz

  common/eqnn/nr,nz,nppn,nppt
  common/psival/psinv
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/pareq1/psia
  common/pareq1t/psiant,psinop
  common/pareq2/rmaxis,zmaxis,rmnm,rmxm,zmnm,zmxm,dr,dz
  common/coffeq1/ccspl,cc01,cc10
  common/coffeq2/cc20,cc02,cc11
  common/coffeqt/tr,tz,tfp,tffp
  common/coffeqn/nsrt,nszt,nsft,nsftt
  common/cofffp/cfp,cffp
  common/fpas/fpolas

  interface
    subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: nx, ny, kx, ky, mx, my
      integer, intent(out) :: lx(mx), ly(my)
      real(kind=r8), intent(in)  :: tx(nx), ty(ny), c((nx-kx-1)*(ny-ky-1)), &
          x(mx), y(my)
      real(kind=r8), intent(out) :: z(mx*my)
      real(kind=r8), intent(out) :: wx(mx,kx+1), wy(my,ky+1)
    end subroutine fpbisp
    subroutine splev(t,n,c,k,x,y,m,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)  :: t(n), c(n), x(m)
      real(kind=r8), intent(out) :: y(m)
    end subroutine splev
  end interface

!      psinv=-1.0d0
  psinv=5.0d0
  dpsidr=0.0d0
  dpsidz=0.0d0
  ddpsidrr=0.0d0
  ddpsidzz=0.0d0
  ddpsidrz=0.0d0
!     here lengths are measured in meters
  fpolv=fpolas
  ffpv=0.0d0
  if (rpsim>=rmxm.or.rpsim<=rmnm) return
  if (zpsim>=zmxm.or.zpsim<=zmnm) return
  rrs(1)=rpsim
  zzs(1)=zpsim
  nrs=1
  nzs=1
  nsr=nsrt
  nsz=nszt
  call fpbisp(tr,nsr,tz,nsz,ccspl,3,3, &
   rrs,nrs,zzs,nzs,ffspl,wrk(1),wrk(5),iwrk(1),iwrk(2))
  psinv=(ffspl(1)-psinop)/psiant
  nur=1
  nuz=0
  kkr=3-nur
  kkz=3-nuz
  iwr=1+(nr-nur-4)*(nz-nuz-4)
  iwz=iwr+4-nur
  nsr=nsrt
  nsz=nszt
  call fpbisp(tr(nur+1),nsr-2*nur,tz(nuz+1),nsz-2*nuz,cc10,kkr,kkz, &
   rrs,nrs,zzs,nzs,ffspl,cc10(iwr),cc10(iwz),iwrk(1),iwrk(2))
  dpsidr= ffspl(1)*psia
  nur=0
  nuz=1
  kkr=3-nur
  kkz=3-nuz
  iwr=1+(nr-nur-4)*(nz-nuz-4)
  iwz=iwr+4-nur
  call fpbisp(tr(nur+1),nsr-2*nur,tz(nuz+1),nsz-2*nuz,cc01,kkr,kkz, &
   rrs,nrs,zzs,nzs,ffspl,cc01(iwr),cc01(iwz),iwrk(1),iwrk(2))
  dpsidz= ffspl(1)*psia
  nur=2
  nuz=0
  kkr=3-nur
  kkz=3-nuz
  iwr=1+(nr-nur-4)*(nz-nuz-4)
  iwz=iwr+4-nur
  call fpbisp(tr(nur+1),nsr-2*nur,tz(nuz+1),nsz-2*nuz,cc20,kkr,kkz, &
   rrs,nrs,zzs,nzs,ffspl,cc20(iwr),cc20(iwz),iwrk(1),iwrk(2))
  ddpsidrr= ffspl(1)*psia
  nur=0
  nuz=2
  kkr=3-nur
  kkz=3-nuz
  iwr=1+(nr-nur-4)*(nz-nuz-4)
  iwz=iwr+4-nur
  call fpbisp(tr(nur+1),nsr-2*nur,tz(nuz+1),nsz-2*nuz,cc02,kkr,kkz, &
   rrs,nrs,zzs,nzs,ffspl,cc02(iwr),cc02(iwz),iwrk(1),iwrk(2))
  ddpsidzz= ffspl(1)*psia
  nur=1
  nuz=1
  kkr=3-nur
  kkz=3-nuz
  iwr=1+(nr-nur-4)*(nz-nuz-4)
  iwz=iwr+4-nur
  call fpbisp(tr(nur+1),nsr-2*nur,tz(nuz+1),nsz-2*nuz,cc11,kkr,kkz, &
   rrs,nrs,zzs,nzs,ffspl,cc11(iwr),cc11(iwz),iwrk(1),iwrk(2))
  ddpsidrz= ffspl(1)*psia
  if(psinv<=1.0d0) then
    rrs(1)=psinv
    call splev(tfp,nsft,cfp,3,rrs,ffspl,1,ier)
    fpolv=ffspl(1)
    call splev(tffp,nsftt,cffp,3,rrs,ffspl,1,ier)
    ffpv=ffspl(1)
  end if
end subroutine equinum

subroutine sub_xg
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUB_XG                                                               C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use phys_constants, only : cxg19
  implicit none
  real(kind=r8) xgcn19,psinv,psia,psdbnd,xg,dxgdpsi,dens,ddenspsin,fghz

  common/psival/psinv
  common/pareq1/psia
  common/densbnd/psdbnd
  common/xgxg/xg
  common/dxgdps/dxgdpsi
  common/dens/dens,ddenspsin
  common/fghz/fghz

  interface
    subroutine density(arg)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: arg
    end subroutine density
  end interface

!  xg=0.0d0
!  dxgdpsi=0.0d0
!  if(psinv<=psdbnd) then
    call density(psinv)
    xgcn19=cxg19/fghz**2
    xg=xgcn19*dens
    dxgdpsi=xgcn19*ddenspsin/psia
!  end if
end subroutine sub_xg

subroutine density(arg)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DENSITY(arg)                                                         C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : zero,one
  use spece_data, only : npmx
  implicit none
  real(kind=r8), intent(in) :: arg
  integer, parameter :: npest = npmx+4
  real(kind=r8) xxs(1),ffs(1)
  real(kind=r8) tfn(npest),cfn(npest),wrkfd(npest)
  real(kind=r8) psdbnd,psnpp,aad,bbd,ccd,dens,ddens,dpsib
  real(kind=r8) fp,fh,dfp,dfh,tt
  real(kind=r8) ad(0:3)
  integer nn,nn1,nn2,ier,nu,nsft
  logical sharp_dbnd

  common/densbnd/psdbnd
  common/denspp/psnpp
  common/cdenspp/ad
  common/dens/dens,ddens
  common/coffdt/tfn
  common/coffdnst/nsft
  common/cofffn/cfn
  common/dbndtype/sharp_dbnd

  interface
    subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, nu, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)    :: t(n), c(n), x(m)
      real(kind=r8), intent(out)   :: y(m), wrk(n)
    end subroutine splder
    subroutine splev(t,n,c,k,x,y,m,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, k, m
      integer, intent(out) :: ier
      real(kind=r8), intent(in)  :: t(n), c(n), x(m)
      real(kind=r8), intent(out) :: y(m)
    end subroutine splev
  end interface

!     computation of density [10^19 m^-3] and derivative wrt psi
  dens=zero
  ddens=zero
  if(arg>psdbnd.or.arg<zero) return
  if(arg>psnpp) then
!     cubic interpolation for 1 < psi < psdbnd
    if (sharp_dbnd) then
      dpsib=arg-psnpp
      dens = ad(0) + dpsib*(      ad(1) + dpsib*(ad(2) + dpsib*ad(3)))
      ddens= ad(1) + dpsib*(2._r8*ad(2) + dpsib*3._r8*ad(3))
    else
      dpsib=arg-psnpp
      tt=dpsib/(psdbnd-psnpp)
      fp = ad(0)  + ad(1) *dpsib + 0.5_r8*ad(2) *dpsib**2.0_r8
      dfp=          ad(1)        +        ad(2) *dpsib
      fh = (one-tt)**3.0_r8 * (one + 3.0_r8*tt + 6.0_r8*tt**2_r8)
      dfh= -30.0_r8 * (one-tt)**2.0_r8 * tt**2.0_r8 / (psdbnd-psnpp)
      dens=fp*fh
      ddens=dfp*fh+fp*dfh
    end if
  else
    xxs(1)=arg
    ier=0
    call splev(tfn,nsft,cfn,3,xxs,ffs,1,ier)
    dens=ffs(1)
    nu=1
    ier=0
    call splder(tfn,nsft,cfn,3,nu,xxs,ffs,1,wrkfd,ier)
    ddens=ffs(1)
    if(ier>0) print*,ier
  end if
  if(dens<zero) then
    print*,'DENSITY NEGATIVE',dens
    dens=zero
  end if
end subroutine density

function temp(arg)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TEMP(arg)                                                            C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : npmx
  implicit none
  real(kind=r8) temp
  real(kind=r8), intent(in) :: arg
  integer k,nr,nz,nppn,nppt
  real(kind=r8) psradt(npmx),terad(npmx),ct(npmx,4)
  real(kind=r8) temp0,dps
  real(kind=r8) psntpp

  common/crad/psradt,terad
  common/temppp/psntpp
  common/coffte/ct
  common/eqnn/nr,nz,nppn,nppt

  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    subroutine intlin(x1,y1,x2,y2,x,y)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: x1, y1, x2, y2, x
      real(kind=r8), intent(out) :: y
    end subroutine intlin
    function spli(cspli,k,dx)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: spli
      integer      , intent(in) :: k
      real(kind=r8), intent(in) :: cspli(:,:), dx
    end function spli
  end interface

!     spline interpolation for n=0,N-2. Linear for n>=N-1. Te(N)=0
  if (arg>=psntpp) then
    temp0=0.0d0
  else if (arg<=0.0d0) then
    temp0=terad(1)
  else
    k=locate(psradt(:nppt),arg)
    if(k==0) k=1
    if (k<nppt-1) then
      dps=arg-psradt(k)
      temp0=spli(ct,k,dps)
    else
      k=nppt-1
      call intlin(psradt(k),terad(k),psradt(nppt),terad(nppt),arg,temp0)
    endif
    if (temp0<0.0d0) temp0=0.0d0
  endif
  temp=temp0  
end function temp

subroutine icrtrn(j,k,st0,ierr)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ICRTRN(j,k) computes initial conditions for multiple ray tracing     C
!             gaussian far-field approximation                         C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi
  use spece_data, only : nmx
  implicit none
  integer, intent(in) :: j, k
  real(kind=r8), intent(out) :: st0
  integer, intent(out) :: ierr
  integer, parameter :: ndim = 6
  integer nray,ktx
  real(kind=r8) ywrk0(ndim),ypwrk0(ndim)
  real(kind=r8) rmx,ak0,w0xt,pw0xt,w0yt,pw0yt,awr
  real(kind=r8) x00,y00,z00,anx0c,any0c,anz0c,anx0,any0,anz0
  real(kind=r8) csth,snth,csps,snps,csawr,snawr,dd,an2s0,an20
  real(kind=r8) alx,aly,dzx,dzy,wxm,wym,dr,da,z0t,rj,alfak
  real(kind=r8) x0tr,y0tr,x0t,y0t,dx0,dy0,dz0,x0,y0,z0
  real(kind=r8) ppxr,ppyr,anxtr,anytr,anztr,anxt,anyt,anzt
  real(kind=r8) conex,coney
  real(kind=r8) :: psinv,yg,anpl,anpl2,anpr,xg
  real(kind=r8), dimension(nmx) :: ygv, psjki, rrr, anplv, stv, zzv, &
                                   phiv, xgv, anprv
  integer iproj

  common/nrayktx/nray,ktx
  common/rhomx/rmx
  common/ak0/ak0
  common/parbeam/w0xt,pw0xt,w0yt,pw0yt,awr
  common/anic/anx0c,any0c,anz0c
  common/wrk/ywrk0,ypwrk0
  common/mirr/x00,y00,z00
  common/ddd/dd,an2s0,an20
  common/ygv/ygv
  common/psjki/psjki
  common/psival/psinv
  common/rrr/rrr
  common/ygyg/yg
  common/nplr/anpl,anpl2,anpr
  common/nplv/anplv
  common/xgxg/xg
  common/traj/stv,zzv,phiv,xgv,anprv
  common/iproj/iproj

  interface
    subroutine goto_bound(x0, y0, z0, nx, ny, nz, st, ierr)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(inout) :: x0, y0, z0
      real(kind=r8), intent(in)    :: nx, ny, nz
      real(kind=r8), intent(out)   :: st
      integer,       intent(out)   :: ierr
    end subroutine goto_bound
    subroutine fwork(y,dery)
      use kinds, only : r8
      implicit none
      integer, parameter :: ndim = 6
      real(kind=r8), intent(in)  :: y(ndim)
      real(kind=r8), intent(out) :: dery(ndim)
    end subroutine fwork
    integer function islhfs(r,z)
      real(kind=8), intent(in) :: r, z
    end function islhfs
  end interface
! vector refraction index N directed towards the antenna
  csth=-anz0c
  snth=sqrt(1.0d0-csth**2)
  csps=1.0d0
  snps=0.0d0
  if(snth>0.0d0) then
    csps=-any0c/snth
    snps=-anx0c/snth
  end if
  csawr=cos(awr*pi/180.0d0)
  snawr=sin(awr*pi/180.0d0)
  z0t=0.0d0
  alx=ak0*w0xt**2/2.0d0
  aly=ak0*w0yt**2/2.0d0
  conex=w0xt/alx
  coney=w0yt/aly
  dzx=-pw0xt+z0t
  dzy=-pw0yt+z0t
  wxm=conex*dzx
  wym=coney*dzy
  dr=1.0d0
  if(nray>1) dr=rmx/dble(nray-1)
  da=2.0d0*pi/dble(ktx)
  z0t=0.0d0
  rj=(j-1)*dr
  alfak=(k-1)*da
  x0tr=rj*cos(alfak)*wxm
  y0tr=rj*sin(alfak)*wym
  x0t=x0tr*csawr-y0tr*snawr
  y0t=x0tr*snawr+y0tr*csawr
  dx0= x0t*csps+snps*(y0t*csth+z0t*snth)
  dy0=-x0t*snps+csps*(y0t*csth+z0t*snth)
  dz0= z0t*csth-y0t*snth
  x0=x00+dx0
  y0=y00+dy0
  z0=z00+dz0
  ppxr=rj*cos(alfak)*conex
  ppyr=rj*sin(alfak)*coney
  anztr=1.0d0/sqrt(1.0d0+ppxr**2+ppyr**2)
  anxtr=ppxr*anztr
  anytr=ppyr*anztr
  anxt=anxtr*csawr-anytr*snawr
  anyt=anxtr*snawr+anytr*csawr
  anzt=anztr
  anx0=-( anxt*csps+snps*(anyt*csth+anzt*snth))
  any0=-(-anxt*snps+csps*(anyt*csth+anzt*snth))
  anz0=-( anzt*csth-anyt*snth)
  ywrk0(1)=x0
  ywrk0(2)=y0
  ywrk0(3)=z0
  ywrk0(4)=anx0
  ywrk0(5)=any0
  ywrk0(6)=anz0
  st0=0.0d0
  ierr=0
!  call goto_bound(x0,y0,z0,anx0,any0,anz0,st0,ierr)
!  ywrk0(1)=x0
!  ywrk0(2)=y0
!  ywrk0(3)=z0
  an20=1.0d0
  an2s0=anx0**2+any0**2+anz0**2
!  dd=an2s0-an20
  call fwork(ywrk0,ypwrk0)
  rrr(1)=sqrt(x0**2+y0**2)
  zzv(1)=z0*1.0d-2
  stv(1)=st0*1.0d-2
  phiv(1)=acos(x0/rrr(1))
  if(y0<0.0d0) phiv(1)=-phiv(1)
  phiv(1)=phiv(1)*180.0d0/pi
  if(iproj.eq.1) then
    if(islhfs(rrr(1)*1.d-2,zzv(1))==1) then
      psjki(1)=psinv
    else
      psjki(1)=-psinv
    end if
  else
    psjki(1)=psinv
  end if
  xgv(1)=xg
  ygv(1)=yg
  anplv(1)=anpl
  anprv(1)=anpr
end subroutine icrtrn

subroutine pweight
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PWEIGHT                                                              C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : jmx
  implicit none
  integer j,nray,ktx
  real(kind=r8) q(jmx)
  real(kind=r8) rmx,r1,r2,dr,summ

  common/qw/q
  common/nrayktx/nray,ktx
  common/rhomx/rmx

  dr=1.0d0
  if(nray>1) dr=rmx/dble(nray-1)
  r1=0.0d0
  summ=0.0d0
  q(1)=1.0d0
  if(nray>1) then
    do j=1,nray
        r2=(dble(j)-0.5d0)*dr
      q(j)=(exp(-2.0d0*r1**2)-exp(-2.0d0*r2**2))
      r1=r2
      summ=summ+q(j)
    end do
    q(1)=q(1)/summ
    do j=2,nray
      q(j)=q(j)/dble(ktx)/summ
    end do
  end if
end subroutine pweight

subroutine emission(i)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! EMISSION(i)                                                          C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx, ncmx
  implicit none
  integer, intent(in) :: i
  integer :: nharm,iovlp,iharm
  real(kind=r8) tauv(nmx),alphav(nmx)
  real(kind=r8) tev(nmx),tta(nmx),rrr(nmx)
  real(kind=r8) tekev(ncmx),terad,ttar1,ttar2,dst,alpha
  real(kind=r8) tau0,dtau,gg0,gg1,ggv(nmx)
  real(kind=r8) vgm,derdnm,ds
  integer nharmv(nmx)

  common/atjki/tauv,alphav,tev
  common/tta/tta
  common/nharm/nharm
  common/harm/iovlp,iharm
  common/harmv/nharmv
  common/tete/tekev,terad
  common/rrr/rrr
  common/ttar/ttar1,ttar2
  common/dsds/dst
  common/absa/alpha,tau0
  common/atemt/ggv
  common/vgv/vgm,derdnm

  interface
    subroutine ecrh(i)
      integer, intent(in) :: i
    end subroutine ecrh
  end interface
  
  ds=dst*vgm
  tau0=tauv(i-1)
  call ecrh(i)
  alphav(i)=alpha
  tev(i)=terad
  dtau=(alphav(i)+alphav(i-1))*ds/2.0d0
  tauv(i)=tauv(i-1)+dtau
  gg1=alphav(i)*tev(i)*exp(-tauv(i))
  gg0=alphav(i-1)*tev(i-1)*exp(-tauv(i-1))
  ggv(i)=gg1
  tta(i)=tta(i-1)+(gg0+gg1)*ds/2.0d0
  ttar1=ttar1+(gg0*rrr(i-1)+gg1*rrr(i))*ds/2.0d0
  ttar2=ttar2+(gg0*rrr(i-1)**2+gg1*rrr(i)**2)*ds/2.0d0
  nharmv(i)=nharm
  if(alpha>1d-8) then
    if(iharm==99) iharm=nharm
    if(iharm>nharm) then
      iharm=nharm
      iovlp=iovlp+1
    end if
  end if
end subroutine emission

subroutine ecrh(i)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ECRH(i)                                                              C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : ncmx
  implicit none
  integer, intent(in) :: i
  logical first,okcall
  integer iwarm,ilarm,nharm,nh,nnh,lrm,ncmp,j,icall,nhblk
  real(kind=r8), parameter :: xxcr = 100.0d0, eps = 1.d-8
  integer, parameter :: lmx = 10
  real(kind=r8) yg,anpl,anpl2,anpr,vgm,derdnm,ak0,anprre,anprim,alpha
  real(kind=r8) alphai,tau,psinv
  real(kind=r8) tekev(ncmx),terad,eta(ncmx),amu(ncmx),ygn,ygc
  real(kind=r8) rdu,u1,u2,uu2,argex,g1,g2,gg,argexp
  real(kind=r8) ygl,anpll,dnll,udr(ncmx),gam(ncmx),gdoppl(ncmx)
  complex(kind=r8) exex,eyey,ezez,exey,exez,eyez
  complex(kind=r8) ea11(0:ncmx),ea22(0:ncmx),ea33(0:ncmx)
  complex(kind=r8) ea12(0:ncmx),ea13(0:ncmx),ea23(0:ncmx)

  common/ndrf/ncmp
  common/nharm/nhblk
  common/warm/iwarm,ilarm
  common/ygyg/yg
  common/nplr/anpl,anpl2,anpr
  common/vgv/vgm,derdnm
  common/ak0/ak0
  common/nprw/anprre,anprim
  common/absa/alpha,tau
  common/epolar/exex,eyey,ezez,exey,exez,eyez
  common/eanti/ea11,ea22,ea33,ea12,ea13,ea23
  common/psival/psinv
  common/tete/tekev,terad
  common/multc/eta,amu
  common/udrf/udr,gam
  common/doppl/gdoppl
  common/icall/icall

  interface
    subroutine fdd_composition(psinv)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: psinv
    end subroutine fdd_composition
    subroutine dispersion(lrm,first)
      implicit none
      integer, intent(in) :: lrm
      logical, intent(in) :: first
    end subroutine dispersion
  end interface
!     absorption computation
  alpha=0.0d0
  terad=0.0d0
  lrm=ilarm
  okcall=.false.
  call fdd_composition(psinv)
  do j=1,ncmp
    gdoppl(j)=gam(j)-anpl*udr(j)
    if (tekev(j)>0.0d0) then
      ygl=yg/gdoppl(j)
      anpll=(gam(j)*anpl-udr(j))/gdoppl(j)
      dnll=1.0d0-anpll*anpll
      if(iwarm==1) then
        ygc=1.0d0-anpll*anpll/2.0d0
      else
        ygc=sqrt(dnll)
      end if
      nharm=int(ygc/ygl-eps)+1
      ygn=nharm*ygl
      if(iwarm==1) then
        rdu=sqrt(2.0d0*(ygn-ygc))
        u1=anpll+rdu
        u2=anpll-rdu
        uu2=min(u1*u1,u2*u2)
        argex=amu(j)*uu2/2.0d0
      else
        rdu=sqrt(ygn**2-ygc**2)
        g1=(ygn+anpll*rdu)/dnll
        g2=(ygn-anpll*rdu)/dnll
        gg=min(g1,g2)
        argex=amu(j)*(gg-1.0d0)
      end if
!      if (nharm>0.and.argex<=xxcr) then
      nh=nharm+1
      ygn=nh*ygl
      if(iwarm==1) then
        rdu=sqrt(2.0d0*(ygn-ygc))
        u1=anpll+rdu
        u2=anpll-rdu
        uu2=min(u1*u1,u2*u2)
        argexp=amu(j)*uu2/2.0d0
      else
        rdu=sqrt(ygn**2-ygc**2)
        g1=(ygn+anpll*rdu)/dnll
        g2=(ygn-anpll*rdu)/dnll
        gg=min(g1,g2)
        argexp=amu(j)*(gg-1.0d0)
      end if
      if(argexp<argex) then
        nnh=nh
        argex=argexp
      else
        nnh=nharm
      end if
      if (j==1) nhblk=nnh
      if (argex<=xxcr) then
        lrm = max(lrm,nnh+ilarm)
        okcall=.true.
      end if
    end if
  end do
  if (lrm.gt.lmx) lrm=lmx
  if (okcall) then
    if (i-icall==1) then
      first=.false.
    else
      first=.true.
    end if
    call dispersion(lrm,first)
    icall=i
    alpha=4.0d0*ak0*anprim*anpr/derdnm
    if (alpha<0.0d0) alpha=0.0d0
    if (ncmp>1) then
      do j=1,ncmp
        alphai=dble(exex*ea11(j)+eyey*ea22(j)+ezez*ea33(j)+         &
                2.0d0*(exey*ea12(j)+exez*ea13(j)+eyez*ea23(j)))
        terad=terad+alphai*tekev(j)/gdoppl(j)
      end do
      alphai=dble(exex*ea11(0)+eyey*ea22(0)+ezez*ea33(0)+           &
              2.0d0*(exey*ea12(0)+exez*ea13(0)+eyez*ea23(0)))
      if (alphai>0.0d0) then
        terad=terad/alphai
      else
        terad=tekev(1)
      end if
    else
      terad=tekev(1)
    end if
  else
!    write(33,*) ' '
    terad=tekev(1)
  end if
end subroutine ecrh

subroutine dispersion(lrm,first)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DISPERSION(lrm,first)                                                C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi, huge_r8, tiny_r8
  use spece_data, only : ncmx
  implicit none
  integer, intent(in) :: lrm
  logical, intent(in) :: first
  integer versign,versign0
  integer imx,imxx,iwarm,ilarm,i,j,k,ilrm,ncmp,nc,ierr
  complex(kind=r8) cc0,cc2,cc4,rr,dc0,dc2,dc4,derdd
  complex(kind=r8) anpr2a,anpr20
  complex(kind=r8) anpr,anpr2,ex,ey,ez,den
  complex(kind=r8) epsl0(3,3,lrm,ncmx),epslt(3,3,lrm)
  complex(kind=r8) ssepsl(3,3),sepsl(3,3,ncmx),sepsl0(3,3,ncmx)
  complex(kind=r8) sdepsl(3,3),depsl(3,3,ncmx),depsl0(3,3,ncmx)
  complex(kind=r8) e330(ncmx),se330,e11,e22,e12,e13,e23,a13,a23,a33
  complex(kind=r8) exex,eyey,ezez,exey,exez,eyez
  complex(kind=r8) ea11(0:ncmx),ea22(0:ncmx),ea33(0:ncmx)
  complex(kind=r8) ea12(0:ncmx),ea13(0:ncmx),ea23(0:ncmx)
  real(kind=r8) xg,yg,anpl,anprf,sox,anprr,anpri,errnpr,anpl2,sgn,anpr2f
  real(kind=r8) eta(ncmx),amu(ncmx),udr(ncmx),gam(ncmx),gdoppl(ncmx)
  real(kind=r8) amul,xgl,ygl,anpll,anpl2l
  real(kind=r8) ey2,ez2,ex2
  real(kind=r8) fghz
  real(kind=r8) angdrr,sgn0
  complex(kind=r8) rr0

  common/ndrf/ncmp
  common/ierr/ierr
  common/xgxg/xg
  common/ygyg/yg
  common/nplr/anpl,anpl2,anprf
  common/mode/sox
  common/fghz/fghz
  common/warm/iwarm,ilarm
  common/idisp/imx
  common/nprw/anprr,anpri
  common/epolar/exex,eyey,ezez,exey,exez,eyez
  common/eanti/ea11,ea22,ea33,ea12,ea13,ea23
  common/multc/eta,amu
  common/udrf/udr,gam
  common/doppl/gdoppl
  common/muloc/amul
  common/xloc/xgl
  common/yloc/ygl
  common/npllc/anpll,anpl2l
  common/vs/versign0
  save sgn0,rr0

  interface
    subroutine diel_tens_fr(e330,epsl,lrm)
      use kinds, only : r8
      implicit none
      integer         , intent(in)  :: lrm
      complex(kind=r8), intent(out) :: e330, epsl(3,3,lrm)
    end subroutine diel_tens_fr
    subroutine diel_tens_wr(ce330,cepsl,lrm)
      use kinds, only : r8
      implicit none
      integer         , intent(in)  :: lrm
      complex(kind=r8), intent(out) :: ce330, cepsl(3,3,lrm)
    end subroutine diel_tens_wr
  end interface

  se330=cmplx(1.0d0-ncmp,0.0d0)
  do nc=1,ncmp
    amul=amu(nc)
    xgl=xg*eta(nc)/(gam(nc)*gdoppl(nc)*gdoppl(nc))
    ygl=yg/gdoppl(nc)
    anpll=(gam(nc)*anpl-udr(nc))/gdoppl(nc)
    anpl2l=anpll*anpll
    if (iwarm==1) then
      call diel_tens_wr(e330(nc),epslt,lrm)
    else
      call diel_tens_fr(e330(nc),epslt,lrm)
    end if
    do j=1,3
      do k=1,j
        do ilrm=1,lrm
          epsl0(k,j,ilrm,nc)=epslt(k,j,ilrm)
        end do
      end do
    end do
    se330=se330+e330(nc)
  end do
!  errnpr=1.0d0
  anpr2f=anprf*anprf
  anpr2a=anpr2f
  if (versign0==2 .and. first) then
    versign=1
  else
    versign=versign0
  end if
  imxx=abs(imx)
  do ! loop to disable convergence if failure detected
  do i=1,imxx
    do nc=1,ncmp
      anpr20=anpr2a/(gdoppl(nc)*gdoppl(nc))
      do j=1,3
        do k=1,j
          sepsl0(k,j,nc)=cmplx(0.0d0,0.0d0)
          depsl0(k,j,nc)=cmplx(0.0d0,0.0d0)
          do ilrm=1,lrm
            sepsl0(k,j,nc)=sepsl0(k,j,nc)+epsl0(k,j,ilrm,nc)*anpr20**(ilrm-1)
            if (ilrm.gt.1) depsl0(k,j,nc)=depsl0(k,j,nc) &
                    +(ilrm-1)*epsl0(k,j,ilrm,nc)*anpr20**(ilrm-2)
          end do
        end do
      end do
      sepsl(1,1,nc)=1.0d0+gdoppl(nc)*gdoppl(nc)*(sepsl0(1,1,nc)-1.0d0)
      sepsl(1,2,nc)=gdoppl(nc)*gdoppl(nc)*sepsl0(1,2,nc)
      sepsl(2,2,nc)=1.0d0+gdoppl(nc)*gdoppl(nc)*(sepsl0(2,2,nc)-1.0d0)
      sepsl(1,3,nc)=sepsl0(1,3,nc)+udr(nc)*gdoppl(nc)*(sepsl0(1,1,nc)-1.0d0)
      sepsl(2,3,nc)=sepsl0(2,3,nc)-udr(nc)*gdoppl(nc)*sepsl0(1,2,nc) !eps12=-eps21
      sepsl(3,3,nc)=sepsl0(3,3,nc)/(gdoppl(nc)*gdoppl(nc))          &
                    +sepsl0(1,3,nc)*2.0d0*udr(nc)/gdoppl(nc)        &
                    +udr(nc)*udr(nc)*(sepsl0(1,1,nc)-1.0d0)
      depsl(1,1,nc) = depsl0(1,1,nc)
      depsl(1,2,nc) = depsl0(1,2,nc)
      depsl(2,2,nc) = depsl0(2,2,nc)
      depsl(1,3,nc) = (depsl0(1,3,nc)/gdoppl(nc) + udr(nc)*depsl0(1,1,nc))/gdoppl(nc)
      depsl(2,3,nc) = (depsl0(2,3,nc)/gdoppl(nc) - udr(nc)*depsl0(1,2,nc))/gdoppl(nc)
      depsl(3,3,nc) = (depsl0(3,3,nc)/(gdoppl(nc)*gdoppl(nc))        &
                      +depsl0(1,3,nc)*2.0d0*udr(nc)/gdoppl(nc)       &
                      +udr(nc)*udr(nc)*depsl0(1,1,nc))/(gdoppl(nc)*gdoppl(nc))
    end do
    do j=1,3
      do k=1,j
        if((k==j).and.(k/=3)) then
          ssepsl(k,j)=cmplx(1.0d0-ncmp,0.0d0)
        else
          ssepsl(k,j)=cmplx(0.0d0,0.0d0)
        end if
        sdepsl(k,j)=cmplx(0.0d0,0.0d0)
        do nc=1,ncmp
          ssepsl(k,j)=ssepsl(k,j)+sepsl(k,j,nc)
          sdepsl(k,j)=sdepsl(k,j)+depsl(k,j,nc)
        end do
      end do
    end do
    e11=ssepsl(1,1)
    e22=ssepsl(2,2)
    e12=ssepsl(1,2)
    a33=ssepsl(3,3)
    a13=ssepsl(1,3)
    a23=ssepsl(2,3)
    cc4 = (e11-anpl2)*(1.0d0-a33) + (a13+anpl)*(a13+anpl)
    cc2 = e12*e12*(a33-1.0d0) + 2.0d0*a23*e12*(a13+anpl)            &
        + (a23*a23 - se330 + (e22-anpl2)*(a33-1.0d0))*(e11-anpl2)   &
        + (a13+anpl)*(a13+anpl)*(anpl2-e22)
    cc0 = se330*((e11-anpl2)*(e22-anpl2) + e12*e12)
    rr = cc2*cc2-4.0d0*cc0*cc4
    select case (versign)
      case (1) ! ---------- original sign choice --------------
        if(yg>1.0d0) then
          sgn=sox
          if(aimag(rr)<=0.0d0) sgn=-sgn
        else
          sgn=-sox
          if(aimag(rr)>=0.d0) sgn=-sgn
        end if
        anpr2=(-cc2 + sgn*sqrt(rr))/(2.0d0*cc4)
      case (2) ! ***** follows previous determination *****
        angdrr=abs(atan2(aimag(rr),dble(rr))-atan2(aimag(rr0),dble(rr0)))
        if (angdrr<pi) then
          sgn=sgn0
        else
          angdrr=2.0d0*pi-angdrr
          sgn=-sgn0
        end if
        anpr2=(-cc2 + sgn*sqrt(rr))/(2.0d0*cc4)
      case default !---- solution found with Newton's method ---
        sgn=sox
        dc0 = se330*( (e22-anpl2)*sdepsl(1,1) + 2.0d0*e12*sdepsl(1,2)    &
                     +(e11-anpl2)*sdepsl(2,2) )
        dc2 = (-se330 + a23**2 + (a33-1.0d0)*(e22-anpl2) )*sdepsl(1,1)   &
              + 2.0d0*( a23*(a13+anpl) + e12*(a33-1.0d0) )*sdepsl(1,2)   &
              + 2.0d0*( e12*a23 - (a13+anpl)*(e22-anpl2) )*sdepsl(1,3)   &
            + ( -(a13+anpl)**2 + (a33-1.0d0)*(e11-anpl2) )*sdepsl(2,2)   &
              + 2.0d0*( e12*(a13+anpl) + a23*(e11-anpl2) )*sdepsl(2,3)   &
                    + ( e12**2 + (e11-anpl2)*(e22-anpl2) )*sdepsl(3,3)
        dc4 = (1.0d0-a33)*sdepsl(1,1) + 2.0d0*(a13+anpl)*sdepsl(1,3)     &
             -(e11-anpl2)*sdepsl(3,3)
        derdd = (dc4*anpr2a + (2.0d0*cc4 + dc2))*anpr2a + cc2 + dc0
        anpr2 = anpr2a - ((cc4*anpr2a + cc2)*anpr2a + cc0)/derdd
    end select
    errnpr=abs(1.0d0-abs(anpr2)/abs(anpr2a))
    if(i>1.and.errnpr<1.0d-5) exit
    anpr2a=anpr2
  end do
  if(i>imxx .and. imxx>1) then
    if (imx<0) then
      write(*,"('f =',f7.2,' X =',f7.4,' Y =',f7.4,' Nperp =',f7.4,"// &
        "': convergence disabled.')") fghz,xg,yg,sqrt(abs(anpr2))
      imxx=1
      anpr2a=anpr2f
    else
      write(*,"('f =',f7.2,' X =',f7.4,' Y =',f7.4,' Nperp =',f7.4,"// &
        "': convergence failed.')") fghz,xg,yg,sqrt(abs(anpr2))
      exit
    end if
  else
    exit
  end if
  end do ! end loop to disable convergence
  if(sqrt(dble(anpr2))<0.0d0 &
    .or. anpr2/=anpr2 .or. abs(anpr2)>=huge_r8 .or. abs(anpr2)<=tiny_r8) then
    write(*,"('f =',f7.2,' X =',f7.4,' Y =',f7.4,' Nperp =',f7.4,'!')") &
        fghz,xg,yg,sqrt(abs(anpr2))
    ierr=99
    anpr2=(0.0d0,0.0d0)
  end if
  rr0=rr
  sgn0=sgn
  anpr=sqrt(anpr2)
  anprr=dble(anpr)
  anpri=aimag(anpr)
!  write(33,'(8(1x,e16.9),i4)') fghz,xg,yg,amu(1),anpl,anprf,anprr,anpri,i
  e13=anpr*a13
  e23=anpr*a23
  ex=cmplx(0.0d0,0.0d0)
  ey=cmplx(0.0d0,0.0d0)
  ez=cmplx(0.0d0,0.0d0)
  ex2=0.0d0
  ey2=0.0d0
  ez2=0.0d0
  if (abs(anpl)>1.0d-6) then
    den=e12*e23-(e13+anpr*anpl)*(e22-anpr2-anpl2)
    ey=-(e12*(e13+anpr*anpl)+(e11-anpl2)*e23)/den
    ez=(e12*e12+(e22-anpr2-anpl2)*(e11-anpl2))/den
    ez2=abs(ez)**2
    ey2=abs(ey)**2
    ex2=1.0d0/(1.0d0+ey2+ez2)
    ex=cmplx(sqrt(ex2),0.0d0)
    ez2=ez2*ex2
    ey2=ey2*ex2
    ez=ez*ex
    ey=ey*ex
  else if (sox<0.0d0) then
    ez=cmplx(1.0d0,0.0d0)
    ez2=abs(ez)**2
  else
    ex2=1.0d0/(1.0d0+abs(-e11/e12)**2)
    ex=sqrt(ex2)
    ey=-ex*e11/e12
    ey2=abs(ey)**2
  end if
  exex=ex2
  eyey=ey2
  ezez=ez2
  exey=conjg(ex)*ey
  exez=conjg(ex)*ez
  eyez=conjg(ey)*ez
  do nc=1,ncmp
    ea11(nc)=aimag(sepsl(1,1,nc))
    ea22(nc)=aimag(sepsl(2,2,nc))
    ea33(nc)=aimag(anpr2*sepsl(3,3,nc)+e330(nc))
    ea12(nc)=cmplx(0.0d0,-1.0d0)*dble(sepsl(1,2,nc))
    ea13(nc)=aimag(anpr*sepsl(1,3,nc))
    ea23(nc)=cmplx(0.0d0,-1.0d0)*dble(anpr*sepsl(2,3,nc))
  end do
  ea11(0)=aimag(e11)
  ea22(0)=aimag(e22)
  ea33(0)=aimag(anpr2*a33+se330)
  ea12(0)=cmplx(0.0d0,-1.0d0)*dble(e12)
  ea13(0)=aimag(e13)
  ea23(0)=cmplx(0.0d0,-1.0d0)*dble(e23)
end subroutine dispersion

subroutine diel_tens_fr(e330,epsl,lrm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DIEL_TENS_FR(e330,epsl,lrm) computes fully relativistic dielectric   C
!        tensor elements. Expansion to lrm-th order in Larmor radius   C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi, pi_sqrt
  implicit none
  integer         , intent(in)  :: lrm
  complex(kind=r8), intent(out) :: e330, epsl(3,3,lrm)
  integer i,j,l,lm,is,k
  complex(kind=r8) ui
  complex(kind=r8) ca11,ca12,ca22,ca13,ca23,ca33
  complex(kind=r8) cq0p,cq0m,cq1p,cq1m,cq2p
  real(kind=r8) rr(-lrm:lrm,0:2,0:lrm),ri(lrm,0:2,lrm)
  real(kind=r8) xg,yg,amu,cr,ci,cmxw
  real(kind=r8) fal,w,asl,bsl

  common/muloc/amu
  common/xloc/xg
  common/yloc/yg
  common/cri/cr,ci

  interface
    function fact(k)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fact
      integer, intent(in) :: k
    end function fact
    subroutine hermitian(rr,lrm)
      use kinds, only : r8
      implicit none
      integer      , intent(in)  :: lrm
      real(kind=r8), intent(out) :: rr(-lrm:lrm,0:2,0:lrm)
    end subroutine hermitian
    subroutine antihermitian(ri,lrm)
      use kinds, only : r8
      implicit none
      integer      , intent(in)  :: lrm
      real(kind=r8), intent(out) :: ri(lrm,0:2,lrm)
    end subroutine antihermitian
  end interface

  ui=cmplx(0.0d0,1.0d0)
  cmxw=1.0d0+15.0d0/(8.0d0*amu)+105.0d0/(128.0d0*amu**2)
  cr=-amu*amu/(pi_sqrt*cmxw)
  ci=sqrt(2.0d0*pi*amu)*amu**2/cmxw
  do i=1,3
    do j=1,3
      do l=1,lrm
        epsl(i,j,l)=cmplx(0.0d0,0.0d0)
      end do
    end do
  end do
  call hermitian(rr,lrm)
  call antihermitian(ri,lrm)
  do l=1,lrm
    lm=l-1
    fal=-0.25d0**l*fact(2*l)/(fact(l)**2*yg**(2*lm))
    ca11=(0.0d0,0.0d0)
    ca12=(0.0d0,0.0d0)
    ca13=(0.0d0,0.0d0)
    ca22=(0.0d0,0.0d0)
    ca23=(0.0d0,0.0d0)
    ca33=(0.0d0,0.0d0)
    do is=0,l
      k=l-is
      w=(-1.0d0)**k
      asl=w/(fact(is+l)*fact(l-is))
      bsl=asl*(is*is+dble(2*k*lm*(l+is))/(2.0d0*l-1.0d0))
      if(is>0) then
        cq0p=rr(is,0,l)+rr(-is,0,l)+ui*ri(is,0,l)
        cq0m=rr(is,0,l)-rr(-is,0,l)+ui*ri(is,0,l)
        cq1p=rr(is,1,l)+rr(-is,1,l)+ui*ri(is,1,l)
        cq1m=rr(is,1,l)-rr(-is,1,l)+ui*ri(is,1,l)
        cq2p=rr(is,2,l)+rr(-is,2,l)+ui*ri(is,2,l)
      else
        cq0p=rr(is,0,l)
        cq0m=rr(is,0,l)
        cq1p=rr(is,1,l)
        cq1m=rr(is,1,l)
        cq2p=rr(is,2,l)
      end if
      ca11=ca11+is**2*asl*cq0p
      ca12=ca12+is*l*asl*cq0m
      ca22=ca22+bsl*cq0p
      ca13=ca13+is*asl*cq1m/yg
      ca23=ca23+l*asl*cq1p/yg
      ca33=ca33+asl*cq2p/yg**2
    end do
    epsl(1,1,l) =  - xg*ca11*fal
    epsl(1,2,l) =  + ui*xg*ca12*fal
    epsl(2,2,l) =  - xg*ca22*fal
    epsl(1,3,l) =  - xg*ca13*fal
    epsl(2,3,l) =  - ui*xg*ca23*fal
    epsl(3,3,l) =  - xg*ca33*fal
  end do
  cq2p=rr(0,2,0)
  e330=1.0d0+xg*cq2p
  epsl(1,1,1) = 1.d0 + epsl(1,1,1)
  epsl(2,2,1) = 1.d0 + epsl(2,2,1)
  do l=1,lrm
    epsl(2,1,l) = - epsl(1,2,l)
    epsl(3,1,l) =   epsl(1,3,l)
    epsl(3,2,l) = - epsl(2,3,l)
  end do
end subroutine diel_tens_fr

subroutine hermitian(rr,lrm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! HERMITIAN(r,lrm)                                                     C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : npts, tmax
  implicit none
  integer      , intent(in)  :: lrm
  real(kind=r8), intent(out) :: rr(-lrm:lrm,0:2,0:lrm)
  integer      , parameter :: nptp1 = npts+1
  integer n,k,m,llm,i
  real(kind=r8) t,x,dt
  real(kind=r8) yg,amu,anpl,cr,ci
  real(kind=r8) upl,upl2,gx,exdx,gr,s,ffe,zm,fe0m
  real(kind=r8) sy1,sy2,sy3,bth2,bth4,bth6,anpl2
  real(kind=r8) rxt,bth,amu2,amu4,amu6,zm2,zm3
  real(kind=r8) ttv(nptp1),extv(nptp1)

  common/ttex/ttv,extv
  common/muloc/amu
  common/yloc/yg
  common/npllc/anpl,anpl2
  common/cri/cr,ci
  
  interface
    function expei(x)
      use kinds, only : r8
      implicit none
      real(kind=r8)              :: expei
      real(kind=r8), intent(in)  :: x
    end function expei
  end interface

  do n=-lrm,lrm
    do k=0,2
      do m=0,lrm
        rr(n,k,m)=0.0d0
      end do
    end do
  end do
  llm=min(3,lrm)
  dt=2.0d0*tmax/dble(npts)
  bth2=2.0d0/amu
  bth=sqrt(bth2)
  amu2=amu*amu
  amu4=amu2*amu2
  amu6=amu4*amu2
  do i = 1, nptp1
    t = ttv(i)
    rxt=sqrt(1.0d0+t*t/(2.0d0*amu))
    x = t*rxt
    upl2=bth2*x**2
    upl=bth*x
    gx=1.0d0+t*t/amu
    exdx=cr*extv(i)*gx/rxt*dt
    do n=1,llm
      gr=anpl*upl+n*yg
      zm=-amu*(gx-gr)
      s=amu*(gx+gr)
      zm2=zm*zm
      zm3=zm2*zm
      fe0m=expei(zm)
      do m=n,llm
        if (m==1) ffe=(1.0d0+s*(1.0d0-zm*fe0m))/amu2
        if (m==2) ffe=(6.0d0-2.0d0*zm+4.0d0*s+s*s*(1.0d0+zm-zm2*fe0m))/amu4
        if (m==3) ffe=(18.0d0*s*(s+4.0d0-zm)+6.0d0*                &
              (20.0d0-8.0d0*zm+zm2)+s**3*(2.0d0+zm+zm2-zm3*fe0m))/amu6
        rr(n,0,m) = rr(n,0,m) + exdx*ffe
        rr(n,1,m) = rr(n,1,m) + exdx*ffe*upl
        rr(n,2,m) = rr(n,2,m) + exdx*ffe*upl2
      end do
    end do
  end do
  sy1=1.0d0+yg
  sy2=1.0d0+yg*2.0d0
  sy3=1.0d0+yg*3.0d0
  bth4=bth2*bth2
  bth6=bth4*bth2
  rr(0,2,0) = -(1.0d0+bth2*(6.0d0*anpl2-5.0d0)/4.0d0                 &
            +(120.d0*anpl2*anpl2-192.0d0*anpl2+55.0d0)*bth4/32.0d0)
  rr(0,1,1) = -anpl*bth2*(1.0d0+0.75d0*(2.0d0*anpl2-3.0d0)*bth2)
  rr(0,2,1) = -bth2-0.5d0*bth4*(3.0d0*anpl2-1.0d0)
  rr(-1,0,1) = -2.0d0/sy1*(1.0d0+(-5.0d0*sy1+2.0d0*anpl2)            &
                *bth2/sy1**2/4.0d0+(-15.0d0/sy1-10.0d0*              &
                (7.0d0+2.0d0*anpl2)/sy1**2+24.0d0*anpl2*anpl2/sy1**4 &
                 -84.0d0*anpl2/sy1**3)*bth4/32.0d0)
  rr(-1,1,1) = -anpl*bth2/sy1**2*(1.0d0+bth2*(6.0d0*anpl2            &
                +5.0d0*sy1**2-14.0d0*sy1)/sy1**2/4.0d0)
  rr(-1,2,1) = -bth2/sy1*(1.0d0+bth2*(6.0d0*anpl2+5.0d0*sy1**2       &
               -7.0d0*sy1)/sy1**2/4.0d0)
  if (llm.eq.1) return
  rr(0,0,2) = -(4.0d0*bth2-2.0d0*bth4*(1.0d0-anpl2)                  &
                +1.5d0*bth6*(3.0d0-5.0d0*anpl2+2.0d0*anpl2*anpl2))
  rr(0,1,2) = -2.0d0*anpl*bth4-3.0d0*anpl*(1.0d0-anpl2)*bth6
  rr(0,2,2) = -2.0d0*bth4-1.5d0*bth6*(1.0d0+2.0d0*anpl2)
  rr(-1,0,2) = -4.0d0*bth2/sy1*(1.0d0+bth2*                          &
                (2.0d0*anpl2+5.0d0*sy1**2-7.0d0*sy1)/sy1**2/4.0d0)
  rr(-1,1,2) = -2.0d0*bth4*anpl/sy1**2*(1.0d0+1.5d0*bth2*            &
                (anpl2-3.0d0*sy1+2.0d0*sy1**2)/sy1**2)
  rr(-1,2,2) = -2.0d0*bth4/sy1*(1.0d0+0.75d0*bth2*                   &
                (2.0d0*anpl2-3.0d0*sy1+4.0d0*sy1**2)/sy1**2)
  rr(-2,0,2) = -4.0d0*bth2/sy2*(1.0d0+bth2*                          &
                (2.0d0*anpl2+5.0d0*sy2**2-7.0d0*sy2)/sy2**2/4.0d0)
  rr(-2,1,2) = -2.0d0*bth4*anpl/sy2**2*(1.0d0+1.5d0*bth2*            &
                (anpl2-3.0d0*sy2+2.0d0*sy2**2)/sy2**2)
  rr(-2,2,2) = -2.0d0*bth4/sy2*(1.0d0+0.75d0*bth2*                   &
                (2.0d0*anpl2-3.0d0*sy2+4.0d0*sy2**2)/sy2**2)
  if (llm.eq.2) return
  rr(0,0,3) = -(12.0d0*bth4+3.0d0*(3.0d0+2.0d0*anpl2)*bth6)
  rr(0,1,3) = -6.0d0*anpl*bth6
  rr(0,2,3) = -6.0d0*bth6
  rr(-1,0,3) = -12.0d0*bth4/sy1*                                     &
               (1.0d0+bth2*(3.0d0+0.5d0*anpl2/sy1**2-2.25d0/sy1))
  rr(-1,1,3) = -6.0d0*anpl*bth6/sy1**2*                              &
               (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy1**2-5.5d0/sy1))
  rr(-1,2,3) = -6.0d0*bth6/sy1*                                      &
               (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy1**2-2.75d0/sy1))
  rr(-2,0,3) = -12.0d0*bth4/sy2*                                     &
                (1.0d0+bth2*(3.0d0+0.5d0*anpl2/sy2**2-2.25d0/sy2))
  rr(-2,1,3) = -6.0d0*anpl*bth6/sy2**2*                              &
                (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy2**2-5.5d0/sy2))
  rr(-2,2,3) = -6.0d0*bth6/sy2*                                      &
                (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy2**2-2.75d0/sy2))
  rr(-3,0,3) = -12.0d0*bth4/sy3*                                     &
                (1.0d0+bth2*(3.0d0+0.5d0*anpl2/sy3**2-2.25d0/sy3))
  rr(-3,1,3) = -6.0d0*anpl*bth6/sy3**2*                              &
                (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy3**2-5.5d0/sy3))
  rr(-3,2,3) = -6.0d0*bth6/sy3*                                      &
                (1.0d0+bth2*(5.25d0+1.5d0*anpl2/sy3**2-2.75d0/sy3))
end subroutine hermitian

subroutine antihermitian(ri,lrm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ANTIHERMITIAN(r,lrm)                                                 C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use constants, only : pi_sqrt
  implicit none
  integer      , intent(in)  :: lrm
  real(kind=r8), intent(out) :: ri(lrm,0:2,lrm)
  integer :: n, k, m, mm
  integer      , parameter :: lmx = 10, nmx = lmx+2
  real(kind=r8) fsbi(nmx)
  real(kind=r8) yg,amu,anpl,cmu,anpl2,dnl,cr,ci,ygn,rdu2,rdu
  real(kind=r8) du,ub,aa,up,um,gp,gm,xp,xm,eep,eem,ee,cm,cim
  real(kind=r8) fi0p0,fi1p0,fi2p0,fi0m0,fi1m0,fi2m0
  real(kind=r8) fi0p1,fi1p1,fi2p1,fi0m1,fi1m1,fi2m1,fi0m,fi1m,fi2m

  common/muloc/amu
  common/yloc/yg
  common/npllc/anpl,anpl2
  common/uu/ub,cmu
  common/cri/cr,ci

  interface
    function fact(k)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fact
      integer, intent(in) :: k
    end function fact
    subroutine ssbi(zz,n,l,fsbi)
      use kinds, only : r8
      implicit none
      integer      , parameter   :: lmx=10, nmx = lmx+2
      integer      , intent(in)  :: n, l
      real(kind=r8), intent(in)  :: zz
      real(kind=r8), intent(out) :: fsbi(nmx)
    end subroutine ssbi
  end interface

  do n=1,lrm
    do k=0,2
      do m=1,lrm
        ri(n,k,m)=0.0d0
      end do
    end do
  end do
  dnl=1.0d0-anpl2
  cmu=anpl*amu
  do n=1,lrm
    ygn=n*yg
    rdu2=ygn**2-dnl
    if(rdu2>0.0d0) then
      rdu=sqrt(rdu2)
      du=rdu/dnl
      ub=anpl*ygn/dnl
      aa=amu*anpl*du
      if (abs(aa)>1.0d0) then
        up=ub+du
        um=ub-du
        gp=anpl*up+ygn
        gm=anpl*um+ygn
        xp=up+1.0d0/cmu
        xm=um+1.0d0/cmu
        eem=exp(-amu*(gm-1.0d0))
        eep=exp(-amu*(gp-1.0d0))
        fi0p0=-1.0d0/cmu
        fi1p0=-xp/cmu
        fi2p0=-(1.0d0/cmu**2+xp*xp)/cmu
        fi0m0=-1.0d0/cmu
        fi1m0=-xm/cmu
        fi2m0=-(1.0d0/cmu**2+xm*xm)/cmu
        do m=1,lrm
          fi0p1=-2.0d0*m*(fi1p0-ub*fi0p0)/cmu
          fi0m1=-2.0d0*m*(fi1m0-ub*fi0m0)/cmu
          fi1p1=-((1.0d0+2.0d0*m)*fi2p0-2.0d0*(m+1.0d0)*ub*fi1p0    &
                 +up*um*fi0p0)/cmu
          fi1m1=-((1.0d0+2.0d0*m)*fi2m0-2.0d0*(m+1.0d0)*ub*fi1m0    &
                 +up*um*fi0m0)/cmu
          fi2p1=(2.0d0*(1.0d0+m)*fi1p1-2.0d0*m*                     &
               (ub*fi2p0-up*um*fi1p0))/cmu
          fi2m1=(2.0d0*(1.0d0+m)*fi1m1-2.0d0*m*                     &
               (ub*fi2m0-up*um*fi1m0))/cmu
          if(m>=n) then
            ri(n,0,m)=0.5d0*ci*dnl**m*(fi0p1*eep-fi0m1*eem)
            ri(n,1,m)=0.5d0*ci*dnl**m*(fi1p1*eep-fi1m1*eem)
            ri(n,2,m)=0.5d0*ci*dnl**m*(fi2p1*eep-fi2m1*eem)
          end if
          fi0p0=fi0p1
          fi1p0=fi1p1
          fi2p0=fi2p1
          fi0m0=fi0m1
          fi1m0=fi1m1
          fi2m0=fi2m1
        end do
      else
        ee=exp(-amu*(ygn-1.0d0+anpl*ub))
        call ssbi(aa,n,lrm,fsbi)
        do m=n,lrm
          cm=pi_sqrt*fact(m)*du**(2*m+1)
          cim=0.5d0*ci*dnl**m
          mm=m-n+1
          fi0m=cm*fsbi(mm)
          fi1m=-0.5d0*aa*cm*fsbi(mm+1)
          fi2m=0.5d0*cm*(fsbi(mm+1)+0.5d0*aa*aa*fsbi(mm+2))
          ri(n,0,m)=cim*ee*fi0m
          ri(n,1,m)=cim*ee*(du*fi1m+ub*fi0m)
          ri(n,2,m)=cim*ee*(du*du*fi2m+2.0d0*du*ub*fi1m+ub*ub*fi0m)
        end do
      end if
    end if
  end do
end subroutine antihermitian

subroutine ssbi(zz,n,l,fsbi)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SSBI(zz,n,l,fsbi)                                                    C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  integer      , parameter   :: lmx=10, nmx = lmx+2
  integer      , intent(in)  :: n, l
  real(kind=r8), intent(in)  :: zz
  real(kind=r8), intent(out) :: fsbi(nmx)
  integer       :: k, m, mm
  real(kind=r8) :: c0, c1, sbi
  real(kind=r8), parameter   :: eps = 1.0d-10

  interface
    function gammln(xx)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: gammln
      real(kind=r8), intent(in) :: xx
    end function gammln
  end interface

  do m=n,l+2
    c0=1.0d0/exp(gammln(dble(m)+1.5d0))
    sbi=c0
    do k=1,50
      c1=c0*0.25d0*zz**2/(dble(m+k)+0.5d0)/dble(k)
      sbi=sbi+c1
      if(c1/sbi<eps) exit
      c0=c1
    end do
    mm=m-n+1
    fsbi(mm)=sbi
  end do
end subroutine ssbi

subroutine diel_tens_wr(ce330,cepsl,lrm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! DIEL_TENS_WR(ce330,cepsl,lrm) weakly relativistic dielectric tensor  C
!         Krivenki and Orefice, JPP 30,125 (1983)                      C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  integer         , intent(in)  :: lrm
  complex(kind=r8), intent(out) :: ce330, cepsl(3,3,lrm)
  integer l,lm,is,k
  real(kind=r8) xg,yg,anpl,anpl2,amu,fcl,asl,bsl,w
  complex(kind=r8) ui,ca11,ca12,ca13,ca22,ca23,ca33
  complex(kind=r8) cq0p,cq0m,cq1p,cq1m,cq2p
  complex(kind=r8) cefp(0:lrm,0:2),cefm(0:lrm,0:2)

  common/muloc/amu
  common/xloc/xg
  common/yloc/yg
  common/npllc/anpl,anpl2

  interface
    function fact(k)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fact
      integer, intent(in) :: k
    end function fact
    subroutine fsup(cefp,cefm,lrm)
      use kinds, only : r8
      implicit none
      integer         , intent(in)  :: lrm
      complex(kind=r8), intent(out) :: cefp(0:lrm,0:2),cefm(0:lrm,0:2)
    end subroutine fsup
  end interface

  ui=(0.0d0,1.0d0)
  call fsup(cefp,cefm,lrm)
  do l=1,lrm
    lm=l-1
    fcl=0.5d0**l*((1.0d0/yg)**2/amu)**lm*fact(2*l)/fact(l)
    ca11=(0.d0,0.d0)
    ca12=(0.d0,0.d0)
    ca13=(0.d0,0.d0)
    ca22=(0.d0,0.d0)
    ca23=(0.d0,0.d0)
    ca33=(0.d0,0.d0)
    do is=0,l
      k=l-is
      w=(-1.0d0)**k
      asl=w/(fact(is+l)*fact(l-is))
      bsl=asl*(is*is+dble(2*k*lm*(l+is))/(2.0d0*l-1.0d0))
      cq0p=amu*cefp(is,0)
      cq0m=amu*cefm(is,0)
      cq1p=amu*anpl*(cefp(is,0)-cefp(is,1))
      cq1m=amu*anpl*(cefm(is,0)-cefm(is,1))
      cq2p=cefp(is,1)+amu*anpl2*(cefp(is,2)+cefp(is,0)-2.0d0*cefp(is,1))
      ca11=ca11+is**2*asl*cq0p
      ca12=ca12+is*l*asl*cq0m
      ca22=ca22+bsl*cq0p
      ca13=ca13+is*asl*cq1m/yg
      ca23=ca23+l*asl*cq1p/yg
      ca33=ca33+asl*cq2p/yg**2
    end do
    cepsl(1,1,l) =  - xg*ca11*fcl
    cepsl(1,2,l) =  + ui*xg*ca12*fcl
    cepsl(2,2,l) =  - xg*ca22*fcl
    cepsl(1,3,l) =  - xg*ca13*fcl
    cepsl(2,3,l) =  - ui*xg*ca23*fcl
    cepsl(3,3,l) =  - xg*ca33*fcl
  end do
  cq2p=cefp(0,1)+amu*anpl2*(cefp(0,2)+cefp(0,0)-2.0d0*cefp(0,1))
  ce330=1.0d0-xg*amu*cq2p
  cepsl(1,1,1) = 1.d0 + cepsl(1,1,1)
  cepsl(2,2,1) = 1.d0 + cepsl(2,2,1)
  do l=1,lrm
    cepsl(2,1,l) = - cepsl(1,2,l)
    cepsl(3,1,l) =   cepsl(1,3,l)
    cepsl(3,2,l) = - cepsl(2,3,l)
  end do
end subroutine diel_tens_wr

subroutine fsup(cefp,cefm,lrm)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! FSUP(cefp,cefm,lrm)                                                  C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  integer         , intent(in)  :: lrm
  complex(kind=r8), intent(out) :: cefp(0:lrm,0:2),cefm(0:lrm,0:2)
  real(kind=r8), parameter :: soglia = 0.7d0
  integer l,is,iq,ir,iflag
  real(kind=r8) yg,anpl,anpl2,amu,alpha,phi2,phim,psi,apsi
  real(kind=r8) xp,yp,xm,ym,x0,y0,zrp,zip,zrm,zim,zr0,zi0
  complex(kind=r8) ui,czp,czm,cf12,cf32,cphi,cz0,cdz0,cf0,cf1,cf2

  common/muloc/amu
  common/yloc/yg
  common/npllc/anpl,anpl2

  interface
    subroutine zetac(xi, yi, zr, zi, iflag)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: xi, yi
      real(kind=r8), intent(out) :: zr, zi
      integer      , intent(out) :: iflag
    end subroutine zetac
  end interface

  ui=(0.0d0,1.0d0)
  psi=sqrt(0.5d0*amu)*anpl
  apsi=abs(psi)
  do is=0,lrm
    alpha=anpl2/2.0d0+is*yg-1.0d0
    phi2=amu*alpha
    phim=sqrt(abs(phi2))
    if (alpha>=0) then
      xp=psi-phim
      yp=0.0d0
      xm=-psi-phim
      ym=0.0d0
      x0=-phim
      y0=0.0d0
    else
      xp=psi
      yp=phim
      xm=-psi
      ym=phim
      x0=0.0d0
      y0=phim
    end if
    call zetac (xp,yp,zrp,zip,iflag)
    call zetac (xm,ym,zrm,zim,iflag)
    czp=cmplx(zrp,zip)
    czm=cmplx(zrm,zim)
    if (alpha>0) then
      cf12=-(czp+czm)/(2.0d0*phim)
    else if (alpha<0) then
      cf12=-ui*(czp+czm)/(2.0d0*phim)
    else
      cf12=(0.0d0,0.0d0)
    end if
    if(apsi>soglia) then
      cf32=-(czp-czm)/(2.0d0*psi)
    else
      cphi=phim
      if(alpha<0) cphi=-ui*phim
      call zetac (x0,y0,zr0,zi0,iflag)
      cz0=cmplx(zr0,zi0)
      cdz0=2.0d0*(1.0d0-cphi*cz0)
      cf32=cdz0
    end if
    cf0=cf12
    cf1=cf32
    cefp(is,0)=cf32
    cefm(is,0)=cf32
    do l=1,is+2
      iq=l-1
      if(apsi>soglia) then
        cf2=(1.0d0+phi2*cf0-(iq+0.5d0)*cf1)/psi**2
      else
        cf2=(1.0d0+phi2*cf1)/dble(iq+1.5d0)
      end if
      ir=l-is
      if(ir>=0) then
        cefp(is,ir)=cf2
        cefm(is,ir)=cf2
      end if
      cf0=cf1
      cf1=cf2
    end do
    if (is/=0) then
      alpha=anpl2/2.0d0-is*yg-1.0d0
      phi2=amu*alpha
      phim=sqrt(abs(phi2))
      if (alpha>=0) then
        xp=psi-phim
        yp=0.0d0
        xm=-psi-phim
        ym=0.0d0
        x0=-phim
        y0=0.0d0
      else
        xp=psi
        yp=phim
        xm=-psi
        ym=phim
        x0=0.0d0
        y0=phim
      end if
      call zetac (xp,yp,zrp,zip,iflag)
      call zetac (xm,ym,zrm,zim,iflag)
      czp=cmplx(zrp,zip)
      czm=cmplx(zrm,zim)
      if (alpha>0) then
        cf12=-(czp+czm)/(2.0d0*phim)
      else if (alpha<0) then
        cf12=-ui*(czp+czm)/(2.0d0*phim)
      else
        cf12=(0.0d0,0.0d0)
      end if
      if(apsi>soglia) then
        cf32=-(czp-czm)/(2.0d0*psi)
      else
        cphi=phim
        if(alpha<0) cphi=-ui*phim
        call zetac (x0,y0,zr0,zi0,iflag)
        cz0=cmplx(zr0,zi0)
        cdz0=2.0d0*(1.0d0-cphi*cz0)
        cf32=cdz0
      end if
      cf0=cf12
      cf1=cf32
      do l=1,is+2
        iq=l-1
        if(apsi>soglia) then
          cf2=(1.0d0+phi2*cf0-(iq+0.5d0)*cf1)/psi**2
        else
          cf2=(1.0d0+phi2*cf1)/dble(iq+1.5d0)
        end if
        ir=l-is
        if(ir>=0) then
          cefp(is,ir)=cefp(is,ir)+cf2
          cefm(is,ir)=cefm(is,ir)-cf2
        end if
        cf0=cf1
        cf1=cf2
      end do
    end if
  end do
end subroutine fsup

subroutine cutoff
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CUTOFF                                                               C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nnw, nnh
  implicit none
  integer nr,nz,nppn,nppt,k,i
  character(len=255) idstring,prefix
  real(kind=r8) rv(nnw),zv(nnh),psin(nnw,nnh)
  real(kind=r8) psinv,dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  real(kind=r8) dens19,ddenspsin,psdbnd,rpsim,zpsim,r0m,z0m,anpl0
  real(kind=r8) x00,y00,z00,anx0c,any0c,anz0c,bphi,brr,bzz,btot,tte
  real(kind=r8) fpe,fr,fl,fce,fce2,fce3,fce4
  real(kind=r8) frpl,flpl,fcepl,fcepl2,fcepl3,fcepl4,anpl,anpl2,dnl

  common/cpsin/rv,zv,psin
  common/eqnn/nr,nz,nppn,nppt
  common/psival/psinv
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/dens/dens19,ddenspsin
  common/anic/anx0c,any0c,anz0c
  common/mirr/x00,y00,z00
  common/densbnd/psdbnd
  common/spectinfo/idstring,prefix

  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
    function temp(arg)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: temp
      real(kind=r8), intent(in)  :: arg
    end function temp
    subroutine density(arg)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: arg
    end subroutine density
  end interface

  open(14,file=trim(prefix)//'freq'//trim(idstring)//'.txt',action='write')
  open(16,file=trim(prefix)//'profiles'//trim(idstring)//'.txt',action='write')
  write(14,*) nr
  write(14,*)' #R psi fpe fL fR fce 2fec 3fce 4fce'//               &
             ' fLa fRa fcea 2fcea 3fcea 4fcea npl'
  write(16,*) ' #R z psi Te B'
  r0m=sqrt(x00**2+y00**2)*1.0d-2
  anpl0=any0c
  z0m=1.0d-2*z00
  k=locate(zv(1:nnh),z0m)
  zpsim=zv(k)
  if(abs(z0m-zv(k))>abs(z0m-zv(k+1))) zpsim=zv(k+1)
  do i=1,nr
    rpsim=rv(i)
    call equinum(rpsim,zpsim)
    bphi=fpolv/rpsim
    brr=-dpsidz/rpsim
    bzz= dpsidr/rpsim
    btot=sqrt(bphi**2+brr**2+bzz**2)
    fce=28.0d0*Btot
    if(psinv>0.and.psinv<psdbnd) then
      tte=temp(psinv)
      write(16,1) rpsim,zpsim,psinv,tte,btot
      call density(psinv)
      fpe=8.98d+0*sqrt(dens19*10.0d0)
      fr=0.5d0*(fce+sqrt(fce**2+4.0d0*fpe**2))
      fl=0.5d0*(-fce+sqrt(fce**2+4.0d0*fpe**2))
    else
      fpe=0.0d0
      fr=fce
      fl=0.0d0
    end if
    fce2=2.0d0*fce
    fce3=3.0d0*fce
    fce4=4.0d0*fce
    anpl=anpl0*r0m/rpsim
    anpl2=anpl*anpl
    fcepl=fce
    frpl=fcepl
    flpl=0.0d0
    if(anpl2<1.0d0) then
      dnl=1.0d0-anpl2
      fcepl=fce/sqrt(dnl)
      frpl=0.5d0*(fce+sqrt(fce**2+4.0d0*fpe**2/dnl))
      flpl=0.5d0*(-fce+sqrt(fce**2+4.0d0*fpe**2/dnl))
      fcepl2=2.0d0*fcepl
      fcepl3=3.0d0*fcepl
      fcepl4=4.0d0*fcepl
    end if
    write(14,1) rpsim,psinv,fpe,fl,fr,fce,fce2,fce3,fce4, &
                flpl,frpl,fcepl,fcepl2,fcepl3,fcepl4,anpl
  end do
  close(14)
  close(16)
  return
  1 format(20(1x,e12.5))
end subroutine cutoff

subroutine rayenddata
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! RAYENDDATA                                                           C
!   uses                                                               C
!   sets                                                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx
  implicit none
  integer iif,iox
  real(kind=r8) rrr(nmx)
  real(kind=r8) zzv(nmx),stv(nmx),phiv(nmx),xgv(nmx),anprv(nmx)
  real(kind=r8) rfin(2),zfin(2)

  common/iiv/iif
  common/imode/iox
  common/rrr/rrr
  common/traj/stv,zzv,phiv,xgv,anprv
  common/rayend/rfin,zfin

  rfin(iox)=rrr(iif)*1.d-2
  zfin(iox)=zzv(iif)
  return
end subroutine rayenddata

subroutine raysections
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! RAYSECTIONS                                                          C
!   uses                                                               C
!   sets                                                               C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx, nscmx
  implicit none
  integer istart,istep,iif,nsect,iistep
  integer lsect(nscmx)
  real(kind=r8) psjki(nmx),psisect(nmx,nscmx)
  real(kind=r8) ttasect(nmx,nscmx),tta(nmx)

  common/tta/tta
  common/psjki/psjki
  common/raysec/lsect,nsect
  common/pstasc/psisect,ttasect
  common/iiv/iif

  nsect=1
  istart=1
  psisect(1,1)=psjki(1)
  ttasect(1,1)=tta(1)
  do istep=2,iif-1
    iistep=istep-istart+1
    psisect(iistep,nsect)=psjki(istep)
    ttasect(iistep,nsect)=tta(istep)
    if( psjki(istep)*psjki(istep-1) < 0.0d0 .or.                 &
         ( (abs(psjki(istep))-abs(psjki(istep-1))) *                &
           (abs(psjki(istep+1))-abs(psjki(istep))) < 0.0d0 )     &
        ) then
      if (psjki(istep)*psjki(istep-1)<0.0d0) then
        psisect(iistep,nsect)=-psisect(iistep,nsect)
      end if
      if (nsect<nscmx) then
        lsect(nsect)=iistep
        nsect=nsect+1
        istart=istep
        psisect(1,nsect)=psjki(istep)
        ttasect(1,nsect)=tta(istep)
      else
        write(*,*) 'step',istep,'/',iif,': max number of ray'//     &
                   'sections exceeded'
      end if
    end if
  end do
  iistep=iif-istart+1
  psisect(iistep,nsect)=psjki(istep)
  ttasect(iistep,nsect)=tta(istep)
  if (psjki(iif)*psjki(iif-1)<0.0d0) then
    psisect(iistep,nsect)=-psisect(iistep,nsect)
  end if
  lsect(nsect)=iistep
end subroutine raysections

integer function islhfs(r,z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ISLHFS(r,z)                                                          C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nhlmx
  implicit none
  real(kind=r8), intent(in) :: r, z
  integer k,nhl
  real(kind=r8) hlr(nhlmx),hlz(nhlmx),chl(nhlmx,4),delz,rhl

  common/hline/chl,hlr,hlz,nhl

  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    function spli(cspli,k,dx)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: spli
      integer      , intent(in) :: k
      real(kind=r8), intent(in) :: cspli(:,:), dx
    end function spli
  end interface

  if(z<hlz(nhl).and.z>hlz(1)) then
    k=locate(hlz(1:nhl),z)
    if(k==0) k=1
    if(k==nhl) k=nhl-1
    delz=z-hlz(k)
    rhl=spli(chl,k,delz)
  else if (z>=hlz(nhl)) then
!       z > separatrix upper limit!
    rhl=hlr(nhl)
  else
!       BELOW X POINT!
    rhl=hlr(1)
  end if
  if (r>=rhl) then
    islhfs=1
  else
    islhfs=-1
  end if
end function islhfs
      
subroutine taonaxis(j,k,iox)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TAONAXIS(j,k,iox)                                                    C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : jmx, kmx, nmx, nscmx, naxmx
  implicit none
  integer, intent(in) :: j, k, iox
  integer i,ist,isect,nsect,nax
  
  integer iinf,isup,lsect(nscmx),iend
  real(kind=r8) psiax(naxmx),rax(naxmx),psisect(nmx,nscmx),psisjk(nmx)
  real(kind=r8) ttasect(nmx,nscmx),ttaax(naxmx,jmx,kmx,2)
  real(kind=r8) ttainf,ttasup,dtta

  common/raysec/lsect,nsect
  common/pstasc/psisect,ttasect
  common/axis/psiax,rax,nax
  common/taax/ttaax

  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    subroutine intlin(x1,y1,x2,y2,x,y)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: x1, y1, x2, y2, x
      real(kind=r8), intent(out) :: y
    end subroutine intlin
  end interface

  ttaax(1,j,k,iox)=0.0d0
  do i=2,nax
    ttaax(i,j,k,iox)=ttaax(i-1,j,k,iox)
    iend=1
    do isect=1,nsect
      iend=iend+lsect(isect)-1
      do ist=1,lsect(isect)
        psisjk(ist)=psisect(ist,isect)
      end do
      iinf=locate(psisjk(1:lsect(isect)),psiax(i-1))
      if (iinf>=lsect(isect)) then
        ttainf=ttasect(lsect(isect),isect)
      else if (iinf<1) then
        ttainf=ttasect(1,isect)
      else
        call intlin(psisjk(iinf),ttasect(iinf,isect), &
              psisjk(iinf+1),ttasect(iinf+1,isect),   &
              psiax(i-1),ttainf)
      end if
      isup=locate(psisjk(1:lsect(isect)),psiax(i))
      if (isup>=lsect(isect)) then
        ttasup=ttasect(lsect(isect),isect)
      else if (isup<1) then
        ttasup=ttasect(1,isect)
      else
        call intlin(psisjk(isup),ttasect(isup,isect), &
              psisjk(isup+1),ttasect(isup+1,isect),   &
              psiax(i),ttasup)
      end if
      if (psisjk(1)<psisjk(lsect(isect))) then
        dtta=ttainf-ttasup
      else
        dtta=ttasup-ttainf
      end if
      ttaax(i,j,k,iox)=ttaax(i,j,k,iox)+dtta
    end do
  end do
end subroutine taonaxis

subroutine sub_trad(j,k)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUB_TRAD(j,k)                                                        C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use phys_constants, only : cai_mev
  use spece_data, only : nmx
  implicit none
  integer, intent(in) :: j, k
  integer iovlp,iharm,iif,iox,ifilp,kk,ifin,i,l !,nharm
  real(kind=r8) tta(nmx),tradv(2),psinv !,psjki(nmx)
  real(kind=r8) tauvm(2),ygv(nmx),rrr(nmx) !,rmeanv(2),drmeanv(2),rresv(2)
  real(kind=r8) st,fghz,taumn,taumx,rmean,drmean,rres,tradc,ygr,cint
  real(kind=r8) ttar1,ttar2,ggr,r50,rinf,rsup,z50,zinf,zsup,t50,tinf,tsup
  real(kind=r8) taures(5),tres(4),tautemp,ttemp,ps50,psinf,pssup
  real(kind=r8) zzv(nmx)
  real(kind=r8) stv(nmx),phiv(nmx),xgv(nmx),anprv(nmx)
  real(kind=r8) tauv(nmx),alphav(nmx),tev(nmx)
  logical alloutput,goutput

  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    subroutine intlin(x1,y1,x2,y2,x,y)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: x1, y1, x2, y2, x
      real(kind=r8), intent(out) :: y
    end subroutine intlin
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface

  common/traj/stv,zzv,phiv,xgv,anprv
  common/ss/st
  common/tta/tta
  common/ttar/ttar1,ttar2
!  common/psjki/psjki
  common/fghz/fghz
  common/taumnx/taumn,taumx
!      common/nharm/nharm
  common/harm/iovlp,iharm
  common/iiv/iif
  common/imode/iox
  common/tradres/tradv,tauvm
  common/ygv/ygv
  common/rrr/rrr
!      common/rad/rmeanv,drmeanv,rresv
  common/atjki/tauv,alphav,tev
  common/out/alloutput,goutput

  common/psival/psinv

  ifin=iif
  tradv(iox)=tta(ifin)
  tradc=tta(ifin)
  tauvm(iox)=taumn
  if(tradc>0.0d0) then
    rmean=ttar1/tradc
    drmean=sqrt(ttar2/tradc-rmean**2)
! ******* effective radial resolution (see Bartlett EC-10 1997) *******
    if(j==1.and.k==1) then
      ggr=0.025d0*tradc
      kk=locate(tta(1:ifin),ggr)
      call intlin(tta(kk),rrr(kk),tta(kk+1),rrr(kk+1),ggr,rinf)
      call intlin(tta(kk),zzv(kk),tta(kk+1),zzv(kk+1),ggr,zinf)
      call intlin(tta(kk),tev(kk),tta(kk+1),tev(kk+1),ggr,tinf)
      call equinum(0.01_r8*rinf,zinf)
      psinf=psinv
!      call intlin(tta(kk),psjki(kk),tta(kk+1),psjki(kk+1),ggr,psinf)
      ggr=0.5d0*tradc
      kk=locate(tta(1:ifin),ggr)
      call intlin(tta(kk),rrr(kk),tta(kk+1),rrr(kk+1),ggr,r50)
      call intlin(tta(kk),zzv(kk),tta(kk+1),zzv(kk+1),ggr,z50)
      call intlin(tta(kk),tev(kk),tta(kk+1),tev(kk+1),ggr,t50)
      call equinum(0.01_r8*r50,z50)
      ps50=psinv
!      call intlin(tta(kk),psjki(kk),tta(kk+1),psjki(kk+1),ggr,ps50)
      ggr=0.975d0*tradc
      kk=locate(tta(1:ifin),ggr)
      call intlin(tta(kk),rrr(kk),tta(kk+1),rrr(kk+1),ggr,rsup)
      call intlin(tta(kk),zzv(kk),tta(kk+1),zzv(kk+1),ggr,zsup)
      call intlin(tta(kk),tev(kk),tta(kk+1),tev(kk+1),ggr,tsup)
      call equinum(0.01_r8*rsup,zsup)
      pssup=psinv
!      call intlin(tta(kk),psjki(kk),tta(kk+1),psjki(kk+1),ggr,pssup)
      zinf=zinf*1.d2
      z50=z50*1.d2
      zsup=zsup*1.d2
    end if
! *********************************************************************
!ccccccc
    do i=4,1,-1
      ygr=1.0d0/dble(i)
      kk=locate(ygv(1:ifin),ygr)
      if (kk>=ifin) then
        tautemp=taumn
        ttemp=0.0d0
      else if (kk<=0) then
        tautemp=0.0d0
        ttemp=0.0d0
      else
        call intlin(ygv(kk),tauv(kk),ygv(kk+1),tauv(kk+1),ygr,tautemp)
        call intlin(ygv(kk),tev(kk),ygv(kk+1),tev(kk+1),ygr,ttemp)
      end if
      taures(i+1)=tautemp
      tres(i)=ttemp
      do l=i+2,5
        taures(i+1)=taures(i+1)-taures(l)
      end do
    end do
    taures(1)=taumn-taures(2)-taures(3)-taures(4)-taures(5)
!ccccccc
    ygr=1.0d0/dble(iharm)
    kk=locate(ygv(1:ifin),ygr)
    if (kk>=ifin.or.kk<=0) then
      rres=0.0d0
    else
      call intlin(ygv(kk),rrr(kk),ygv(kk+1),rrr(kk+1),ygr,rres)
    end if
  else
    rmean=0.0d0
    drmean=0.0d0
    rres=0.0d0
    taures(1)=0.0d0
    do i=1,4
      taures(i+1)=0.0d0
      tres(i)=0.0d0
    end do
    if(j==1.and.k==1) then
      rinf=0.0d0
      rsup=0.0d0
      r50=0.0d0
      zinf=0.0d0
      zsup=0.0d0
      z50=0.0d0
      tinf=0.0d0
      tsup=0.0d0
      t50=0.0d0
      psinf=5.0d0
      pssup=5.0d0
      ps50=5.0d0
    end if
  end if
!     cint measured in 10^6 eV steradian^_1 m^-2
  cint=cai_mev*fghz**2*tradc
  if(alloutput.and.j==1.and.k==1) then
    ifilp=10+iox
    write(ifilp,99) fghz,st,tradc,cint,rmean,drmean,taumn,            &
            dble(iovlp),dble(iharm),rres,r50,rinf,rsup,z50,zinf,zsup, &
            ps50,psinf,pssup,t50,tinf,tsup,                           &
            tres(2),tres(3),tres(4),taures(2),taures(3),taures(4)
  end if
!      rmeanv(iox)=rmean
!      drmeanv(iox)=drmean
!      rresv(iox)=rres
  99 format(30(1x,e16.8e3))
end subroutine sub_trad

subroutine trad_back
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TRAD_BACK                                                            C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx
  implicit none
  integer i,iif,iox
  real(kind=r8) dst,gg0,gg1,tta,tau,tau0
  real(kind=r8) trbkv(2),tauv(nmx),alphav(nmx),tev(nmx)

  common/dsds/dst
  common/tradbck/trbkv
  common/atjki/tauv,alphav,tev
  common/iiv/iif
  common/imode/iox

  tta=0.0d0
  tau0=0.0d0
  do i=iif,2,-1
    tau=tau0+(alphav(i)+alphav(i-1))*dst/2.0d0
    gg1=alphav(i-1)*tev(i-1)*exp(-tau)
    gg0=alphav(i)*tev(i)*exp(-tau0)
    tta=tta+(gg0+gg1)*dst/2.0d0
    tau0=tau
  end do
  trbkv(iox)=tta
end subroutine trad_back

subroutine trad_tot(j,k)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TRAD_TOT(j,k)                                                        C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use phys_constants, only : cai_mev
  use spece_data, only : jmx, kmx
  implicit none
  integer, intent(in) :: j, k
  real(kind=r8) tradv(2),tauvm(2),trbkv(2) !,rmeanv(2),drmeanv(2),rresv(2)
  real(kind=r8) ttradtxx(jmx,kmx),ttradtyy(jmx,kmx)
  real(kind=r8) citxx(jmx,kmx),cityy(jmx,kmx)
  real(kind=r8) fghz,r,p,trom,trxm,trob,trxb,pto,ptx,pto2,ptx2,coffpol
  real(kind=r8) p1,p1p,pp,r2,dd,trqom,trqxm
  real(kind=r8) rfin(2),zfin(2)
  real(kind=r8) x00,y00,z00,rstart,zstart
  integer co,lhfsst
  logical alloutput,goutput

  common/fghz/fghz
  common/mirr/x00,y00,z00
  common/rp/r,p
  common/tradres/tradv,tauvm
  common/tradbck/trbkv
  common/polres/coffpol
!      common/rad/rmeanv,drmeanv,rresv
  common/trad/ttradtxx,ttradtyy
  common/intens/citxx,cityy
  common/rayend/rfin,zfin
  common/out/alloutput,goutput

  interface
    integer function islhfs(r,z)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: r, z
    end function islhfs
  end interface

  trom=tradv(1)
  trxm=tradv(2)
  trob=trbkv(1)
  trxb=trbkv(2)
!     Check valid only for almost horizontal antenna
  rstart=sqrt(x00**2+y00**2)*1.d-2
  zstart=z00*1.d-2
  lhfsst=islhfs(rstart,zstart)
  if (islhfs(rfin(1),zfin(1))==lhfsst .or.    &
      islhfs(rfin(2),zfin(2))==lhfsst) then
    co=1
!        pp=0.0d0
  else
    co=0
!        pp=p
  end if
  pp=p
  pto=exp(-tauvm(1))
  ptx=exp(-tauvm(2))
  pto2=pto*pto
  ptx2=ptx*ptx
  p1=1.0d0-pp
  p1p=p1-pp
  r2=r*r
  if (co==1) then
    dd=1.0d0-r*p1*(ptx+pto)+r2*pto*ptx*p1p
    trqom=((1.0d0-r*p1*ptx)*trom+r*pp*pto*trxm)/dd
    trqxm=((1.0d0-r*p1*pto)*trxm+r*pp*ptx*trom)/dd
  else
!     quasi-O,-X temperatures considering Trad_in/=Trad_out
    dd=1.0d0+r2*(r2*p1p*p1p*pto2*ptx2                               &
                 -p1*p1*(pto2+ptx2)-2.0d0*pp*pp*pto*ptx)
    trqom=((1.0d0-r2*ptx*(pp*pp*pto+p1*p1*ptx))*trom                &
          +r*pto*(p1*(1.0d0-r2*p1p*ptx2)*trob                       &
                 +r*pp*p1*(pto+ptx)*trxm                            &
                 +pp*(1.0d0+r2*p1p*pto*ptx)*trxb))/dd
    trqxm=((1.0d0-r2*pto*(pp*pp*ptx+p1*p1*pto))*trxm                &
          +r*ptx*(p1*(1.0d0-r2*p1p*pto2)*trxb                       &
                 +r*pp*p1*(pto+ptx)*trom                            &
                 +pp*(1.0d0+r2*p1p*pto*ptx)*trob))/dd
  end if
!     citxx,yy measured in 10^6 eV steradian^_1 m^-2
!     to obtain spectrum in pW Hz^-1 m^-2 st^-1 (see costley)
!     one should multiply cai_mev by 0.16022
  ttradtxx(j,k)=trqom*coffpol+trqxm*(1-coffpol)
  ttradtyy(j,k)=trqom*(1-coffpol)+trqxm*coffpol
  citxx(j,k)=fghz**2*cai_mev*ttradtxx(j,k)
  cityy(j,k)=fghz**2*cai_mev*ttradtyy(j,k)
  if(alloutput.and.j==1.and.k==1) then
    write(25,99) fghz,ttradtxx(j,k),ttradtyy(j,k)                   &
     ,citxx(j,k),cityy(j,k),trqom,trqxm,tradv(1),tradv(2)           &
     ,tauvm(1),tauvm(2),trbkv(1),trbkv(2),coffpol,dble(co)
  end if
  99 format(20(1x,e16.8e3))
end subroutine trad_tot

subroutine trad_sum
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TRAD_SUM                                                             C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : jmx, kmx
  implicit none
  integer :: nray, ktx, j, k
  real(kind=r8) q(jmx)
  real(kind=r8) tradxtot(jmx,kmx),tradytot(jmx,kmx)
  real(kind=r8) citxtot(jmx,kmx),citytot(jmx,kmx)
  real(kind=r8) ttradtxx(jmx,kmx),ttradtyy(jmx,kmx)
  real(kind=r8) citxx(jmx,kmx),cityy(jmx,kmx)
  real(kind=r8) sumcitx,sumcity,sumtradx,sumtrady,fghz

  common/qw/q
  common/nrayktx/nray,ktx
  common/trad/ttradtxx,ttradtyy
  common/intens/citxx,cityy
  common/fghz/fghz

  tradxtot(1,1)=q(1)*ttradtxx(1,1)
  tradytot(1,1)=q(1)*ttradtyy(1,1)
  citxtot(1,1)=q(1)*citxx(1,1)
  citytot(1,1)=q(1)*cityy(1,1)
  sumtradx=tradxtot(1,1)
  sumtrady=tradytot(1,1)
  sumcitx=citxtot(1,1)
  sumcity=citytot(1,1)
  if(nray>1) then
    do j=2,nray
      do k=1,ktx
        tradxtot(j,k)=q(j)*ttradtxx(j,k)
        tradytot(j,k)=q(j)*ttradtyy(j,k)
        citxtot(j,k)=q(j)*citxx(j,k)
        citytot(j,k)=q(j)*cityy(j,k)
        sumtradx=sumtradx+tradxtot(j,k)
        sumtrady=sumtrady+tradytot(j,k)
        sumcitx=sumcitx+citxtot(j,k)
        sumcity=sumcity+citytot(j,k)
      end do
    end do
  end if
  write(15,99) fghz,sumtradx,sumtrady,sumcitx,sumcity
  99 format(20(1x,e16.8e3))
end subroutine trad_sum

subroutine taaxis_sum(goutput)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! TAAXIS_SUM                                                           C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : jmx, kmx, naxmx, nnw
  implicit none
  logical, intent(in) :: goutput
  integer nax,nray,ktx,i,j,k,kk
  real(kind=r8) q(jmx),psiax(naxmx),rax(naxmx),ttaax(naxmx,jmx,kmx,2)
  real(kind=r8) ggr,sumtaaxo(naxmx),sumtaaxx(naxmx),tradax
  real(kind=r8) tradaxmix,sumaxqo(naxmx),sumaxqx(naxmx)
  real(kind=r8) r25o,r50o,r975o,r25x,r50x,r975x
  real(kind=r8) r25qo,r50qo,r975qo,r25qx,r50qx,r975qx
  real(kind=r8) psi50o,psi25o,psi975o,psi50x,psi25x,psi975x
  real(kind=r8) rhot50o,rhot25o,rhot975o,rhot50x,rhot25x,rhot975x
  real(kind=r8) fghz,coffpol,drh
  integer i_sup,i_inf
  real(kind=r8) denom_der,gax_o,gax_x,gax_qo,gax_qx
  integer nr,nz,nppn,nppt,il
  real(kind=r8), dimension(nnw) :: psinr
  real(kind=r8), dimension(nnw,4) :: crhot

  common/axis/psiax,rax,nax
  common/taax/ttaax
  common/qw/q
  common/nrayktx/nray,ktx
  common/fghz/fghz
  common/polres/coffpol
  common/eqnn/nr,nz,nppn,nppt
  common/psinr/psinr
  common/coffrhot/crhot
  
  interface
    function locate(xx,x)
      use kinds, only : r8
      implicit none
      integer locate
      real(kind=r8), intent(in)  :: xx(:), x
    end function locate
    subroutine intlin(x1,y1,x2,y2,x,y)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)  :: x1, y1, x2, y2, x
      real(kind=r8), intent(out) :: y
    end subroutine intlin
    function spli(cspli,k,dx)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: spli
      integer      , intent(in) :: k
      real(kind=r8), intent(in) :: cspli(:,:), dx
    end function spli
  end interface

  if (goutput) write(26,*) ' '
  do i=1,nax
    sumtaaxo(i)=q(1)*ttaax(i,1,1,1)
    sumtaaxx(i)=q(1)*ttaax(i,1,1,2)
    if(nray>1) then
      do j=2,nray
        do k=1,ktx
          sumtaaxo(i)=sumtaaxo(i)+q(j)*ttaax(i,j,k,1)
          sumtaaxx(i)=sumtaaxx(i)+q(j)*ttaax(i,j,k,2)
        end do
      end do
    end if
    sumaxqo(i)=sumtaaxo(i)*coffpol+sumtaaxx(i)*(1.d0-coffpol)
    sumaxqx(i)=sumtaaxx(i)*coffpol+sumtaaxo(i)*(1.d0-coffpol)
    if (goutput) then
      i_sup=min(i+1,nax)
      i_inf=max(i-1,1)
      denom_der=rax(i_sup)-rax(i_inf)
      gax_o=(sumtaaxo(i_sup)-sumtaaxo(i_inf))/denom_der
      gax_x=(sumtaaxx(i_sup)-sumtaaxx(i_inf))/denom_der
      gax_qo=(sumaxqo(i_sup)-sumaxqo(i_inf))/denom_der
      gax_qx=(sumaxqx(i_sup)-sumaxqx(i_inf))/denom_der
!      write(26,99) fghz,rax(i),gax_o,gax_x,gax_qo,gax_qx
      write(26,99) fghz,rax(i),sumtaaxo(i),sumtaaxx(i),sumaxqo(i),sumaxqx(i)
    end if
  end do

  tradax=sumtaaxo(nax)
  tradaxmix=sumaxqo(nax)
  if (tradax>0.0d0) then
    ggr=0.025d0*tradax
    kk=locate(sumtaaxo(1:nax),ggr)
    call intlin(sumtaaxo(kk),rax(kk),sumtaaxo(kk+1),rax(kk+1),ggr,r25o)
    call intlin(sumtaaxo(kk),psiax(kk),sumtaaxo(kk+1),psiax(kk+1),ggr,psi25o)
    ggr=0.5d0*tradax
    kk=locate(sumtaaxo(1:nax),ggr)
    call intlin(sumtaaxo(kk),rax(kk),sumtaaxo(kk+1),rax(kk+1),ggr,r50o)
    call intlin(sumtaaxo(kk),psiax(kk),sumtaaxo(kk+1),psiax(kk+1),ggr,psi50o)
    ggr=0.975d0*tradax
    kk=locate(sumtaaxo(1:nax),ggr)
    call intlin(sumtaaxo(kk),rax(kk),sumtaaxo(kk+1),rax(kk+1),ggr,r975o)
    call intlin(sumtaaxo(kk),psiax(kk),sumtaaxo(kk+1),psiax(kk+1),ggr,psi975o)
  else
    r25o=0.0d0
    r50o=0.0d0
    r975o=0.0d0
    psi25o=0.0d0
    psi50o=0.0d0
    psi975o=0.0d0
  end if
  if (tradaxmix>0.0d0) then
    ggr=0.025d0*tradaxmix
    kk=locate(sumaxqo(1:nax),ggr)
    call intlin(sumaxqo(kk),rax(kk),sumaxqo(kk+1),rax(kk+1),ggr,r25qo)
    ggr=0.5d0*tradaxmix
    kk=locate(sumaxqo(1:nax),ggr)
    call intlin(sumaxqo(kk),rax(kk),sumaxqo(kk+1),rax(kk+1),ggr,r50qo)
    ggr=0.975d0*tradaxmix
    kk=locate(sumaxqo(1:nax),ggr)
    call intlin(sumaxqo(kk),rax(kk),sumaxqo(kk+1),rax(kk+1),ggr,r975qo)
  else
    r25qo=0.0d0
    r50qo=0.0d0
    r975qo=0.0d0
  end if

  tradax=sumtaaxx(nax)
  tradaxmix=sumaxqx(nax)
  if (tradax>0.0d0) then
    ggr=0.025d0*tradax
    kk=locate(sumtaaxx(1:nax),ggr)
    call intlin(sumtaaxx(kk),rax(kk),sumtaaxx(kk+1),rax(kk+1),ggr,r25x)
    call intlin(sumtaaxx(kk),psiax(kk),sumtaaxx(kk+1),psiax(kk+1),ggr,psi25x)
    ggr=0.5d0*tradax
    kk=locate(sumtaaxx(1:nax),ggr)
    call intlin(sumtaaxx(kk),rax(kk),sumtaaxx(kk+1),rax(kk+1),ggr,r50x)
    call intlin(sumtaaxx(kk),psiax(kk),sumtaaxx(kk+1),psiax(kk+1),ggr,psi50x)
    ggr=0.975d0*tradax
    kk=locate(sumtaaxx(1:nax),ggr)
    call intlin(sumtaaxx(kk),rax(kk),sumtaaxx(kk+1),rax(kk+1),ggr,r975x)
    call intlin(sumtaaxx(kk),psiax(kk),sumtaaxx(kk+1),psiax(kk+1),ggr,psi975x)
  else
    r25x=0.0d0
    r50x=0.0d0
    r975x=0.0d0
    psi25x=0.0d0
    psi50x=0.0d0
    psi975x=0.0d0
  end if
  if (tradaxmix>0.0d0) then
    ggr=0.025d0*tradaxmix
    kk=locate(sumaxqx(1:nax),ggr)
    call intlin(sumaxqx(kk),rax(kk),sumaxqx(kk+1),rax(kk+1),ggr,r25qx)
    ggr=0.5d0*tradaxmix
    kk=locate(sumaxqx(1:nax),ggr)
    call intlin(sumaxqx(kk),rax(kk),sumaxqx(kk+1),rax(kk+1),ggr,r50qx)
    ggr=0.975d0*tradaxmix
    kk=locate(sumaxqx(1:nax),ggr)
    call intlin(sumaxqx(kk),rax(kk),sumaxqx(kk+1),rax(kk+1),ggr,r975qx)
  else
    r25qx=0.0d0
    r50qx=0.0d0
    r975qx=0.0d0
  end if
! 
! rhotor
  il=locate(psinr(:nr),abs(psi50o))
  il=min(max(1,il),nr-1)
  drh=abs(psi50o)-psinr(il)
  rhot50o=abs(psi50o)/psi50o*spli(crhot,il,drh)

  il=locate(psinr(:nr),abs(psi25o))
  il=min(max(1,il),nr-1)
  drh=abs(psi25o)-psinr(il)
  rhot25o=abs(psi25o)/psi25o*spli(crhot,il,drh)

  il=locate(psinr(:nr),abs(psi975o))
  il=min(max(1,il),nr-1)
  drh=abs(psi975o)-psinr(il)
  rhot975o=abs(psi975o)/psi975o*spli(crhot,il,drh)

  il=locate(psinr(:nr),abs(psi50x))
  il=min(max(1,il),nr-1)
  drh=abs(psi50x)-psinr(il)
  rhot50x=abs(psi50x)/psi50x*spli(crhot,il,drh)

  il=locate(psinr(:nr),abs(psi25x))
  il=min(max(1,il),nr-1)
  drh=abs(psi25x)-psinr(il)
  rhot25x=abs(psi25x)/psi25x*spli(crhot,il,drh)

  il=locate(psinr(:nr),abs(psi975x))
  il=min(max(1,il),nr-1)
  drh=abs(psi975x)-psinr(il)
  rhot975x=abs(psi975x)/psi975x*spli(crhot,il,drh)

  write(23,99) fghz,r50o,r25o,r975o,r50qo,r25qo,r975qo, &
               psi50o,psi25o,psi975o,rhot50o,rhot25o,rhot975o
  write(24,99) fghz,r50x,r25x,r975x,r50qx,r25qx,r975qx, &
               psi50x,psi25x,psi975x,rhot50x,rhot25x,rhot975x
  99 format(20(1x,e16.8e3))
end subroutine taaxis_sum

subroutine pol_limit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! POL_LIMIT                                                            C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  implicit none
  real(kind=r8) bv(3)
  real(kind=r8) anx,any,anz,anpl,anpr,yg,xe2om !,ye2om,xe2xm,ye2xm
  real(kind=r8) an2,anxypr,sngam,csgam,csg2,sng2,ffo,ffo2 !,ffx,ffx2
  real(kind=r8) anpl2,dnl,del0,coffpol,deno !,denx

  common/anv/anx,any,anz
  common/nplr/anpl,anpl2,anpr
  common/resah/dnl
  common/ygyg/yg
  common/bb/bv
  common/polres/coffpol

  del0=sqrt(dnl**2+4.0d0*anpl2/yg**2)
  ffo=0.5d0*yg*(dnl+del0)
!  ffx=0.5d0*yg*(dnl-del0)
  an2=anx*anx+any*any+anz*anz
  anxypr=sqrt(anx*anx+any*any)*anpr
  sngam=(anz*anpl-an2*bv(3))/anxypr
  csgam=(any*bv(1)-anx*bv(2))*sqrt(an2)/anxypr
  csg2=csgam**2
  sng2=sngam**2
  ffo2=ffo*ffo
!  ffx2=ffx*ffx
  deno=ffo2+anpl2
!  denx=ffx2+anpl2
  xe2om=(ffo2*csg2+anpl2*sng2)/deno
!  ye2om=(ffo2*sng2+anpl2*csg2)/deno
!  xe2xm=(ffx2*csg2+anpl2*sng2)/denx
!  ye2xm=(ffx2*sng2+anpl2*csg2)/denx
  coffpol=xe2om
end subroutine pol_limit
!-------------------------------------------------------
subroutine printtraj
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! PRINTTRAJ                                                            C
!   uses                                                               C
!   sets                                                               C
!   calls                                                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  use kinds, only : r8
  use spece_data, only : nmx
  implicit none
  integer iox,iif,i,iengy,itraj,nnzero,iprtr,iprfr
  integer nharmv(nmx),iprn(nmx)
  real(kind=r8) rrr(nmx),ygv(nmx)
  real(kind=r8) eresmn(nmx),eresmx(nmx),fghz
  real(kind=r8) anplv(nmx),ggv(nmx)
  real(kind=r8) stv(nmx),zzv(nmx),phiv(nmx),xgv(nmx),anprv(nmx),psjki(nmx)
  real(kind=r8) tauv(nmx),alphav(nmx),tev(nmx),ddv(nmx),ddrmxv(nmx)
  real(kind=r8) ygn,umnnpl2,gcentr,a

  common/atemt/ggv
  common/ygv/ygv
  common/imode/iox
  common/fghz/fghz
  common/iiv/iif
  common/rrr/rrr
  common/nplv/anplv
  common/harmv/nharmv
  common/traj/stv,zzv,phiv,xgv,anprv
  common/atjki/tauv,alphav,tev
  common/psjki/psjki
  common/iprint/iprtr,iprfr
  common/ddv/ddv,ddrmxv

  iengy=20+iox
  itraj=8+iox
  if (iprtr>0) write(itraj,*) ' '
  nnzero=0
  do i=1,iif
    if (alphav(i)>0.0d0) then
      nnzero=nnzero+1               ! sets iprn to print data
      iprn(nnzero)=i                ! for resonance in this point
      umnnpl2=1.0d0-anplv(i)*anplv(i)
      ygn=ygv(i)*nharmv(i)
      gcentr=ygn/umnnpl2
      a=abs(anplv(i))*sqrt(ygn*ygn-umnnpl2)/umnnpl2
      eresmn(nnzero)=(gcentr-a-1.0d0)*511.053d0  ! min and max
      eresmx(nnzero)=(gcentr+a-1.0d0)*511.053d0  ! resonant energy 
    end if
    if (iprtr>0.and.mod(i,iprtr)==0)                           &
                write(itraj,99) fghz,stv(i),rrr(i)*1.0d-2,zzv(i),   &
                phiv(i),xgv(i),ygv(i),anprv(i),anplv(i),psjki(i),ddv(i),ddrmxv(i)
  end do
  do i=1,nnzero
    if (i==1.or.(i>1.and.(iprn(i)-1)/=iprn(i-1))) then
      write(iengy,*) ' '
    end if
    write(iengy,99) fghz,stv(iprn(i)),rrr(iprn(i)),anplv(iprn(i)),  &
     tauv(iprn(i)),alphav(iprn(i)),tev(iprn(i)),ggv(iprn(i)),       &
     eresmn(i),eresmx(i),dble(nharmv(iprn(i)))
  end do
  return
  99 format(20(1x,e16.8e3))
end subroutine printtraj
!-------------------------------------------------------
subroutine points_ox(rz,zz,rf,zf,psinvf,info)
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
  real(kind=r8), intent(in)  :: rz, zz
  real(kind=r8), intent(out) :: rf, zf, psinvf
  integer      , intent(out) :: info
! local variables
  integer, parameter :: n = 2, ldfjac = n, lwa = (n*(n+13))/2
  real(kind=r8) :: xvec(n),fvec(n),fjac(ldfjac,n),wa(lwa),tol
! common variables
  real(kind=r8) psinv
  common/psival/psinv
! external subroutines/functions
  interface
    subroutine fcnox(n,x,fvec,fjac,ldfjac,iflag)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: n, ldfjac, iflag
      real(kind=r8), intent(in)    :: x(n)
      real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fcnox
    subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: n, ldfjac, lwa
      integer      , intent(out)   :: info
      real(kind=r8), intent(in)    :: tol
      real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), wa(lwa)
      real(kind=r8), intent(inout) :: x(n)
      external fcn
    end subroutine hybrj1
  end interface
! executable statements
  xvec(1)=rz
  xvec(2)=zz
  tol = sqrt(eps_r8)
  call hybrj1(fcnox,n,xvec,fvec,fjac,ldfjac,tol,info,wa,lwa)
  rf=xvec(1)
  zf=xvec(2)
  psinvf=psinv
end subroutine points_ox
!-------------------------------------------------------
subroutine fcnox(n,x,fvec,fjac,ldfjac,iflag)
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in) :: n, ldfjac, iflag
  real(kind=r8), intent(in)    :: x(n)
  real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
! common variables
  real(kind=r8) :: dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv,psia
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/pareq1/psia
! external subroutines
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface
! executable statements
  call equinum(x(1),x(2))
  if (iflag/=2) then
    fvec(1) = dpsidr/psia
    fvec(2) = dpsidz/psia
  else
    fjac(1,1) = ddpsidrr/psia
    fjac(1,2) = ddpsidrz/psia
    fjac(2,1) = ddpsidrz/psia
    fjac(2,2) = ddpsidzz/psia
  end if
end subroutine fcnox
!-------------------------------------------------------
subroutine points_tgo(rz,zz,rf,zf,psin,info)
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
! arguments
  real(kind=r8), intent(in)  :: rz,zz,psin
  real(kind=r8), intent(out) :: rf,zf
  integer      , intent(out) :: info
! local variables
  integer, parameter :: n = 2, ldfjac = n, lwa = (n*(n+13))/2
  real(kind=r8) :: xvec(n),fvec(n),fjac(ldfjac,n),wa(lwa),tol
! common variables
  real(kind=r8) :: h
  common/cnpsi/h
! external subroutines/functions
  interface
    subroutine fcntgo(n,x,fvec,fjac,ldfjac,iflag)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, ldfjac, iflag
      real(kind=r8), intent(in)    :: x(n)
      real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fcntgo
    subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: n, ldfjac, lwa
      integer      , intent(out)   :: info
      real(kind=r8), intent(in)    :: tol
      real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), wa(lwa)
      real(kind=r8), intent(inout) :: x(n)
      external fcn
    end subroutine hybrj1
  end interface
! executable statements
  h=psin
  xvec(1)=rz
  xvec(2)=zz
  tol = sqrt(eps_r8)
  call hybrj1(fcntgo,n,xvec,fvec,fjac,ldfjac,tol,info,wa,lwa)
  rf=xvec(1)
  zf=xvec(2)
end subroutine points_tgo
!-------------------------------------------------------
subroutine fcntgo(n,x,fvec,fjac,ldfjac,iflag)
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in) :: n, ldfjac, iflag
  real(kind=r8), intent(in)    :: x(n)
  real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
! common variables
  real(kind=r8) :: dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  real(kind=r8) :: psinv,h,psia
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/psival/psinv
  common/cnpsi/h
  common/pareq1/psia
! external subroutines
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface
! executable statements
  call equinum(x(1),x(2))
  if (iflag/=2) then
    fvec(1) = psinv-h
    fvec(2) = dpsidr/psia
  else
    fjac(1,1) = dpsidr/psia
    fjac(1,2) = dpsidz/psia
    fjac(2,1) = ddpsidrr/psia
    fjac(2,2) = ddpsidrz/psia
  end if
end subroutine fcntgo
!-------------------------------------------------------
subroutine points_mino(rz,zz,rf,psinf,info)
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
! arguments
  real(kind=r8), intent(in)  :: rz, zz
  real(kind=r8), intent(out) :: rf, psinf
  integer      , intent(out) :: info
! local variables
  integer, parameter :: n = 2, ldfjac = n, lwa = (n*(n+13))/2
  real(kind=r8) :: xvec(n), fvec(n), fjac(ldfjac,n), wa(lwa), tol
! common variables
  real(kind=r8) h, psinv
  common/cnz/h
  common/psival/psinv
! external subroutines
  interface
    subroutine fcnstazr(n,x,fvec,fjac,ldfjac,iflag)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, ldfjac, iflag
      real(kind=r8), intent(in)    :: x(n)
      real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fcnstazr
    subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: n, ldfjac, lwa
      integer      , intent(out)   :: info
      real(kind=r8), intent(in)    :: tol
      real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), wa(lwa)
      real(kind=r8), intent(inout) :: x(n)
      external fcn
    end subroutine hybrj1
  end interface
! executable statements
  h=zz
  xvec(1)=rz
  xvec(2)=zz
  tol = sqrt(eps_r8)
  call hybrj1(fcnstazr,n,xvec,fvec,fjac,ldfjac,tol,info,wa,lwa)
  rf=xvec(1)
  psinf=psinv
end subroutine points_mino
!-------------------------------------------------------
subroutine fcnstazr(n,x,fvec,fjac,ldfjac,iflag)
  use kinds, only : r8
  implicit none
! arguments
  integer      , intent(in)    :: n, ldfjac, iflag
  real(kind=r8), intent(in)    :: x(n)
  real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
! common variables
  real(kind=r8) :: dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv,h,psia
  common/derip/dpsidr,dpsidz,ddpsidrr,ddpsidzz,ddpsidrz,fpolv,ffpv
  common/cnz/h
  common/pareq1/psia
! external subroutines
  interface
    subroutine equinum(rpsim,zpsim)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in) :: rpsim, zpsim
    end subroutine equinum
  end interface
! executable statements
  call equinum(x(1),x(2))
  if (iflag/=2) then
    fvec(1) = x(2)-h
    fvec(2) = dpsidr/psia
  else
    fjac(1,1) = 0.0d0
    fjac(1,2) = 1.0d0
    fjac(2,1) = ddpsidrr/psia
    fjac(2,2) = ddpsidrz/psia
  end if
end subroutine fcnstazr
!-------------------------------------------------------
subroutine rhotor
  use kinds, only : r8
  use spece_data, only : nnw
  implicit none
  integer :: iopt,ier,k
  real(r8) :: dx
! common variables
  integer :: nr,nz,nppn,nppt
  real(kind=r8), dimension(nnw) :: psinr,rhotnr,phit
  real(kind=r8), dimension(nnw,4) :: cq,crhot
  common/eqnn/nr,nz,nppn,nppt
  common/psinr/psinr
  common/coffq/cq
  common/coffrhot/crhot
! external subroutines
  interface
   subroutine difcs(xx,yy,iopt,cc,ier)
     use kinds, only : r8
     implicit none
     integer      , intent(in)    :: iopt
     real(kind=r8), intent(in)    :: xx(:),yy(:)
     real(kind=r8), intent(inout) :: cc(:,:)
     integer      , intent(out)   :: ier
   end subroutine difcs
  end interface
! executable statements
  phit(1)=0.0_r8
  do k=1,nr-1
    dx=psinr(k+1)-psinr(k)
    phit(k+1)=phit(k)+dx*(cq(k,1)+dx*(cq(k,2)/2.0_r8+dx*(cq(k,3)/3.0_r8+dx*cq(k,4)/4.0_r8)))
  end do
  rhotnr(1:nr)=sqrt(phit(1:nr)/phit(nr))
! rhotor spline interpolation
  iopt=0
  call difcs(psinr(:nr),rhotnr(:nr),iopt,crhot,ier)
end subroutine
