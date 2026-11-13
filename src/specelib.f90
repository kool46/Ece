integer function length(string)
! returns length of string ignoring trailing blanks
  implicit none
! arguments
  character(len=*), intent(in) :: string
! local variables
  integer :: i, l
! intrinsic functions
!  integer, intrinsic :: len
! executable statements
  l = len(string)
  do i = l, 1, -1
    if(string(i:i) /= ' ') exit
  end do
  length = i
end function length
!-------------------------------------------------------
integer function leadblank(string) 
! returns number of leading blanks 
  implicit none
! arguments
  character(len=*), intent(in) :: string
! local variables
  integer :: i,lstr
! external functions
  interface
    integer function length(string)
      character(len=*), intent(in) :: string
    end function length
  end interface
! executable statements
  i = 1
  lstr = length(string)
  do while (string(i:i) == ' ' .and. i < lstr)
    i = i+1
  end do
  leadblank = i-1
end function leadblank
!-------------------------------------------------------
!-------------------------------------------------------
function locate(xx,x)
! Given an array xx(n), and a value x 
! returns a value j such that xx(j) < x < xx(j+1)
! xx(n) must be monotonic, either increasing or decreasing. 
! j=0 or j=n  indicate that x is out of range  (Numerical Recipes)
  use kinds, only : r8
  implicit none
  integer locate
! arguments
  real(kind=r8), intent(in)  :: xx(:), x
! local variables
  integer :: n, jl, ju, jm
  logical :: incrx
! executable statements
  n = size(xx)
  jl = 0
  ju = n+1
  incrx = (xx(n) > xx(1))
! incr == .true. for increasing xx, .false. otherwise
  do while ((ju-jl) > 1)
    jm = (ju+jl)/2
    if (incrx .eqv. (x > xx(jm))) then
      jl = jm
    else
      ju = jm
    endif
  end do
  locate = jl
end function locate
!-------------------------------------------------------
subroutine intlin(x1,y1,x2,y2,x,y)
! linear interpolation
! (x1,y1) < (x,y) < (x2,y2)
  use kinds, only : r8
  use constants, only : tiny_r8
  implicit none
! arguments
  real(kind=r8), intent(in)  :: x1, y1, x2, y2, x
  real(kind=r8), intent(out) :: y
! local variables
  real(kind=r8) :: a, b
! executable statements
  if (x2 /= x1) then
    a = (x2-x)/(x2-x1)
    b = 1._r8 - a
    y = a*y1 + b*y2
  else
    if (x==x1 .and. y2==y1) then
      y = y1
    else
      print*,'Error in INTLIN'
    end if
  endif
end subroutine intlin
!-------------------------------------------------------
subroutine swap(a,b)
  use kinds, only : r8
  implicit none
  real(kind=r8), intent(inout) :: a, b
  real(kind=r8) :: temp
  temp=a
  a=b
  b=temp
end subroutine swap
!-------------------------------------------------------
!-------------------------------------------------------
function fact(k)
! factorial function
  use kinds, only : r8
  implicit none
  real(kind=r8) fact
! arguments
  integer, intent(in) :: k
! local variables
  integer :: i
! executable statements
  if(k<0) then
    fact=0.0d0
  else
    fact=1.0d0
    do i=2,k
      fact=fact*i
    end do
  end if
end function fact
!-------------------------------------------------------
function gammln(xx)
! Returns the value ln[Gamma(xx)] for xx > 0.
  use kinds, only : r8
  implicit none
  real(kind=r8) gammln
! arguments
  real(kind=r8), intent(in) :: xx
! parameters
  real(kind=r8), parameter :: stp = 2.5066282746310005d0
  real(kind=r8), dimension(6), parameter :: cof=(/76.18009172947146d0,  &
        -86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,  &
        .1208650973866179d-2,-.5395239384953d-5/)
! local variables
  integer       :: j
  real(kind=r8) :: ser, tmp, x, y
! executable statements
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  end do
  gammln=tmp+log(stp*ser/x)
end function gammln
!-------------------------------------------------------
subroutine zetac(xi, yi, zr, zi, iflag)
! PLASMA DISPERSION FUNCTION Z of complex argument
! Z(z) = i sqrt(pi) w(z)
! Function w(z) from:
! algorithm 680, collected algorithms from acm.
! this work published in transactions on mathematical software,
! vol. 16, no. 1, pp. 47.
!
!  given a complex number z = (xi,yi), this subroutine computes
!  the value of the faddeeva-function w(z) = exp(-z**2)*erfc(-i*z),
!  where erfc is the complex complementary error-function and i
!  means sqrt(-1).
!  the accuracy of the algorithm for z in the 1st and 2nd quadrant
!  is 14 significant digits; in the 3rd and 4th it is 13 significant
!  digits outside a circular region with radius 0.126 around a zero
!  of the function.
!  all real variables in the program are double precision.
!
!
!  the code contains a few compiler-dependent parameters :
! rmaxreal = the maximum value of rmaxreal equals the root of
!            rmax = the largest number which can still be
!            implemented on the computer in double precision
!            floating-point arithmetic
! rmaxexp  = ln(rmax) - ln(2)
! rmaxgoni = the largest possible argument of a double precision
!            goniometric function (dcos, dsin, ...)
!  the reason why these parameters are needed as they are defined will
!  be explained in the code by means of comments
!
!
!  parameter list
! xi     = real      part of z
! yi     = imaginary part of z
! u      = real      part of w(z)
! v      = imaginary part of w(z)
! iflag   = an error flag indicating whether overflow will
!          occur or not; type integer;
!          the values of this variable have the following
!          meaning :
!          iflag=0 : no error condition
!          iflag=1 : overflow will occur, the routine
!                         becomes inactive
!  xi, yi      are the input-parameters
!  u, v, iflag  are the output-parameters
!
!  furthermore the parameter factor equals 2/sqrt(pi)
!
!  the routine is not underflow-protected but any variable can be
!  put to 0 upon underflow;
!
!  reference - gpm poppe, cmj wijers; more efficient computation of
!  the complex error-function, acm trans. math. software.
!
  use kinds, only : r8
  implicit none
! arguments
  real(kind=r8), intent(in)  :: xi, yi
  real(kind=r8), intent(out) :: zr, zi
  integer      , intent(out) :: iflag
! local variables
  integer       :: n, j, i, kapn, nu, np1
  real(kind=r8) :: u, v, x, y, xabs, yabs, qrho, xabsq, xquad, yquad,   &
                   xsum, ysum, xaux, u1, v1, daux, u2, v2, h, h2,       &
                   qlambda, rx, ry, sx, sy, tx, ty, c, w1
! parameters
  real(kind=r8), parameter ::      factor   = 1.12837916709551257388d0, &
     rpi     = 2.0d0/factor,       rmaxreal = 0.5d+154,                 &
     rmaxexp = 708.503061461606d0, rmaxgoni = 3.53711887601422d+15
!
  iflag = 0
  xabs = abs(xi)
  yabs = abs(yi)
  x    = xabs/6.3d0
  y    = yabs/4.4d0
!
! the following if-statement protects
! qrho = (x**2 + y**2) against overflow
!
  if ((xabs>rmaxreal).or.(yabs>rmaxreal)) then
    iflag = 1
    return
  end if
  qrho = x**2 + y**2
  xabsq = xabs**2
  xquad = xabsq - yabs**2
  yquad = 2*xabs*yabs
  if (qrho<0.085264d0) then
!
!  if (qrho<0.085264d0) then the faddeeva-function is evaluated
!  using a power-series (abramowitz/stegun, equation (7.1.5), p.297)
!  n is the minimum number of terms needed to obtain the required
!  accuracy
!
    qrho  = (1-0.85d0*y)*sqrt(qrho)
    n     = idnint(6 + 72*qrho)
    j     = 2*n+1
    xsum  = 1.0d0/j
    ysum  = 0.0d0
    do i=n, 1, -1
      j    = j - 2
      xaux = (xsum*xquad - ysum*yquad)/i
      ysum = (xsum*yquad + ysum*xquad)/i
      xsum = xaux + 1.0d0/j
    end do
    u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0d0
    v1   =  factor*(xsum*xabs - ysum*yabs)
    daux =  dexp(-xquad)
    u2   =  daux*dcos(yquad)
    v2   = -daux*dsin(yquad)
    u    = u1*u2 - v1*v2
    v    = u1*v2 + v1*u2
  else
!
!  if (qrho>1.o) then w(z) is evaluated using the laplace
!  continued fraction
!  nu is the minimum number of terms needed to obtain the required
!  accuracy
!
!  if ((qrho>0.085264d0).and.(qrho<1.0)) then w(z) is evaluated
!  by a truncated taylor expansion, where the laplace continued fraction
!  is used to calculate the derivatives of w(z)
!  kapn is the minimum number of terms in the taylor expansion needed
!  to obtain the required accuracy
!  nu is the minimum number of terms of the continued fraction needed
!  to calculate the derivatives with the required accuracy
!
    if (qrho>1.0d0) then
      h    = 0.0d0
      kapn = 0
      qrho = sqrt(qrho)
      nu   = idint(3 + (1442/(26*qrho+77)))
    else
      qrho = (1-y)*sqrt(1-qrho)
      h    = 1.88d0*qrho
      h2   = 2*h
      kapn = idnint(7  + 34*qrho)
      nu   = idnint(16 + 26*qrho)
    endif
    if (h>0.0d0) qlambda = h2**kapn
    rx = 0.0d0
    ry = 0.0d0
    sx = 0.0d0
    sy = 0.0d0
    do n=nu, 0, -1
      np1 = n + 1
      tx  = yabs + h + np1*rx
      ty  = xabs - np1*ry
      c   = 0.5d0/(tx**2 + ty**2)
      rx  = c*tx
      ry  = c*ty
      if ((h>0.0d0).and.(n<=kapn)) then
        tx = qlambda + sx
        sx = rx*tx - ry*sy
        sy = ry*tx + rx*sy
        qlambda = qlambda/h2
      endif
    end do
    if (h==0.0d0) then
      u = factor*rx
      v = factor*ry
    else
      u = factor*sx
      v = factor*sy
    end if
    if (yabs==0.0d0) u = dexp(-xabs**2)
  end if
!
!  evaluation of w(z) in the other quadrants
!
  if (yi<0.0d0) then
    if (qrho<0.085264d0) then
      u2    = 2*u2
      v2    = 2*v2
    else
      xquad =  -xquad
!
!     the following if-statement protects 2*exp(-z**2)
!     against overflow
!
      if ((yquad>rmaxgoni).or.(xquad>rmaxexp)) then
        iflag=1
        return
      end if
      w1 =  2.0d0*dexp(xquad)
      u2  =  w1*dcos(yquad)
      v2  = -w1*dsin(yquad)
    end if
    u = u2 - u
    v = v2 - v
    if (xi>0.0d0) v = -v
  else
    if (xi<0.0d0) v = -v
  end if
  zr = -v*rpi
  zi =  u*rpi
end subroutine zetac
!-------------------------------------------------------
function expei(x)
! derived from calcei for use with  int=3 only
! result ->  y = expei(x)
!                               
! computes the exponential integral exp(-x)*ei(x)  for real arguments x
! where
!
!       integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
!  ei(x) =
!      -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
!
!  and where the first integral is a principal value integral.
!
  use kinds, only : r8
  implicit none
  real(kind=r8) :: expei
! arguments
  real(kind=r8), intent(in)  :: x
! local variables
  integer       :: i, ip1, im1
  real(kind=r8) :: frac, sump, sumq, t, w, xmx0, y, ysq
  real(kind=r8), dimension(10) :: px, qx
  real(kind=r8), parameter :: zero=0.0d0, p037=0.037d0, half=0.5d0, &
       one=1.0d0, two=2.0d0, three=3.0d0, four=4.0d0, six=6.0d0,    &
       twelve=12.0d0, two4=24.0d0, x01=381.5d0, x11=1024.0d0,       &
       x02=-5.1182968633365538008d-5, x0=3.7250741078136663466d-1
! machine-dependent constants
  real(kind=r8), parameter :: xinf=1.79d+308
! coefficients  for -1.0 <= x < 0.0
  real(kind=r8), dimension(7), parameter ::                         &
     a=(/1.1669552669734461083368d2, 2.1500672908092918123209d3,    &
         1.5924175980637303639884d4, 8.9904972007457256553251d4,    &
         1.5026059476436982420737d5,-1.4815102102575750838086d5,    &
         5.0196785185439843791020d0/)
  real(kind=r8), dimension(6), parameter ::                         &
     b=(/4.0205465640027706061433d1, 7.5043163907103936624165d2,    &
         8.1258035174768735759855d3, 5.2440529172056355429883d4,    &
         1.8434070063353677359298d5, 2.5666493484897117319268d5/)
! coefficients for -4.0 <= x < -1.0
  real(kind=r8), dimension(9), parameter ::                         &
     c=(/3.828573121022477169108d-1, 1.107326627786831743809d+1,    &
         7.246689782858597021199d+1, 1.700632978311516129328d+2,    &
         1.698106763764238382705d+2, 7.633628843705946890896d+1,    &
         1.487967702840464066613d+1, 9.999989642347613068437d-1,    &
         1.737331760720576030932d-8/),                              &
     d=(/8.258160008564488034698d-2, 4.344836335509282083360d+0,    &
         4.662179610356861756812d+1, 1.775728186717289799677d+2,    &
         2.953136335677908517423d+2, 2.342573504717625153053d+2,    &
         9.021658450529372642314d+1, 1.587964570758947927903d+1,    &
         1.000000000000000000000d+0/)
! coefficients for x < -4.0
  real(kind=r8), dimension(10), parameter ::                        &
     e=(/1.3276881505637444622987d+2,3.5846198743996904308695d+4,   &
         1.7283375773777593926828d+5,2.6181454937205639647381d+5,   &
         1.7503273087497081314708d+5,5.9346841538837119172356d+4,   &
         1.0816852399095915622498d+4,1.0611777263550331766871d03,   &
         5.2199632588522572481039d+1,9.9999999999999999087819d-1/), &
     f=(/3.9147856245556345627078d+4,2.5989762083608489777411d+5,   &
         5.5903756210022864003380d+5,5.4616842050691155735758d+5,   &
         2.7858134710520842139357d+5,7.9231787945279043698718d+4,   &
         1.2842808586627297365998d+4,1.1635769915320848035459d+3,   &
         5.4199632588522559414924d+1,1.0d0/)
!  coefficients for rational approximation to ln(x/a), |1-x/a| < .1
  real(kind=r8), dimension(4), parameter ::                         &
     plg=(/-2.4562334077563243311d+01,2.3642701335621505212d+02,    &
           -5.4989956895857911039d+02,3.5687548468071500413d+02/),  &
     qlg=(/-3.5553900764052419184d+01,1.9400230218539473193d+02,    &
           -3.3442903192607538956d+02,1.7843774234035750207d+02/)
! coefficients for  0.0 < x < 6.0,
!  ratio of chebyshev polynomials
  real(kind=r8), dimension(10), parameter ::                         &
     p=(/-1.2963702602474830028590d01,-1.2831220659262000678155d03,  &
         -1.4287072500197005777376d04,-1.4299841572091610380064d06,  &
         -3.1398660864247265862050d05,-3.5377809694431133484800d08,  &
          3.1984354235237738511048d08,-2.5301823984599019348858d10,  &
          1.2177698136199594677580d10,-2.0829040666802497120940d11/),&
     q=(/ 7.6886718750000000000000d01,-5.5648470543369082846819d03,  &
          1.9418469440759880361415d05,-4.2648434812177161405483d06,  &
          6.4698830956576428587653d07,-7.0108568774215954065376d08,  &
          5.4229617984472955011862d09,-2.8986272696554495342658d10,  &
          9.8900934262481749439886d10,-8.9673749185755048616855d10/)
! j-fraction coefficients for 6.0 <= x < 12.0
  real(kind=r8), dimension(10), parameter ::                        &
     r=(/-2.645677793077147237806d00,-2.378372882815725244124d00,   &
         -2.421106956980653511550d01, 1.052976392459015155422d01,   &
          1.945603779539281810439d01,-3.015761863840593359165d01,   &
          1.120011024227297451523d01,-3.988850730390541057912d00,   &
          9.565134591978630774217d00, 9.981193787537396413219d-1/)
  real(kind=r8), dimension(9), parameter ::                         &
     s=(/ 1.598517957704779356479d-4, 4.644185932583286942650d00,   &
          3.697412299772985940785d02,-8.791401054875438925029d00,   &
          7.608194509086645763123d02, 2.852397548119248700147d01,   &
          4.731097187816050252967d02,-2.369210235636181001661d02,   &
          1.249884822712447891440d00/)
! j-fraction coefficients for 12.0 <= x < 24.0
  real(kind=r8), dimension(10), parameter ::                        &
     p1=(/-1.647721172463463140042d00,-1.860092121726437582253d01,  &
          -1.000641913989284829961d01,-2.105740799548040450394d01,  &
          -9.134835699998742552432d-1,-3.323612579343962284333d01,  &
           2.495487730402059440626d01, 2.652575818452799819855d01,  &
          -1.845086232391278674524d00, 9.999933106160568739091d-1/)
  real(kind=r8), dimension(9), parameter ::                         &
     q1=(/ 9.792403599217290296840d01, 6.403800405352415551324d01,  &
           5.994932325667407355255d01, 2.538819315630708031713d02,  &
           4.429413178337928401161d01, 1.192832423968601006985d03,  &
           1.991004470817742470726d02,-1.093556195391091143924d01,  &
           1.001533852045342697818d00/)
! j-fraction coefficients for  x >= 24.0
  real(kind=r8), dimension(10), parameter ::                        &
     p2=(/ 1.75338801265465972390d02,-2.23127670777632409550d02,    &
          -1.81949664929868906455d01,-2.79798528624305389340d01,    &
          -7.63147701620253630855d00,-1.52856623636929636839d01,    &
          -7.06810977895029358836d00,-5.00006640413131002475d00,    &
          -3.00000000320981265753d00, 1.00000000000000485503d00/)
  real(kind=r8), dimension(9), parameter ::                         &
     q2=(/ 3.97845977167414720840d04, 3.97277109100414518365d00,    &
           1.37790390235747998793d02, 1.17179220502086455287d02,    &
           7.04831847180424675988d01,-1.20187763547154743238d01,    &
          -7.99243595776339741065d00,-2.99999894040324959612d00,    &
           1.99999999999048104167d00/)
!
  if (x == zero) then
    expei = -xinf
  else if (x < zero) then
!
! calculate ei for negative argument or for e1.
!
    y = abs(x)
    if (y <= one) then
      sump = a(7) * y + a(1)
      sumq = y + b(1)
      do i = 2, 6
        sump = sump * y + a(i)
        sumq = sumq * y + b(i)
      end do
      expei = (log(y) - sump / sumq ) * exp(y)
    else if (y <= four) then
      w = one / y
      sump = c(1)
      sumq = d(1)
      do i = 2, 9
        sump = sump * w + c(i)
        sumq = sumq * w + d(i)
      end do
      expei = - sump / sumq
    else
      w = one / y
      sump = e(1)
      sumq = f(1)
      do i = 2, 10
        sump = sump * w + e(i)
        sumq = sumq * w + f(i)
      end do
      expei = w * (w * sump / sumq - one)
    end if
  else if (x < six) then
!
!  to improve conditioning, rational approximations are expressed
!in terms of chebyshev polynomials for 0 <= x < 6, and in
!continued fraction form for larger x.
!
    t = (x + x) / three - two
    px(1) = zero
    qx(1) = zero
    px(2) = p(1)
    qx(2) = q(1)
    do i = 2, 9
      ip1 = i+1
      im1 = i-1
      px(ip1) = t * px(i) - px(im1) + p(i)
      qx(ip1) = t * qx(i) - qx(im1) + q(i)
    end do
    sump = half * t * px(10) - px(9) + p(10)
    sumq = half * t * qx(10) - qx(9) + q(10)
    frac = sump / sumq
    xmx0 = (x - x01/x11) - x02
    if (abs(xmx0) >= p037) then
      expei = exp(-x) * ( log(x/x0) + xmx0 * frac )
    else
!
! special approximation to  ln(x/x0)  for x close to x0
!
      y = xmx0 / (x + x0)
      ysq = y*y
      sump = plg(1)
      sumq = ysq + qlg(1)
      do i = 2, 4
        sump = sump*ysq + plg(i)
        sumq = sumq*ysq + qlg(i)
      end do
      expei = exp(-x) * (sump / (sumq*(x+x0)) + frac) * xmx0
    end if
  else if (x < twelve) then
    frac = zero
    do i = 1, 9
      frac = s(i) / (r(i) + x + frac)
    end do
    expei = (r(10) + frac) / x
  else if (x <= two4) then
    frac = zero
    do i = 1, 9
      frac = q1(i) / (p1(i) + x + frac)
    end do
    expei = (p1(10) + frac) / x
  else
    y = one / x
    frac = zero
    do i = 1, 9
      frac = q2(i) / (p2(i) + x + frac)
    end do
    frac = p2(10) + frac
    expei = y + y * y * frac
  end if
end function expei
!-------------------------------------------------------
!-------------------------------------------------------
function spli(cspli,k,dx)
! spline evaluation at x distance dx from k-th knot
! coefficient matrix must have dimension (n,4), n >= k
  use kinds, only : r8
  implicit none
  real(kind=r8) spli
! arguments
  integer      , intent(in) :: k
  real(kind=r8), intent(in) :: cspli(:,:), dx
! executable statements
  spli = cspli(k,1) + dx*(cspli(k,2) + dx*(cspli(k,3) + dx*cspli(k,4)))
end function spli
!-------------------------------------------------------
subroutine difcs(x,y,iopt,c,ier)
  use kinds, only : r8
  implicit none
! arguments
  integer      , intent(in)    :: iopt
  real(kind=r8), intent(in)    :: x(:),y(:)
  real(kind=r8), intent(inout) :: c(:,:)
  integer      , intent(out)   :: ier
! local variables
  integer       :: n, iol, ioh, j, i, j1, ii
  real(kind=r8) :: xb, xc, ya, yb, h, a, r, dyb, dy2, dya
! executable statements
  n = size(x)
  if (n>1) then
!
!   initialization
!
    xc = x(1)
    yb = y(1)
    h = 0.0d0
    a = 0.0d0
    r = 0.0d0
    dyb = 0.0d0
!
!   iol=0 - given derivative at first point
!   ioh=0 - given derivative at last point
!
    if (iopt==3) then
      iol = 0
      ioh = 0
    else
      iol = iopt-1
      ioh = iopt-2
    end if
    dy2 = c(2,1)
!
!   form the system of linear equations
!   and eliminate subsequentially
!
    j = 1
    ier = 0
    do i = 1,n
      a = h*(2.0d0-a)
      dya = dyb+h*r
      if (i<n) then
        j = j+1
        xb = xc
        xc = x(j)
        h = xc-xb
!
!       ii=0 - increasing abscissae
!       ii=1 - decreasing abscissae
!
        if (h==0.0d0) then
          ier = 1
          exit
        else if (h<0.0d0) then
          ii = 1
        else
          ii = 0
        end if
        ya = yb
        yb = y(j)
        dyb = (yb-ya)/h
        if (i==1) then
          j1 = ii
          if (iol==0) dya = c(1,1)
        end if
      else
!
!       set derivative dy2 at last point
!
        h = 0.0d0
        if (ioh==0) then
          dyb = dy2
        else
          dyb = dya
        end if
      end if
      if ((i>1.or.iol==0) .and. &
          (i<n.or.ioh==0)) then
        if (j1/=ii) then
          ier = 1
          exit
        end if
        a = 1.0d0/(h+h+a)
      end if
      r = a*(dyb-dya)
      c(i,3) = r
      a = h*a
      c(i,2) = a
      c(i,1) = dyb
    end do
    if (ier==0) then
!
!     back substitution of the system of linear equations
!     and computation of the other coefficients
!
      a = 1.0d0
      do i = n,1,-1
        xb = x(i)
        h = xc-xb
        xc = xb
        a = a+h
        yb = r
        r = c(i,3)-r*c(i,2)
        ya = r+r
        c(i,3) = ya+r
        c(i,2) = c(i,1)-h*(ya+yb)
        if (i<n) c(i+ii,4) = (yb-r)/a
        if (i==n) then
          if (ii==0) then
            c(n,4) = (yb-r)/a
          else
            c(1,4) = (yb-r)/a
          end if
        end if
        c(i,1) = y(i)
        a = 0.0d0
      end do
    end if
  end if
end subroutine difcs
!-------------------------------------------------------
!-------------------------------------------------------
! spline routines from netlib/DIERCKX
!  
subroutine fpback(a,z,n,k,c,nest)
!  subroutine fpback calculates the solution of the system of
!  equations a*c = z with a a n x n upper triangular matrix
!  of bandwidth k.
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)          :: n, k, nest
  real(kind=r8), intent(in)    :: a(nest,k), z(n)
  real(kind=r8), intent(inout) :: c(n)
! local variables
  real(kind=r8) :: store
  integer       :: i, i1, j, k1, l, m
! executable statements
  k1 = k-1
  c(n) = z(n)/a(n,1)
  i = n-1
  if(i==0) return
  do j=2,n
    store = z(i)
    i1 = k1
    if(j<=k1) i1 = j-1
    m = i
    do l=1,i1
      m = m+1
      store = store-c(m)*a(i,l+1)
    end do
    c(i) = store/a(i,1)
    i = i-1
  end do
end subroutine fpback
!-------------------------------------------------------
subroutine fpbspl(t,n,k,x,l,h)
!  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
!  degree k at t(l) <= x < t(l+1) using the stable recurrence
!  relation of de boor and cox.
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in) :: n, k, l
  real(kind=r8), intent(in) :: x, t(n)
  real(kind=r8), intent(out) :: h(6)
! local variables
  integer :: i, j, li, lj
  real(kind=r8) :: f, hh(5)
! parameters
  real(kind=r8), parameter :: one = 1.0d0
! executable statements
  h(1) = one
  do j=1,k
    do i=1,j
      hh(i) = h(i)
    end do
    h(1) = 0.0d0
    do i=1,j
      li = l+i
      lj = li-j
      f = hh(i)/(t(li)-t(lj))
      h(i) = h(i)+f*(t(li)-x)
      h(i+1) = f*(x-t(lj))
    end do
  end do
end subroutine fpbspl
!-------------------------------------------------------
subroutine fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wx,wy,lx,ly)
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)  :: nx, ny, kx, ky, mx, my
  integer, intent(out) :: lx(mx), ly(my)
  real(kind=r8), intent(in)  :: tx(nx), ty(ny), c((nx-kx-1)*(ny-ky-1)), &
      x(mx), y(my)
  real(kind=r8), intent(out) :: z(mx*my)
  real(kind=r8), intent(out) :: wx(mx,kx+1), wy(my,ky+1)
! local variables
  integer       :: kx1, ky1, l, l1, l2, m, nkx1, nky1, i, j, i1, j1
  real(kind=r8) :: arg, sp, tb, te, h(6)
! external subroutines
  interface
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
  end interface
! executable statements
  kx1 = kx+1
  nkx1 = nx-kx1
  tb = tx(kx1)
  te = tx(nkx1+1)
  l = kx1
  l1 = l+1
  do i=1,mx
    arg = x(i)
    if(arg<tb) arg = tb
    if(arg>te) arg = te
    do while (arg>=tx(l1) .and. l/=nkx1)
      l = l1
      l1 = l+1
    end do
    call fpbspl(tx,nx,kx,arg,l,h)
    lx(i) = l-kx1
    do j=1,kx1
      wx(i,j) = h(j)
    end do
  end do
  ky1 = ky+1
  nky1 = ny-ky1
  tb = ty(ky1)
  te = ty(nky1+1)
  l = ky1
  l1 = l+1
  do i=1,my
    arg = y(i)
    if(arg<tb) arg = tb
    if(arg>te) arg = te
    do while (arg>=ty(l1) .and. l/=nky1)
      l = l1
      l1 = l+1
    end do
    call fpbspl(ty,ny,ky,arg,l,h)
    ly(i) = l-ky1
    do j=1,ky1
      wy(i,j) = h(j)
    end do
  end do
  m = 0
  do i=1,mx
    l = lx(i)*nky1
    do i1=1,kx1
      h(i1) = wx(i,i1)
    end do
    do j=1,my
      l1 = l+ly(j)
      sp = 0.0d0
      do i1=1,kx1
        l2 = l1
        do j1=1,ky1
          l2 = l2+1
          sp = sp+c(l2)*h(i1)*wy(j,j1)
        end do
        l1 = l1+nky1
      end do
      m = m+1
      z(m) = sp
    end do
  end do
end subroutine fpbisp
!-------------------------------------------------------
subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,iwrk,kwrk,ier)
!  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
!  b-spline representation.
!
!  calling sequence:
! call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
!* iwrk,kwrk,ier)
!
!  input parameters:
!   tx    : real array, length nx, which contains the position of the
!       knots in the x-direction.
!   nx    : integer, giving the total number of knots in the x-direction
!   ty    : real array, length ny, which contains the position of the
!       knots in the y-direction.
!   ny    : integer, giving the total number of knots in the y-direction
!   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!       b-spline coefficients.
!   kx,ky : integer values, giving the degrees of the spline.
!   x     : real array of dimension (mx).
!       before entry x(i) must be set to the x co-ordinate of the
!       i-th grid point along the x-axis.
!       tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
!   mx    : on entry mx must specify the number of grid points along
!       the x-axis. mx >=1.
!   y     : real array of dimension (my).
!       before entry y(j) must be set to the y co-ordinate of the
!       j-th grid point along the y-axis.
!       ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
!   my    : on entry my must specify the number of grid points along
!       the y-axis. my >=1.
!   wrk   : real array of dimension lwrk. used as workspace.
!   lwrk  : integer, specifying the dimension of wrk.
!       lwrk >= mx*(kx+1)+my*(ky+1)
!   iwrk  : integer array of dimension kwrk. used as workspace.
!   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!
!  output parameters:
!   z     : real array of dimension (mx*my).
!       on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
!       at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
!   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!
!  other subroutines required:
!fpbisp,fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!            6 (1972) 50-62.
!cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!            applics 10 (1972) 134-149.
!dierckx p. : curve and surface fitting with splines, monographs on
!             numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: nx, ny, kx, ky, mx, my, lwrk, kwrk
  integer, intent(out)   :: ier
  integer, intent(inout) :: iwrk(kwrk)
  real(kind=r8), intent(in)    :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
  real(kind=r8), intent(out)   :: z(mx*my)
  real(kind=r8), intent(inout) :: x(mx), y(my), wrk(lwrk)
! local variables
  integer :: i, iw, lwest
! external subroutines
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
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
  ier = 10
  lwest = (kx+1)*mx+(ky+1)*my
  if (lwrk<lwest) return
  if (kwrk<(mx+my)) return
  if (mx<1) return
  do i=2,mx
    if(x(i)<x(i-1)) return
  end do
  if (my<1) return
  do i=2,my
    if(y(i)<y(i-1)) return
  end do
  ier = 0
  iw = mx*(kx+1)+1
  call fpbisp(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk(:iw-1),wrk(iw:),&
   iwrk(:mx),iwrk(mx+1:))
end subroutine bispev
!-------------------------------------------------------
subroutine fpdisc(t,n,k2,b,nest)
!  subroutine fpdisc calculates the discontinuity jumps of the kth
!  derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1)
  use kinds, only : r8
  implicit none
! arguments
  integer      , intent(in)  :: n, k2, nest
  real(kind=r8), intent(in)  :: t(n)
  real(kind=r8), intent(out) :: b(nest,k2)
! local variables
  real(kind=r8) :: an, fac, prod
  integer       :: i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint
  real(kind=r8), dimension(12) :: h
! executable statements
  k1 = k2-1
  k = k1-1
  nk1 = n-k1
  nrint = nk1-k
  an = nrint
  fac = an/(t(nk1+1)-t(k1))
  do l=k2,nk1
    lmk = l-k1
    do j=1,k1
      ik = j+k1
      lj = l+j
      lk = lj-k2
      h(j) = t(l)-t(lk)
      h(ik) = t(l)-t(lj)
    end do
    lp = lmk
    do j=1,k2
      jk = j
      prod = h(j)
      do i=1,k
        jk = jk+1
        prod = prod*h(jk)*fac
      end do
      lk = lp+k1
      b(lmk,j) = (t(lk)-t(lp))/prod
      lp = lp+1
    end do
  end do
end subroutine fpdisc
!-------------------------------------------------------
subroutine fpgivs(piv,ww,cs,sn)
!  subroutine fpgivs calculates the parameters of a givens
!  transformation .
  use kinds, only : r8
  implicit none
! arguments
  real(kind=r8), intent(in)    :: piv
  real(kind=r8), intent(out)   :: cs, sn
  real(kind=r8), intent(inout) :: ww
! local variables
  real(kind=r8) :: dd, store
! parameters
  real(kind=r8), parameter :: one = 1._r8
! intrinsic functions
!  real, intrinsic :: sqrt
! executable statements
  store = abs(piv)
  if(store>=ww) then
    dd = store*sqrt(one+(ww/piv)**2)
  else
    dd = ww*sqrt(one+(piv/ww)**2)
  end if
  cs = ww/dd
  sn = piv/dd
  ww = dd
end subroutine fpgivs
!-------------------------------------------------------
subroutine fprota(cs,sn,a,b)
!  subroutine fprota applies a givens rotation to a and b.
  use kinds, only : r8
  implicit none
! arguments
  real(kind=r8), intent(in)    :: cs, sn
  real(kind=r8), intent(inout) :: a, b
! local variables
  real(kind=r8) :: stor1,stor2
! executable statements
  stor1 = a
  stor2 = b
  b = cs*stor2+sn*stor1
  a = cs*stor1-sn*stor2
end subroutine fprota
!-------------------------------------------------------
subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,ny,p,c,nc,&
  fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,ax,ay,bx,by,nrx,nry)
!  the b-spline coefficients of the smoothing spline are calculated as
!  the least-squares solution of the over-determined linear system of
!  equations  (ay) c (ax)' = q       where
!
!           |   (spx)    |            |   (spy)    |
!    (ax) = | ---------- |     (ay) = | ---------- |
!           | (1/p) (bx) |            | (1/p) (by) |
!
!                            | z  ' 0 |
!                        q = | ------ |
!                            | 0  ' 0 |
!
!  with c      : the (ny-ky-1) x (nx-kx-1) matrix which contains the
!            b-spline coefficients.
!   z      : the my x mx matrix which contains the function values.
!   spx,spy: the mx x (nx-kx-1) and  my x (ny-ky-1) observation
!            matrices according to the least-squares problems in
!            the x- and y-direction.
!   bx,by  : the (nx-2*kx-1) x (nx-kx-1) and (ny-2*ky-1) x (ny-ky-1)
!            matrices which contain the discontinuity jumps of the
!            derivatives of the b-splines in the x- and y-direction.
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: kx, kx1, kx2, ky, ky1, ky2, mm, mx, my, mz, &
                            mynx, nx, ny, nc
  integer, intent(inout) :: ifsx, ifsy, ifbx, ifby, nrx(mx), nry(my)
  real(kind=r8), intent(in)    :: p, x(mx), y(my), z(mz), tx(nx), ty(ny)
  real(kind=r8), intent(out)   :: fp, c(nc), q(mynx), ax(nx,kx2), ay(ny,ky2), fpx(nx), fpy(ny)
  real(kind=r8), intent(inout) :: spx(mx,kx1), spy(my,ky1), bx(nx,kx2), by(ny,ky2)
  real(kind=r8), intent(out)   :: right(mm)
! local variables
  real(kind=r8) :: arg, cs, fac, pinv, piv, sn, term
  integer       :: i, ibandx, ibandy, ic, iq, irot, it, iz,i1, i2, i3,  &
   j, k, k1, k2, l, l1, l2, ncof, nk1x, nk1y, nrold, nroldx, nroldy,    &
   number, numx, numx1, numy, numy1, n1
  real(kind=r8), dimension(7) :: h
! parameters
  real(kind=r8), parameter :: one = 1.0d0, half = 0.5d0
! external subroutines
  interface
    subroutine fpback(a,z,n,k,c,nest)
      use kinds, only : r8
      implicit none
      integer, intent(in)          :: n, k, nest
      real(kind=r8), intent(in)    :: a(nest,k), z(n)
      real(kind=r8), intent(inout) :: c(n)
    end subroutine fpback
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
    subroutine fpdisc(t,n,k2,b,nest)
      use kinds, only : r8
      implicit none
      integer      , intent(in)  :: n, k2, nest
      real(kind=r8), intent(in)  :: t(n)
      real(kind=r8), intent(out) :: b(nest,k2)
    end subroutine fpdisc
    subroutine fpgivs(piv,ww,cs,sn)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: piv
      real(kind=r8), intent(out)   :: cs, sn
      real(kind=r8), intent(inout) :: ww
    end subroutine fpgivs
    subroutine fprota(cs,sn,a,b)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: cs, sn
      real(kind=r8), intent(inout) :: a, b
    end subroutine fprota
  end interface
! executable statements
  nk1x = nx-kx1
  nk1y = ny-ky1
  if (p>0.0d0) pinv = one/p
!  it depends on the value of the flags ifsx,ifsy,ifbx and ifby and on
!  the value of p whether the matrices (spx),(spy),(bx) and (by) still
!  must be determined.
  if (ifsx==0) then
!  calculate the non-zero elements of the matrix (spx) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the x-direction.
    l = kx1
    l1 = kx2
    number = 0
    do it=1,mx
      arg = x(it)
      do while (arg>=tx(l1) .and. l/=nk1x)
        l = l1
        l1 = l+1
        number = number+1
      end do
      call fpbspl(tx,nx,kx,arg,l,h(:6))
      do i=1,kx1
        spx(it,i) = h(i)
      end do
      nrx(it) = number
    end do
    ifsx = 1
  end if
  if (ifsy==0) then
!  calculate the non-zero elements of the matrix (spy) which is the
!  observation matrix according to the least-squares spline approximat-
!  ion problem in the y-direction.
    l = ky1
    l1 = ky2
    number = 0
    do it=1,my
      arg = y(it)
      do while (arg>=ty(l1) .and. l/=nk1y)
        l = l1
        l1 = l+1
        number = number+1
      end do
      call fpbspl(ty,ny,ky,arg,l,h)
      do i=1,ky1
        spy(it,i) = h(i)
      end do
      nry(it) = number
    end do
    ifsy = 1
  end if
  if (p>0.0d0) then
!  calculate the non-zero elements of the matrix (bx).
    if (ifbx==0 .and. nx/=2*kx1) then
      call fpdisc(tx,nx,kx2,bx,nx)
      ifbx = 1
    end if
!  calculate the non-zero elements of the matrix (by).
    if (ifby==0 .and. ny/=2*ky1) then
      call fpdisc(ty,ny,ky2,by,ny)
      ifby = 1
    end if
!  reduce the matrix (ax) to upper triangular form (rx) using givens
!  rotations. apply the same transformations to the rows of matrix q
!  to obtain the my x (nx-kx-1) matrix g.
!  store matrix (rx) into (ax) and g into q.
  end if
  l = my*nk1x
!  initialization.
  do i=1,l
    q(i) = 0.0d0
  end do
  do i=1,nk1x
    do j=1,kx2
      ax(i,j) = 0.0d0
    end do
  end do
  l = 0
  nrold = 0
!  ibandx denotes the bandwidth of the matrices (ax) and (rx).
  ibandx = kx1
  do it=1,mx
    number = nrx(it)
    do
      if (nrold==number .or. p>0.0d0) then
        if(nrold==number) then
!  fetch a new row of matrix (spx).
          h(ibandx) = 0.0d0
          do j=1,kx1
            h(j) = spx(it,j)
          end do
!  find the appropriate column of q.
          do j=1,my
            l = l+1
            right(j) = z(l)
          end do
          irot = number
        else
          ibandx = kx2
!  fetch a new row of matrix (bx).
          n1 = nrold+1
          do j=1,kx2
            h(j) = bx(n1,j)*pinv
          end do
!  find the appropriate column of q.
          do j=1,my
            right(j) = 0.0d0
          end do
          irot = nrold
        end if
!  rotate the new row of matrix (ax) into triangle.
        do i=1,ibandx
          irot = irot+1
          piv = h(i)
          if (piv/=0.0d0) then
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,ax(irot,1),cs,sn)
!  apply that transformation to the rows of matrix q.
            iq = (irot-1)*my
            do j=1,my
              iq = iq+1
              call fprota(cs,sn,right(j),q(iq))
            end do
!  apply that transformation to the columns of (ax).
            if(i==ibandx) exit
            i2 = 1
            i3 = i+1
            do j=i3,ibandx
              i2 = i2+1
              call fprota(cs,sn,h(j),ax(irot,i2))
            end do
          end if
        end do
        if(nrold==number) exit
      end if
      nrold = nrold+1
    end do
  end do
!  reduce the matrix (ay) to upper triangular form (ry) using givens
!  rotations. apply the same transformations to the columns of matrix g
!  to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
!  store matrix (ry) into (ay) and h into c.
  ncof = nk1x*nk1y
!  initialization.
  do i=1,ncof
    c(i) = 0.0d0
  end do
  do i=1,nk1y
    do j=1,ky2
      ay(i,j) = 0.0d0
    end do
  end do
  nrold = 0
!  ibandy denotes the bandwidth of the matrices (ay) and (ry).
  ibandy = ky1
  do it=1,my
    number = nry(it)
    do
      if (nrold==number .or. p>0.0d0) then
        if (nrold==number) then
!  fetch a new row of matrix (spy)
          h(ibandy) = 0.0d0
          do j=1,ky1
            h(j) = spy(it,j)
          end do
!  find the appropiate row of g.
          l = it
          do j=1,nk1x
            right(j) = q(l)
            l = l+my
          end do
          irot = number
        else
          ibandy = ky2
!  fetch a new row of matrix (by).
          n1 = nrold+1
          do j=1,ky2
            h(j) = by(n1,j)*pinv
          end do
!  find the appropiate row of g.
          do j=1,nk1x
            right(j) = 0.0d0
          end do
          irot = nrold
        end if
!  rotate the new row of matrix (ay) into triangle.
        do i=1,ibandy
          irot = irot+1
          piv = h(i)
          if (piv/=0.0d0) then
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,ay(irot,1),cs,sn)
!  apply that transformation to the colums of matrix g.
            ic = irot
            do j=1,nk1x
              call fprota(cs,sn,right(j),c(ic))
              ic = ic+nk1y
            end do
!  apply that transformation to the columns of matrix (ay).
            if(i==ibandy) exit
            i2 = 1
            i3 = i+1
            do j=i3,ibandy
              i2 = i2+1
              call fprota(cs,sn,h(j),ay(irot,i2))
            end do
          end if
        end do
        if(nrold==number) exit
      end if
      nrold = nrold+1
    end do
  end do
!  backward substitution to obtain the b-spline coefficients as the
!  solution of the linear system    (ry) c (rx)' = h.
!  first step: solve the system  (ry) (c1) = h.
  k = 1
  do i=1,nk1x
    call fpback(ay,c(k),nk1y,ibandy,c(k),ny)
    k = k+nk1y
  end do
!  second step: solve the system  c (rx)' = (c1).
  k = 0
  do j=1,nk1y
    k = k+1
    l = k
    do i=1,nk1x
      right(i) = c(l)
      l = l+nk1y
    end do
    call fpback(ax,right,nk1x,ibandx,right,nx)
    l = k
    do i=1,nk1x
      c(l) = right(i)
      l = l+nk1y
    end do
  end do
!  calculate the quantities
!res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
!fp = sumi=1,mx(sumj=1,my(res(i,j)))
!fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
!              tx(r+kx) <= x(i) <= tx(r+kx+1)
!fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
!              ty(r+ky) <= y(j) <= ty(r+ky+1)
  fp = 0.0d0
  do i=1,nx
    fpx(i) = 0.0d0
  end do
  do i=1,ny
    fpy(i) = 0.0d0
  end do
  nk1y = ny-ky1
  iz = 0
  nroldx = 0
!  main loop for the different grid points.
  do i1=1,mx
    numx = nrx(i1)
    numx1 = numx+1
    nroldy = 0
    do i2=1,my
      numy = nry(i2)
      numy1 = numy+1
      iz = iz+1
!  evaluate s(x,y) at the current grid point by making the sum of the
!  cross products of the non-zero b-splines at (x,y), multiplied with
!  the appropiate b-spline coefficients.
      term = 0.0d0
      k1 = numx*nk1y+numy
      do l1=1,kx1
        k2 = k1
        fac = spx(i1,l1)
        do l2=1,ky1
          k2 = k2+1
          term = term+fac*spy(i2,l2)*c(k2)
        end do
        k1 = k1+nk1y
      end do
!  calculate the squared residual at the current grid point.
      term = (z(iz)-term)**2
!  adjust the different parameters.
      fp = fp+term
      fpx(numx1) = fpx(numx1)+term
      fpy(numy1) = fpy(numy1)+term
      fac = term*half
      if(numy/=nroldy) then
        fpy(numy1) = fpy(numy1)-fac
        fpy(numy) = fpy(numy)+fac
      end if
      nroldy = numy
      if(numx/=nroldx) then
        fpx(numx1) = fpx(numx1)-fac
        fpx(numx) = fpx(numx)+fac
      end if
    end do
    nroldx = numx
  end do
end subroutine fpgrre
!-------------------------------------------------------
subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
!  subroutine fpknot locates an additional knot for a spline of degree
!  k and adjusts the corresponding parameters,i.e.
!t     : the position of the knots.
!n     : the number of knots.
!nrint : the number of knotintervals.
!fpint : the sum of squares of residual right hand sides
!        for each knot interval.
!nrdata: the number of data points inside each knot interval.
!  istart indicates that the smallest data point at which the new knot
!  may be added is x(istart+1)
  use kinds, only : r8
  implicit none
! arguments
  integer      , intent(in)    :: m, nest, istart
  integer      , intent(inout) :: n, nrint, nrdata(nest)
  real(kind=r8), intent(in)    :: x(m)
  real(kind=r8), intent(inout) :: t(nest), fpint(nest)
! local variables
  real(kind=r8) :: an, am, fpmax
  integer       :: ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg, maxpt,  &
                   next, nrx, number
! executable statements
  k = (n-nrint-1)/2
!  search for knot interval t(number+k) <= x <= t(number+k+1) where
!  fpint(number) is maximal on the condition that nrdata(number)
!  not equals zero.
  fpmax = 0.0d0
  jbegin = istart
  do j=1,nrint
    jpoint = nrdata(j)
    if(fpmax<fpint(j) .and. jpoint/=0) then
      fpmax = fpint(j)
      number = j
      maxpt = jpoint
      maxbeg = jbegin
    end if
    jbegin = jbegin+jpoint+1
  end do
!  let coincide the new knot t(number+k+1) with a data point x(nrx)
!  inside the old knot interval t(number+k) <= x <= t(number+k+1).
  ihalf = maxpt/2+1
  nrx = maxbeg+ihalf
  next = number+1
!  adjust the different parameters.
  do j=next,nrint
    jj = next+nrint-j
    fpint(jj+1) = fpint(jj)
    nrdata(jj+1) = nrdata(jj)
    jk = jj+k
    t(jk+1) = t(jk)
  end do
  nrdata(number) = ihalf-1
  nrdata(next) = maxpt-ihalf
  am = maxpt
  an = nrdata(number)
  fpint(number) = fpmax*an/am
  an = nrdata(next)
  fpint(next) = fpmax*an/am
  jk = next+k
  t(jk) = x(nrx)
  n = n+1
  nrint = nrint+1
end subroutine fpknot
!-------------------------------------------------------
function fprati(p1,f1,p2,f2,p3,f3)
!  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
!  gives the value of p such that the rational interpolating function
!  of the form r(p) = (u*p+v)/(p+w) equals zero at p.
  use kinds, only : r8
  implicit none
  real(kind=r8) fprati
! arguments
  real(kind=r8), intent(in)    :: p2, f2
  real(kind=r8), intent(inout) :: p1, f1, p3, f3
! local variables
  real(kind=r8) :: h1, h2, h3, p
! executable statements
  if(p3>0.0d0) then
!  value of p in case p3 ^= infinity.
    h1 = f1*(f2-f3)
    h2 = f2*(f3-f1)
    h3 = f3*(f1-f2)
    p = -(p1*p2*h3+p2*p3*h1+p3*p1*h2)/(p1*h1+p2*h2+p3*h3)
  else
!  value of p in case p3 = infinity.
    p = (p1*(f1-f3)*f2-p2*(f2-f3)*f1)/((f1-f2)*f3)
  end if
!  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
  if(f2<0.0d0) then
    p3 = p2
    f3 = f2
  else
    p1 = p2
    f1 = f2
  end if
  fprati = p
end function fprati
!-------------------------------------------------------
subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,tol,&
   maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,reducy,fpintx,fpinty,lastdi,&
   nplusx,nplusy,nrx,nry,nrdatx,nrdaty,wrk,lwrk,ier)
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: kx,ky,mx,my,mz,nxest,nyest,maxit,nc,lwrk,iopt
  integer, intent(inout) :: nx,ny,lastdi,nplusx,nplusy,                   &
                            nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
  integer, intent(out)   :: ier
  real(kind=r8), intent(in)    :: xb,xe,yb,ye,s,tol,x(mx),y(my),z(mz)
  real(kind=r8), intent(out)   :: fp,c(nc)
  real(kind=r8), intent(inout) :: fp0,fpold,reducx,reducy,tx(nxest),      &
                    ty(nyest),fpintx(nxest),fpinty(nyest),wrk(lwrk)
! local variables
  real(kind=r8) :: acc,fpms,f1,f2,f3,p,p1,p2,p3,rn
  integer :: i,ich1,ich3,ifbx,ifby,ifsx,ifsy,iter,j,kx1,kx2,ky1,ky2,k3,   &
   l,lax,lay,lbx,lby,lq,lri,lsx,lsy,mk1,mm,mpm,mynx,                      &
   nmaxx,nmaxy,nminx,nminy,nplx,nply,npl1,nrintx,nrinty,nxe,nxk,nye
! parameters
  real(kind=r8), parameter :: one=1.0d0, half=0.5d0, con1=0.1d0,          &
                              con9=0.9d0, con4=0.4d-1
! intrinsic functions
!  integer, intrinsic :: max, min
!  real, intrinsic :: abs
! external functions/subroutines
  interface
    subroutine fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,ny,p,c,nc,&
      fp,fpx,fpy,mm,mynx,kx1,kx2,ky1,ky2,spx,spy,right,q,ax,ay,bx,by,nrx,nry)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: kx, kx1, kx2, ky, ky1, ky2, mm, mx, my, mz, &
                                mynx, nx, ny, nc
      integer, intent(inout) :: ifsx, ifsy, ifbx, ifby, nrx(mx), nry(my)
      real(kind=r8), intent(in)    :: p, x(mx), y(my), z(mz), tx(nx), ty(ny)
      real(kind=r8), intent(out)   :: fp, c(nc), q(mynx), ax(nx,kx2), ay(ny,ky2), fpx(nx), fpy(ny)
      real(kind=r8), intent(inout) :: spx(mx,kx1), spy(my,ky1), bx(nx,kx2), by(ny,ky2)
      real(kind=r8), intent(out)   :: right(mm)
    end subroutine fpgrre
    subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: m, nest, istart
      integer      , intent(inout) :: n, nrint, nrdata(nest)
      real(kind=r8), intent(in)    :: x(m)
      real(kind=r8), intent(inout) :: t(nest), fpint(nest)
    end subroutine fpknot
    function fprati(p1,f1,p2,f2,p3,f3)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fprati
      real(kind=r8), intent(in)    :: p2, f2
      real(kind=r8), intent(inout) :: p1, f1, p3, f3
    end function fprati
  end interface
! executable statements
!  we partition the working space.
  kx1 = kx+1
  ky1 = ky+1
  kx2 = kx1+1
  ky2 = ky1+1
  lsx = 1
  lsy = lsx+mx*kx1
  lri = lsy+my*ky1
  mm = max(nxest,my)
  lq = lri+mm
  mynx = nxest*my
  lax = lq+mynx
  nxk = nxest*kx2
  lbx = lax+nxk
  lay = lbx+nxk
  lby = lay+nyest*ky2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
!  given a set of knots we compute the least-squares spline sinf(x,y), c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1  sinf(x,y) is the requested approximation.               c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!if fp <=s we will continue with the current set of knots.             c
!if fp > s we will increase the number of knots and compute the        c
!   corresponding least-squares spline until finally fp<=s.            c
!the initial choice of knots depends on the value of s and iopt.       c
!if s=0 we have spline interpolation; in that case the number of       c
!knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.                   c
!if s>0 and                                                            c
! *iopt=0 we first compute the least-squares polynomial of degree      c
!  kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.           c
! *iopt=1 we start with the knots found at the last call of the        c
!  routine, except for the case that s > fp0; then we can compute      c
!  the least-squares polynomial directly.                              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine the number of knots for polynomial approximation.
  nminx = 2*kx1
  nminy = 2*ky1
  if (iopt>=0) then
!  acc denotes the absolute tolerance for the root of f(p)=s.
    acc = tol*s
!  find nmaxx and nmaxy which denote the number of knots in x- and y-
!  direction in case of spline interpolation.
    nmaxx = mx+kx1
    nmaxy = my+ky1
!  find nxe and nye which denote the maximum number of knots
!  allowed in each direction
    nxe = min(nmaxx,nxest)
    nye = min(nmaxy,nyest)
    if (s>0.0d0) then
!  if s > 0 our initial choice of knots depends on the value of iopt.
      if (iopt==0 .or. fp0<=s) then
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial of degree kx in x and ky in y (which is a spline without
!  interior knots).
        nx = nminx
        ny = nminy
        nrdatx(1) = mx-2
        nrdaty(1) = my-2
        lastdi = 0
        nplusx = 0
        nplusy = 0
        fp0 = 0.0d0
        fpold = 0.0d0
        reducx = 0.0d0
        reducy = 0.0d0
      else
!  if iopt=1 and fp0 > s we start computing the least- squares spline
!  according to the set of knots found at the last call of the routine.
!  we determine the number of grid coordinates x(i) inside each knot
!  interval (tx(l),tx(l+1)).
        l = kx2
        j = 1
        nrdatx(1) = 0
        mpm = mx-1
        do i=2,mpm
          nrdatx(j) = nrdatx(j)+1
          if (x(i)>=tx(l)) then
            nrdatx(j) = nrdatx(j)-1
            l = l+1
            j = j+1
            nrdatx(j) = 0
          end if
        end do
!  we determine the number of grid coordinates y(i) inside each knot
!  interval (ty(l),ty(l+1)).
        l = ky2
        j = 1
        nrdaty(1) = 0
        mpm = my-1
        do i=2,mpm
          nrdaty(j) = nrdaty(j)+1
          if (y(i)>=ty(l)) then
            nrdaty(j) = nrdaty(j)-1
            l = l+1
            j = j+1
            nrdaty(j) = 0
          end if
        end do
      end if
    else
!  if s = 0, s(x,y) is an interpolating spline.
      nx = nmaxx
      ny = nmaxy
!  test whether the required storage space exceeds the available one.
      if (ny>nyest .or. nx>nxest) then
        ier = 1
        return
      end if
!  find the position of the interior knots in case of interpolation.
!  the knots in the x-direction.
      mk1 = mx-kx1
      if (mk1/=0) then
        k3 = kx/2
        i = kx1+1
        j = k3+2
        if (k3*2==kx) then
          do l=1,mk1
            tx(i) = (x(j)+x(j-1))*half
            i = i+1
            j = j+1
          end do
        else
          do l=1,mk1
            tx(i) = x(j)
            i = i+1
            j = j+1
          end do
        end if
      end if
!  the knots in the y-direction.
      mk1 = my-ky1
      if (mk1/=0) then
        k3 = ky/2
        i = ky1+1
        j = k3+2
        if (k3*2==ky) then
          do l=1,mk1
            ty(i) = (y(j)+y(j-1))*half
            i = i+1
            j = j+1
          end do
        else
          do l=1,mk1
            ty(i) = y(j)
            i = i+1
            j = j+1
          end do
        end if
      end if
    end if
  end if
  mpm = mx+my
  ifsx = 0
  ifsy = 0
  ifbx = 0
  ifby = 0
  p = -one
!  main loop for the different sets of knots.mpm=mx+my is a save upper
!  bound for the number of trials.
  do iter=1,mpm
    if (nx==nminx .and. ny==nminy) ier = -2
!  find nrintx (nrinty) which is the number of knot intervals in the
!  x-direction (y-direction).
    nrintx = nx-nminx+1
    nrinty = ny-nminy+1
!  find ncof, the number of b-spline coefficients for the current set
!  of knots.
!    nk1x = nx-kx1
!    nk1y = ny-ky1
!    ncof = nk1x*nk1y
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(x,y).
    i = nx
    do j=1,kx1
      tx(j) = xb
      tx(i) = xe
      i = i-1
    end do
    i = ny
    do j=1,ky1
      ty(j) = yb
      ty(i) = ye
      i = i-1
    end do
!  find the least-squares spline sinf(x,y) and calculate for each knot
!  interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
!  of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
!  ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
!  value belonging to that interval.
!  fp gives the total sum of squared residuals.
    call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,&
    ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx:lsy-1),&
    wrk(lsy:lri-1),wrk(lri:lq-1),wrk(lq:lax-1),wrk(lax:lbx-1),&
    wrk(lay:lby-1),wrk(lbx:lay-1),wrk(lby:),nrx,nry)
    if (ier==(-2)) fp0 = fp
!  test whether the least-squares spline is an acceptable solution.
    if (iopt<0) return
    fpms = fp-s
    if (abs(fpms) < acc) return
!  if f(p=inf) < s, we accept the choice of knots.
    if (fpms<0.0d0) exit
!  if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
    if (nx==nmaxx .and. ny==nmaxy) then
      ier = -1
      fp = 0.0d0
      return
    end if
!  increase the number of knots.
!  if nx=nxe and ny=nye we cannot further increase the number of knots
!  because of the storage capacity limitation.
    if (nx==nxe .and. ny==nye) then
      ier = 1
      return
    end if
    ier = 0
!  adjust the parameter reducx or reducy according to the direction
!  in which the last added knots were located.
    if (lastdi<0) then
      reducx = fpold-fp
    else if (lastdi>0) then
      reducy = fpold-fp
    end if
!  store the sum of squared residuals for the current set of knots.
    fpold = fp
!  find nplx, the number of knots we should add in the x-direction.
    nplx = 1
    if (nx/=nminx) then
      npl1 = nplusx*2
      rn = nplusx
      if (reducx>acc) npl1 = rn*fpms/reducx
      nplx = min(nplusx*2,max(npl1,nplusx/2,1))
    end if
!  find nply, the number of knots we should add in the y-direction.
    nply = 1
    if (ny/=nminy) then
      npl1 = nplusy*2
      rn = nplusy
      if (reducy>acc) npl1 = rn*fpms/reducy
      nply = min(nplusy*2,max(npl1,nplusy/2,1))
    end if
    if (ny==nye .or. (nx/=nxe .and. ny/=nye .and.&
        (nplx<nply .or. (nplx==nply .and. lastdi>=0)))) then
!  addition in the x-direction.
      lastdi = -1
      nplusx = nplx
      ifsx = 0
      do l=1,nplusx
!  add a new knot in the x-direction
        call fpknot(x,mx,tx,nx,fpintx,nrdatx,nrintx,nxest,1)
!  test whether we cannot further increase the number of knots in the
!  x-direction.
        if (nx==nxe) exit
      end do
    else
!  addition in the y-direction.
      lastdi = 1
      nplusy = nply
      ifsy = 0
      do l=1,nplusy
!  add a new knot in the y-direction.
        call fpknot(y,my,ty,ny,fpinty,nrdaty,nrinty,nyest,1)
!  test whether we cannot further increase the number of knots in the
!  y-direction.
        if (ny==nye) exit
      end do
    end if
!  restart the computations with the new set of knots.
  end do
!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
  if (ier==(-2)) return
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(x,y)                c
! *****************************************************                c
!  we have determined the number of knots and their position. we now   c
!  compute the b-spline coefficients of the smoothing spline sp(x,y).  c
!  this smoothing spline varies with the parameter p in such a way thatc
!f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)                 c
!  is a continuous, strictly decreasing function of p. moreover the    c
!  least-squares polynomial corresponds to p=0 and the least-squares   c
!  spline to p=infinity. iteratively we then have to determine the     c
!  positive value of p such that f(p)=s. the process which is proposed c
!  here makes use of rational interpolation. f(p) is approximated by a c
!  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  c
!  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)c
!  are used to calculate the new value of p such that r(p)=s.          c
!  convergence is guaranteed by taking f1 > 0 and f3 < 0.              c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  initial value for p.
  p1 = 0.0d0
  f1 = fp0-s
  p3 = -one
  f3 = fpms
  p = one
  ich1 = 0
  ich3 = 0
!  iteration process to find the root of f(p)=s.
  do iter = 1,maxit
!  find the smoothing spline sp(x,y) and the corresponding sum of
!  squared residuals fp.
    call fpgrre(ifsx,ifsy,ifbx,ifby,x,mx,y,my,z,mz,kx,ky,tx,nx,ty,&
    ny,p,c,nc,fp,fpintx,fpinty,mm,mynx,kx1,kx2,ky1,ky2,wrk(lsx:lsy-1),&
    wrk(lsy:lri-1),wrk(lri:lq-1),wrk(lq:lax-1),wrk(lax:lbx-1),&
    wrk(lay:lby-1),wrk(lbx:lay-1),wrk(lby:),nrx,nry)
!  test whether the approximation sp(x,y) is an acceptable solution.
    fpms = fp-s
    if(abs(fpms)<acc) return
!  test whether the maximum allowable number of iterations has been
!  reached.
    if(iter==maxit) exit
!  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3==0) then
      if ((f2-f3)<=acc) then
!  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if(p<=p1) p = p1*con9 + p2*con1
        cycle
      end if
      if (f2<0.0d0) ich3 = 1
    end if
    if (ich1==0) then
      if ((f1-f2)<=acc) then
!  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if (p3<0.0d0) cycle
        if (p>=p3) p = p2*con1 + p3*con9
        cycle
      end if
!  test whether the iteration process proceeds as theoretically
!  expected.
      if (f2>0.0d0) ich1 = 1
    end if
    if (f2>=f1 .or. f2<=f3) then
      ier = 2
      return
    end if
!  find the new value of p.
    p = fprati(p1,f1,p2,f2,p3,f3)
  end do
  ier = 3
end subroutine fpregr
!-------------------------------------------------------
subroutine fpchec(x,m,t,n,k,ier)
!  subroutine fpchec verifies the number and the position of the knots
!  t(j),j=1,2,...,n of a spline of degree k, in relation to the number
!  and the position of the data points x(i),i=1,2,...,m. if all of the
!  following conditions are fulfilled, the error parameter ier is set
!  to zero. if one of the conditions is violated ier is set to ten.
!  1) k+1 <= n-k-1 <= m
!  2) t(1) <= t(2) <= ... <= t(k+1)
!     t(n-k) <= t(n-k+1) <= ... <= t(n)
!  3) t(k+1) < t(k+2) < ... < t(n-k)
!  4) t(k+1) <= x(i) <= t(n-k)
!  5) the conditions specified by schoenberg and whitney must hold
!     for at least one subset of data points, i.e. there must be a
!     subset of data points y(j) such that
!         t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1
!  ..
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)       :: m, n, k
  real(kind=r8), intent(in) :: x(m), t(n)
  integer, intent(out)      :: ier
! local variables
  integer       :: i, j, k1, k2, l, nk1, nk2, nk3
  real(kind=r8) :: tj, tl
! executable statements
  k1 = k+1
  k2 = k1+1
  nk1 = n-k1
  nk2 = nk1+1
  ier = 10
!  check condition no 1
  if (nk1<k1 .or. nk1>m) return
!  check condition no 2
  j = n
  do i=1,k
    if (t(i)>t(i+1) .or. t(j)<t(j-1)) return
    j = j-1
  end do
!  check condition no 3
  do i=k2,nk2
    if (t(i)<=t(i-1)) return
  end do
!  check condition no 4
  if (x(1)<t(k1) .or. x(m)>t(nk2)) return
!  check condition no 5
  if (x(1)>=t(k2) .or. x(m)<=t(nk1)) return
  i = 1
  l = k2
  nk3 = nk1-1
  do j=2,nk3
    tj = t(j)
    l = l+1
    tl = t(l)
    do
      i = i+1
      if (i>=m) return
      if (x(i)>tj) exit
    end do
    if (x(i)>=tl) return
  end do
  ier = 0
end subroutine fpchec
!-------------------------------------------------------
  subroutine scatterspl(x,y,z,w,n,kspl,sspl,xmin,xmax,ymin,ymax, &
                        tx,nknt_x,ty,nknt_y,coeff,ierr)
    use kinds, only : r8
    implicit none
! arguments
    integer, intent(in) :: n
    real(kind=r8), dimension(n), intent(in) :: x, y, z
    real(kind=r8), dimension(n), intent(in) :: w
    integer, intent(in) :: kspl
    real(kind=r8), intent(in) :: sspl
    real(kind=r8), intent(in) :: xmin, xmax, ymin, ymax
    real(kind=r8), dimension(nknt_x), intent(inout) :: tx
    real(kind=r8), dimension(nknt_y), intent(inout) :: ty
    integer, intent(inout) :: nknt_x, nknt_y
    real(kind=r8), dimension(nknt_x*nknt_y), intent(out) :: coeff
    integer, intent(out) :: ierr
! local variables
    integer :: iopt
    real(kind=r8) :: resid
    integer :: u,v,km,ne,b1,b2,lwrk1,lwrk2,kwrk,nxest,nyest
    real(kind=r8), dimension(:), allocatable :: wrk1, wrk2
    integer, dimension(:), allocatable :: iwrk

    interface
      subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
          nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
        use kinds, only : r8
        implicit none
        real(kind=r8) xb,xe,yb,ye,s,eps,fp
        integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
        real(kind=r8) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax), &
          c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
        integer iwrk(kwrk)
      end subroutine surfit
    end interface

    nxest=nknt_x
    nyest=nknt_y
    ne = max(nxest,nyest)

    km = kspl+1
    u = nxest-km
    v = nyest-km
    b1 = kspl*min(u,v)+kspl+1
    b2 = (kspl+1)*min(u,v)+1
    lwrk1 = u*v*(2+b1+b2)+2*(u+v+km*(n+ne)+ne-2*kspl)+b2+1
    lwrk2 = u*v*(b2+1)+b2
    kwrk = n+(nknt_x-2*kspl-1)*(nknt_y-2*kspl-1)
    allocate(wrk1(lwrk1),wrk2(lwrk2),iwrk(kwrk))

    iopt=0
    call surfit(iopt,n,x,y,z,w,xmin,xmax,ymin,ymax,kspl,kspl, &
        sspl,nxest,nyest,ne,EPSILON(1.0d0),nknt_x,tx,nknt_y,ty, &
        coeff,resid,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ierr)

    deallocate(wrk1,wrk2,iwrk)

  end subroutine scatterspl
!-------------------------------------------------------
  subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
      nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
! given the set of data points (x(i),y(i),z(i)) and the set of positive
! numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
! iate spline approximation s(x,y) of degrees kx and ky on the rect-
! angle xb <= x <= xe, yb <= y <= ye.
! if iopt = -1 surfit calculates the weighted least-squares spline
! according to a given set of knots.
! if iopt >= 0 the total numbers nx and ny of these knots and their
! position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
! ally by the routine. the smoothness of s(x,y) is then achieved by
! minimalizing the discontinuity jumps in the derivatives of s(x,y)
! across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
! the amounth of smoothness is determined by the condition that f(p) =
! sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
! ative constant, called the smoothing factor.
! the fit is given in the b-spline representation (b-spline coefficients
! c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
! uated by means of subroutine bispev.
!
! calling sequence:
! call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
!    *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
!
! parameters:
!  iopt  : integer flag. on entry iopt must specify whether a weighted
!      least-squares spline (iopt=-1) or a smoothing spline (iopt=0
!      or 1) must be determined.
!      if iopt=0 the routine will start with an initial set of knots
!      tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
!      1,...,ky+1. if iopt=1 the routine will continue with the set
!      of knots found at the last call of the routine.
!      attention: a call with iopt=1 must always be immediately pre-
!                 ceded by another call with iopt=1 or iopt=0.
!      unchanged on exit.
!  m     : integer. on entry m must specify the number of data points.
!      m >= (kx+1)*(ky+1). unchanged on exit.
!  x     : real array of dimension at least (m).
!  y     : real array of dimension at least (m).
!  z     : real array of dimension at least (m).
!      before entry, x(i),y(i),z(i) must be set to the co-ordinates
!      of the i-th data point, for i=1,...,m. the order of the data
!      points is immaterial. unchanged on exit.
!  w     : real array of dimension at least (m). before entry, w(i) must
!      be set to the i-th value in the set of weights. the w(i) must
!      be strictly positive. unchanged on exit.
!  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
!  yb,ye   aries of the rectangular approximation domain.
!      xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
!  kx,ky : integer values. on entry kx and ky must specify the degrees
!      of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
!      (kx=ky=3) splines. unchanged on exit.
!  s     : real. on entry (in case iopt>=0) s must specify the smoothing
!      factor. s >=0. unchanged on exit.
!      for advice on the choice of s see further comments
!  nxest : integer. unchanged on exit.
!  nyest : integer. unchanged on exit.
!      on entry, nxest and nyest must specify an upper bound for the
!      number of knots required in the x- and y-directions respect.
!      these numbers will also determine the storage space needed by
!      the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
!      in most practical situation nxest = kx+1+sqrt(m/2), nyest =
!      ky+1+sqrt(m/2) will be sufficient. see also further comments.
!  nmax  : integer. on entry nmax must specify the actual dimension of
!      the arrays tx and ty. nmax >= nxest, nmax >=nyest.
!      unchanged on exit.
!  eps   : real.
!      on entry, eps must specify a threshold for determining the
!      effective rank of an over-determined linear system of equat-
!      ions. 0 < eps < 1.  if the number of decimal digits in the
!      computer representation of a real number is q, then 10**(-q)
!      is a suitable value for eps in most practical applications.
!      unchanged on exit.
!  nx    : integer.
!      unless ier=10 (in case iopt >=0), nx will contain the total
!      number of knots with respect to the x-variable, of the spline
!      approximation returned. if the computation mode iopt=1 is
!      used, the value of nx should be left unchanged between sub-
!      sequent calls.
!      in case iopt=-1, the value of nx should be specified on entry
!  tx    : real array of dimension nmax.
!      on succesful exit, this array will contain the knots of the
!      spline with respect to the x-variable, i.e. the position of
!      the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
!      position of the additional knots tx(1)=...=tx(kx+1)=xb and
!      tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
!      if the computation mode iopt=1 is used, the values of tx(1),
!      ...,tx(nx) should be left unchanged between subsequent calls.
!      if the computation mode iopt=-1 is used, the values tx(kx+2),
!      ...tx(nx-kx-1) must be supplied by the user, before entry.
!      see also the restrictions (ier=10).
!  ny    : integer.
!      unless ier=10 (in case iopt >=0), ny will contain the total
!      number of knots with respect to the y-variable, of the spline
!      approximation returned. if the computation mode iopt=1 is
!      used, the value of ny should be left unchanged between sub-
!      sequent calls.
!      in case iopt=-1, the value of ny should be specified on entry
!  ty    : real array of dimension nmax.
!      on succesful exit, this array will contain the knots of the
!      spline with respect to the y-variable, i.e. the position of
!      the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
!      position of the additional knots ty(1)=...=ty(ky+1)=yb and
!      ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
!      if the computation mode iopt=1 is used, the values of ty(1),
!      ...,ty(ny) should be left unchanged between subsequent calls.
!      if the computation mode iopt=-1 is used, the values ty(ky+2),
!      ...ty(ny-ky-1) must be supplied by the user, before entry.
!      see also the restrictions (ier=10).
!  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
!      on succesful exit, c contains the coefficients of the spline
!      approximation s(x,y)
!  fp    : real. unless ier=10, fp contains the weighted sum of
!      squared residuals of the spline approximation returned.
!  wrk1  : real array of dimension (lwrk1). used as workspace.
!      if the computation mode iopt=1 is used the value of wrk1(1)
!      should be left unchanged between subsequent calls.
!      on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
!      contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
!      with d(i) the i-th diagonal element of the reduced triangular
!      matrix for calculating the b-spline coefficients. it includes
!      those elements whose square is less than eps,which are treat-
!      ed as 0 in the case of presumed rank deficiency (ier<-2).
!  lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
!      the array wrk1 as declared in the calling (sub)program.
!      lwrk1 must not be too small. let
!        u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
!        ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
!        if(bx<=by) b1 = bx, b2 = b1+v-ky
!        if(bx>by) b1 = by, b2 = b1+u-kx  then
!      lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
!  wrk2  : real array of dimension (lwrk2). used as workspace, but
!      only in the case a rank deficient system is encountered.
!  lwrk2 : integer. on entry lwrk2 must specify the actual dimension of
!      the array wrk2 as declared in the calling (sub)program.
!      lwrk2 > 0 . a save upper boundfor lwrk2 = u*v*(b2+1)+b2
!      where u,v and b2 are as above. if there are enough data
!      points, scattered uniformly over the approximation domain
!      and if the smoothing factor s is not too small, there is a
!      good chance that this extra workspace is not needed. a lot
!      of memory might therefore be saved by setting lwrk2=1.
!      (see also ier > 10)
!  iwrk  : integer array of dimension (kwrk). used as workspace.
!  kwrk  : integer. on entry kwrk must specify the actual dimension of
!      the array iwrk as declared in the calling (sub)program.
!      kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
!  ier   : integer. unless the routine detects an error, ier contains a
!      non-positive value on exit, i.e.
!   ier=0  : normal return. the spline returned has a residual sum of
!        squares fp such that abs(fp-s)/s <= tol with tol a relat-
!        ive tolerance set to 0.001 by the program.
!   ier=-1 : normal return. the spline returned is an interpolating
!        spline (fp=0).
!   ier=-2 : normal return. the spline returned is the weighted least-
!        squares polynomial of degrees kx and ky. in this extreme
!        case fp gives the upper bound for the smoothing factor s.
!   ier<-2 : warning. the coefficients of the spline returned have been
!        computed as the minimal norm least-squares solution of a
!        (numerically) rank deficient system. (-ier) gives the rank.
!        especially if the rank deficiency which can be computed as
!        (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
!        curate. they could also seriously depend on the value of
!        eps.
!   ier=1  : error. the required storage space exceeds the available
!        storage space, as specified by the parameters nxest and
!        nyest.
!        probably causes : nxest or nyest too small. if these param-
!        eters are already large, it may also indicate that s is
!        too small
!        the approximation returned is the weighted least-squares
!        spline according to the current set of knots.
!        the parameter fp gives the corresponding weighted sum of
!        squared residuals (fp>s).
!   ier=2  : error. a theoretically impossible result was found during
!        the iteration proces for finding a smoothing spline with
!        fp = s. probably causes : s too small or badly chosen eps.
!        there is an approximation returned but the corresponding
!        weighted sum of squared residuals does not satisfy the
!        condition abs(fp-s)/s < tol.
!   ier=3  : error. the maximal number of iterations maxit (set to 20
!        by the program) allowed for finding a smoothing spline
!        with fp=s has been reached. probably causes : s too small
!        there is an approximation returned but the corresponding
!        weighted sum of squared residuals does not satisfy the
!        condition abs(fp-s)/s < tol.
!   ier=4  : error. no more knots can be added because the number of
!        b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
!        the number of data points m.
!        probably causes : either s or m too small.
!        the approximation returned is the weighted least-squares
!        spline according to the current set of knots.
!        the parameter fp gives the corresponding weighted sum of
!        squared residuals (fp>s).
!   ier=5  : error. no more knots can be added because the additional
!        knot would (quasi) coincide with an old one.
!        probably causes : s too small or too large a weight to an
!        inaccurate data point.
!        the approximation returned is the weighted least-squares
!        spline according to the current set of knots.
!        the parameter fp gives the corresponding weighted sum of
!        squared residuals (fp>s).
!   ier=10 : error. on entry, the input data are controlled on validity
!        the following restrictions must be satisfied.
!        -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
!        nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
!        xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
!        lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
!        kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
!        if iopt=-1: 2*kx+2<=nx<=nxest
!                    xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
!                    2*ky+2<=ny<=nyest
!                    yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
!        if iopt>=0: s>=0
!        if one of these conditions is found to be violated,control
!        is immediately repassed to the calling program. in that
!        case there is no approximation returned.
!   ier>10 : error. lwrk2 is too small, i.e. there is not enough work-
!        space for computing the minimal least-squares solution of
!        a rank deficient system of linear equations. ier gives the
!        requested value for lwrk2. there is no approximation re-
!        turned but, having saved the information contained in nx,
!        ny,tx,ty,wrk1, and having adjusted the value of lwrk2 and
!        the dimension of the array wrk2 accordingly, the user can
!        continue at the point the program was left, by calling
!        surfit with iopt=1.
!
! further comments:
!  by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating spline if
!   s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in z(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   polynomial and the corresponding upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximation shows more detail) to obtain closer fits.
!   to choose s very small is strongly discouraged. this considerably
!   increases computation time and memory requirements. it may also
!   cause rank-deficiency (ier<-2) and endager numerical stability.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if surfit is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. if the computation mode iopt=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   surfit once more with the selected value for s but now with iopt=0.
!   indeed, surfit may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!   the number of knots may also depend on the upper bounds nxest and
!   nyest. indeed, if at a certain stage in surfit the number of knots
!   in one direction (say nx) has reached the value of its upper bound
!   (nxest), then from that moment on all subsequent knots are added
!   in the other (y) direction. this may indicate that the value of
!   nxest is too small. on the other hand, it gives the user the option
!   of limiting the number of knots the routine locates in any direction
!   for example, by setting nxest=2*kx+2 (the lowest allowable value for
!   nxest), the user can indicate that he wants an approximation which
!   is a simple polynomial of degree kx in the variable x.
!
!  other subroutines required:
!    fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
!    fprank,fporde
!
!  references:
!   dierckx p. : an algorithm for surface fitting with spline functions
!            ima j. numer. anal. 1 (1981) 267-283.
!   dierckx p. : an algorithm for surface fitting with spline functions
!            report tw50, dept. computer science,k.u.leuven, 1980.
!   dierckx p. : curve and surface fitting with splines, monographs on
!            numerical analysis, oxford university press, 1993.
!
!  author:
!    p.dierckx
!    dept. computer science, k.u. leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : may 1979
!  latest update : march 1987
!
!  ..
    use kinds, only : r8
    implicit none
!  ..scalar arguments..
    real(kind=r8) xb,xe,yb,ye,s,eps,fp
    integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
!  ..array arguments..
    real(kind=r8) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax), &
     c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
    integer iwrk(kwrk)
!  ..local scalars..
    real(kind=r8) tol
    integer i,ib1,ib3,jb1,ki,kmax,km1,km2,kn,kwest,kx1,ky1,la,lbx, &
     lby,lco,lf,lff,lfp,lh,lq,lsx,lsy,lwest,maxit,ncest,nest,nek, &
     nminx,nminy,nmx,nmy,nreg,nrint,nxk,nyk
!  ..function references..
    integer max0
!  ..subroutine references..
!    fpsurf
    interface
      subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest, &
         nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest, &
         nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h, &
         idx,nummer,wrk,lwrk,ier)
        use kinds, only : r8
        implicit none
        real(kind=r8) xb,xe,yb,ye,s,eta,tol,fp,fp0
        integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3, &
         nc,intest,nrest,nx0,ny0,lwrk,ier
        real(kind=r8) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest), &
         coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2), &
         by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
        integer idx(nrest),nummer(m)
      end subroutine fpsurf
    end interface
!  ..
!  we set up the parameters tol and maxit.
    maxit = 20
    tol = 0.1e-02
!  before starting computations a data check is made. if the input data
!  are invalid,control is immediately repassed to the calling program.
    ier = 10
    if(eps<=0. .or. eps>=1.) return
    if(kx<=0 .or. kx>5) return
    kx1 = kx+1
    if(ky<=0 .or. ky>5) return
    ky1 = ky+1
    kmax = max0(kx,ky)
    km1 = kmax+1
    km2 = km1+1
    if(iopt<(-1) .or. iopt>1) return
    if(m<(kx1*ky1)) return
    nminx = 2*kx1
    if(nxest<nminx .or. nxest>nmax) return
    nminy = 2*ky1
    if(nyest<nminy .or. nyest>nmax) return
    nest = max0(nxest,nyest)
    nxk = nxest-kx1
    nyk = nyest-ky1
    ncest = nxk*nyk
    nmx = nxest-nminx+1
    nmy = nyest-nminy+1
    nrint = nmx+nmy
    nreg = nmx*nmy
    ib1 = kx*nyk+ky1
    jb1 = ky*nxk+kx1
    ib3 = kx1*nyk+1
    if(ib1>jb1) then
      ib1 = jb1
      ib3 = ky1*nxk+1
    end if
    lwest = ncest*(2+ib1+ib3)+2*(nrint+nest*km2+m*km1)+ib3
    kwest = m+nreg
    if(lwrk1<lwest .or. kwrk<kwest) return
    if(xb>=xe .or. yb>=ye) return
    do i=1,m
      if(w(i)<=0.) return
      if(x(i)<xb .or. x(i)>xe) return
      if(y(i)<yb .or. y(i)>ye) return
    end do
    if(iopt<0) then
      if(nx<nminx .or. nx>nxest) return
      nxk = nx-kx1
      tx(kx1) = xb
      tx(nxk+1) = xe
      do i=kx1,nxk
        if(tx(i+1)<=tx(i)) return
      end do
      if(ny<nminy .or. ny>nyest) return
      nyk = ny-ky1
      ty(ky1) = yb
      ty(nyk+1) = ye
      do i=ky1,nyk
        if(ty(i+1)<=ty(i)) return
      end do
    else
      if(s<0.) return
    end if
    ier = 0
!  we partition the working space and determine the spline approximation
    kn = 1
    ki = kn+m
    lq = 2
    la = lq+ncest*ib3
    lf = la+ncest*ib1
    lff = lf+ncest
    lfp = lff+ncest
    lco = lfp+nrint
    lh = lco+nrint
    lbx = lh+ib3
    nek = nest*km2
    lby = lbx+nek
    lsx = lby+nek
    lsy = lsx+m*km1
    call fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
     eps,tol,maxit,nest,km1,km2,ib1,ib3,ncest,nrint,nreg,nx,tx, &
     ny,ty,c,fp,wrk1(1),wrk1(lfp),wrk1(lco),wrk1(lf),wrk1(lff), &
     wrk1(la),wrk1(lq),wrk1(lbx),wrk1(lby),wrk1(lsx),wrk1(lsy), &
     wrk1(lh),iwrk(ki),iwrk(kn),wrk2,lwrk2,ier)
  end subroutine surfit
!-------------------------------------------------------
  subroutine fpsurf(iopt,m,x,y,z,w,xb,xe,yb,ye,kxx,kyy,s,nxest, &
     nyest,eta,tol,maxit,nmax,km1,km2,ib1,ib3,nc,intest,nrest, &
     nx0,tx,ny0,ty,c,fp,fp0,fpint,coord,f,ff,a,q,bx,by,spx,spy,h, &
     idx,nummer,wrk,lwrk,ier)
!  ..
    use kinds, only : r8
    implicit none
!  ..scalar arguments..
    real(kind=r8) xb,xe,yb,ye,s,eta,tol,fp,fp0
    integer iopt,m,kxx,kyy,nxest,nyest,maxit,nmax,km1,km2,ib1,ib3, &
     nc,intest,nrest,nx0,ny0,lwrk,ier
!  ..array arguments..
    real(kind=r8) x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c(nc),fpint(intest), &
     coord(intest),f(nc),ff(nc),a(nc,ib1),q(nc,ib3),bx(nmax,km2), &
     by(nmax,km2),spx(m,km1),spy(m,km1),h(ib3),wrk(lwrk)
    integer idx(nrest),nummer(m)
!  ..local scalars..
    real(kind=r8) acc,arg,cs,dmax,fac1,fac2,fpmax,fpms,f1,f2,f3,hxi,p,pinv, &
     piv,p1,p2,p3,sigma,sn,sq,store,wi,x0,x1,y0,y1,zi,eps, &
     rn,one,con1,con9,con4,ten
    integer i,iband,iband1,iband3,iband4,ibb,ichang,ich1,ich3,ii, &
     in,irot,iter,i1,i2,i3,j,jrot,jxy,j1,kx,kx1,kx2,ky,ky1,ky2,l, &
     la,lf,lh,lwest,lx,ly,l1,l2,n,ncof,nk1x,nk1y,nminx,nminy,nreg, &
     nrint,num,num1,nx,nxe,nxx,ny,nye,nyy,n1,rank
!  ..local arrays..
    real(kind=r8) hx(6),hy(6)
!  ..function references..
!  real(8) abs,sqrt
!  integer min0
!  ..subroutine references..
!    fpback,fpbspl,fpgivs,fpdisc,fporde,fprank,fprota
  interface
    subroutine fpback(a,z,n,k,c,nest)
      use kinds, only : r8
      implicit none
      integer, intent(in)          :: n, k, nest
      real(kind=r8), intent(in)    :: a(nest,k), z(n)
      real(kind=r8), intent(inout) :: c(n)
    end subroutine fpback
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
    subroutine fpdisc(t,n,k2,b,nest)
      use kinds, only : r8
      implicit none
      integer      , intent(in)  :: n, k2, nest
      real(kind=r8), intent(in)  :: t(n)
      real(kind=r8), intent(out) :: b(nest,k2)
    end subroutine fpdisc
    subroutine fpgivs(piv,ww,cs,sn)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: piv
      real(kind=r8), intent(out)   :: cs, sn
      real(kind=r8), intent(inout) :: ww
    end subroutine fpgivs
    subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,idx,nreg)
      use kinds, only : r8
      implicit none
      integer m,kx,ky,nx,ny,nreg
      real(kind=r8) x(m),y(m),tx(nx),ty(ny)
      integer nummer(m),idx(nreg)
    end subroutine fporde
    subroutine fprota(cs,sn,a,b)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: cs, sn
      real(kind=r8), intent(inout) :: a, b
    end subroutine fprota
    subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
      use kinds, only : r8
      implicit none
      integer n,m,na,rank
      real(kind=r8) tol,sq
      real(kind=r8) a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
    end subroutine fprank
    function fprati(p1,f1,p2,f2,p3,f3)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fprati
      real(kind=r8), intent(in)    :: p2, f2
      real(kind=r8), intent(inout) :: p1, f1, p3, f3
    end function fprati
  end interface
!  ..
!  set constants
    one = 0.1e+01
    con1 = 0.1e0
    con9 = 0.9e0
    con4 = 0.4e-01
    ten = 0.1e+02
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 1: determination of the number of knots and their position.     c
! ****************************************************************     c
! given a set of knots we compute the least-squares spline sinf(x,y),  c
! and the corresponding weighted sum of squared residuals fp=f(p=inf). c
! if iopt=-1  sinf(x,y) is the requested approximation.                c
! if iopt=0 or iopt=1 we check whether we can accept the knots:        c
!   if fp <=s we will continue with the current set of knots.          c
!   if fp > s we will increase the number of knots and compute the     c
!  corresponding least-squares spline until finally  fp<=s.        c
! the initial choice of knots depends on the value of s and iopt.      c
!   if iopt=0 we first compute the least-squares polynomial of degree  c
! kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c
! fp0=f(0) denotes the corresponding weighted sum of squared       c
! residuals                                                        c
!   if iopt=1 we start with the knots found at the last call of the    c
! routine, except for the case that s>=fp0; then we can compute    c
! the least-squares polynomial directly.                           c
! eventually the independent variables x and y (and the corresponding  c
! parameters) will be switched if this can reduce the bandwidth of the c
! system to be solved.                                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ichang denotes whether(1) or not(-1) the directions have been inter-
!  changed.
    ichang = -1
    x0 = xb
    x1 = xe
    y0 = yb
    y1 = ye
    kx = kxx
    ky = kyy
    kx1 = kx+1
    ky1 = ky+1
    nxe = nxest
    nye = nyest
    eps = sqrt(eta)
    if(iopt<0) go to 20
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol*s
    if(iopt==0) go to 10
    if(fp0>s) go to 20
!  initialization for the least-squares polynomial.
    10 continue
    nminx = 2*kx1
    nminy = 2*ky1
    nx = nminx
    ny = nminy
    ier = -2
    go to 30
    20 continue
    nx = nx0
    ny = ny0
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
    30 continue
    do iter=1,m
!  find the position of the additional knots which are needed for the
!  b-spline representation of s(x,y).
      l = nx
      do i=1,kx1
        tx(i) = x0
        tx(l) = x1
        l = l-1
      end do
      l = ny
      do i=1,ky1
        ty(i) = y0
        ty(l) = y1
        l = l-1
      end do
!  find nrint, the total number of knot intervals and nreg, the number
!  of panels in which the approximation domain is subdivided by the
!  intersection of knots.
      nxx = nx-2*kx1+1
      nyy = ny-2*ky1+1
      nrint = nxx+nyy
      nreg = nxx*nyy
!  find the bandwidth of the observation matrix a.
!  if necessary, interchange the variables x and y, in order to obtain
!  a minimal bandwidth.
      iband1 = kx*(ny-ky1)+ky
      l = ky*(nx-kx1)+kx
      if(iband1>l) then
        iband1 = l
        ichang = -ichang
        do i=1,m
          store = x(i)
          x(i) = y(i)
          y(i) = store
        end do
        store = x0
        x0 = y0
        y0 = store
        store = x1
        x1 = y1
        y1 = store
        n = min0(nx,ny)
        do i=1,n
          store = tx(i)
          tx(i) = ty(i)
          ty(i) = store
        end do
        n1 = n+1
        if(nx<ny) then
          do i=n1,ny
            tx(i) = ty(i)
          end do
        else if (nx>ny) then
          do i=n1,nx
            ty(i) = tx(i)
          end do
        end if
        l = nx
        nx = ny
        ny = l
        l = nxe
        nxe = nye
        nye = l
        l = nxx
        nxx = nyy
        nyy = l
        l = kx
        kx = ky
        ky = l
        kx1 = kx+1
        ky1 = ky+1
      end if
      iband = iband1+1
!  arrange the data points according to the panel they belong to.
      call fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,idx,nreg)
!  find ncof, the number of b-spline coefficients.
      nk1x = nx-kx1
      nk1y = ny-ky1
      ncof = nk1x*nk1y
!  initialize the observation matrix a.
      do i=1,ncof
        f(i) = 0.
        do j=1,iband
          a(i,j) = 0.
        end do
      end do
!  initialize the sum of squared residuals.
      fp = 0.
!  fetch the data points in the new order. main loop for the
!  different panels.
      do num=1,nreg
!  fix certain constants for the current panel; jrot records the column
!  number of the first non-zero element in a row of the observation
!  matrix according to a data point of the panel.
        num1 = num-1
        lx = num1/nyy
        l1 = lx+kx1
        ly = num1-lx*nyy
        l2 = ly+ky1
        jrot = lx*nk1y+ly
!  test whether there are still data points in the panel.
        in = idx(num)
        do
          if(in==0) exit
!  fetch a new data point.
          wi = w(in)
          zi = z(in)*wi
!  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
          call fpbspl(tx,nx,kx,x(in),l1,hx)
!  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
          call fpbspl(ty,ny,ky,y(in),l2,hy)
!  store the value of these b-splines in spx and spy respectively.
          do i=1,kx1
            spx(in,i) = hx(i)
          end do
          do i=1,ky1
            spy(in,i) = hy(i)
          end do
!  initialize the new row of observation matrix.
          do i=1,iband
            h(i) = 0.
          end do
!  calculate the non-zero elements of the new row by making the cross
!  products of the non-zero b-splines in x- and y-direction.
          i1 = 0
          do i=1,kx1
            hxi = hx(i)
            j1 = i1
            do j=1,ky1
              j1 = j1+1
              h(j1) = hxi*hy(j)*wi
            end do
            i1 = i1+nk1y
          end do
!  rotate the row into triangle by givens transformations .
          irot = jrot
          do i=1,iband
            irot = irot+1
            piv = h(i)
            if(piv==0.) cycle
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(irot,1),cs,sn)
!  apply that transformation to the right hand side.
            call fprota(cs,sn,zi,f(irot))
            if(i==iband) exit
!  apply that transformation to the left hand side.
            i2 = 1
            i3 = i+1
            do j=i3,iband
              i2 = i2+1
              call fprota(cs,sn,h(j),a(irot,i2))
            end do
          end do
!  add the contribution of the row to the sum of squares of residual
!  right hand sides.
          fp = fp+zi**2
!  find the number of the next data point in the panel.
          in = nummer(in)
        end do
      end do
!  find dmax, the maximum value for the diagonal elements in the reduced
!  triangle.
      dmax = 0.
      do i=1,ncof
        if(a(i,1)<=dmax) cycle
        dmax = a(i,1)
      end do
!  check whether the observation matrix is rank deficient.
      sigma = eps*dmax
      do i=1,ncof
        if(a(i,1)<=sigma) go to 280
      end do
!  backward substitution in case of full rank.
      call fpback(a,f,ncof,iband,c,nc)
      rank = ncof
      do i=1,ncof
        q(i,1) = a(i,1)/dmax
      end do
      go to 300
!  in case of rank deficiency, find the minimum norm solution.
!  check whether there is sufficient working space
      280 continue
      lwest = ncof*iband+ncof+iband
      if(lwrk<lwest) go to 780
      do i=1,ncof
        ff(i) = f(i)
        do j=1,iband
          q(i,j) = a(i,j)
        end do
      end do
      lf =1
      lh = lf+ncof
      la = lh+iband
      call fprank(q,ff,ncof,iband,nc,sigma,c,sq,rank,wrk(la), &
        wrk(lf),wrk(lh))
      do i=1,ncof
        q(i,1) = q(i,1)/dmax
      end do
!  add to the sum of squared residuals, the contribution of reducing
!  the rank.
      fp = fp+sq
      300 continue
      if(ier==(-2)) fp0 = fp
!  test whether the least-squares spline is an acceptable solution.
      if(iopt<0) go to 820
      fpms = fp-s
      if(abs(fpms)<=acc) then
        if(fp>0.0d0) then
          go to 820
        else
          go to 815
        end if
      end if
!  test whether we can accept the choice of knots.
      if(fpms<0.) exit
!  test whether we cannot further increase the number of knots.
      if(ncof>m) go to 790
      ier = 0
!  search where to add a new knot.
!  find for each interval the sum of squared residuals fpint for the
!  data points having the coordinate belonging to that knot interval.
!  calculate also coord which is the same sum, weighted by the position
!  of the data points considered.
      do i=1,nrint
        fpint(i) = 0.
        coord(i) = 0.
      end do
      do num=1,nreg
        num1 = num-1
        lx = num1/nyy
        l1 = lx+1
        ly = num1-lx*nyy
        l2 = ly+1+nxx
        jrot = lx*nk1y+ly
        in = idx(num)
        do
          if(in==0) exit
          store = 0.
          i1 = jrot
          do i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do j=1,ky1
              j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
            end do
            i1 = i1+nk1y
          end do
          store = (w(in)*(z(in)-store))**2
          fpint(l1) = fpint(l1)+store
          coord(l1) = coord(l1)+store*x(in)
          fpint(l2) = fpint(l2)+store
          coord(l2) = coord(l2)+store*y(in)
          in = nummer(in)
        end do
      end do
!  find the interval for which fpint is maximal on the condition that
!  there still can be added a knot.
      do
        l = 0
        fpmax = 0.
        l1 = 1
        l2 = nrint
        if(nx==nxe) l1 = nxx+1
        if(ny==nye) l2 = nxx
        if(l1>l2) go to 810
        do i=l1,l2
          if(fpmax>=fpint(i)) cycle
          l = i
          fpmax = fpint(i)
        end do
!  test whether we cannot further increase the number of knots.
        if(l==0) go to 785
!  calculate the position of the new knot.
        arg = coord(l)/fpint(l)
!  test in what direction the new knot is going to be added.
        if(l<=nxx) then
!  addition in the x-direction.
          jxy = l+kx1
          fpint(l) = 0.
          fac1 = tx(jxy)-arg
          fac2 = arg-tx(jxy-1)
          if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) cycle
          j = nx
          do i=jxy,nx
            tx(j+1) = tx(j)
            j = j-1
          end do
          tx(jxy) = arg
          nx = nx+1
        else
!  addition in the y-direction.
          jxy = l+ky1-nxx
          fpint(l) = 0.
          fac1 = ty(jxy)-arg
          fac2 = arg-ty(jxy-1)
          if(fac1>(ten*fac2) .or. fac2>(ten*fac1)) cycle
          j = ny
          do i=jxy,ny
            ty(j+1) = ty(j)
            j = j-1
          end do
          ty(jxy) = arg
          ny = ny+1
        end if
        exit
      end do
!  restart the computations with the new set of knots.
    end do
!  test whether the least-squares polynomial is a solution of our
!  approximation problem.
    if(ier==(-2)) go to 830
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! part 2: determination of the smoothing spline sp(x,y)                c
! *****************************************************                c
! we have determined the number of knots and their position. we now    c
! compute the b-spline coefficients of the smoothing spline sp(x,y).   c
! the observation matrix a is extended by the rows of a matrix,        c
! expressing that sp(x,y) must be a polynomial of degree kx in x and   c
! ky in y. the corresponding weights of these additional rows are set  c
! to 1./p.  iteratively we than have to determine the value of p       c
! such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
! we already know that the least-squares polynomial corresponds to     c
! p=0  and that the least-squares spline corresponds to p=infinity.    c
! the iteration process which is proposed here makes use of rational   c
! interpolation. since f(p) is a convex and strictly decreasing        c
! function of p, it can be approximated by a rational function r(p)=   c
! (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
! of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
! new value of p such that r(p)=s. convergence is guaranteed by taking c
! f1 > 0 and f3 < 0.                                                   c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    kx2 = kx1+1
!  test whether there are interior knots in the x-direction.
!  and
!  evaluate the discotinuity jumps of the kx-th order derivative of
!  the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
    if(nk1x/=kx1) call fpdisc(tx,nx,kx2,bx,nmax) 
    ky2 = ky1 + 1
!  test whether there are interior knots in the y-direction.
!  and
!  evaluate the discontinuity jumps of the ky-th order derivative of
!  the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
    if(nk1y/=ky1) call fpdisc(ty,ny,ky2,by,nmax)
!  initial value for p.
    p1 = 0.
    f1 = fp0-s
    p3 = -one
    f3 = fpms
    p = 0.
    do i=1,ncof
      p = p+a(i,1)
    end do
    rn = ncof
    p = rn/p
!  find the bandwidth of the extended observation matrix.
    iband3 = kx1*nk1y
    iband4 = iband3 +1
    ich1 = 0
    ich3 = 0
!  iteration process to find the root of f(p)=s.
    do 770 iter=1,maxit
      pinv = one/p
!  store the triangularized observation matrix into q.
      do i=1,ncof
        ff(i) = f(i)
        do j=1,iband
          q(i,j) = a(i,j)
        end do
        ibb = iband+1
        do j=ibb,iband4
          q(i,j) = 0.
        end do
      end do
      if(nk1y/=ky1) then
!  extend the observation matrix with the rows of a matrix, expressing
!  that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
        do i=ky2,nk1y
          ii = i-ky1
          do j=1,nk1x
!  initialize the new row.
            do l=1,iband
              h(l) = 0.
            end do
!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.
            do l=1,ky2
              h(l) = by(ii,l)*pinv
            end do
            zi = 0.
            jrot = (j-1)*nk1y+ii
!  rotate the new row into triangle by givens transformations without
!  square roots.
            do irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband1,ncof-irot)
              if(piv/=0.) then
!  calculate the parameters of the givens transformation.
                call fpgivs(piv,q(irot,1),cs,sn)
!  apply that givens transformation to the right hand side.
                call fprota(cs,sn,zi,ff(irot))
                if(i2==0) exit
!  apply that givens transformation to the left hand side.
                do l=1,i2
                  l1 = l+1
                  call fprota(cs,sn,h(l1),q(irot,l1))
                end do
              else
                if(i2<=0) exit
              end if
              do l=1,i2
                h(l) = h(l+1)
              end do
              h(i2+1) = 0.
            end do
          end do
        end do
      end if
      if(nk1x/=kx1) then
!  extend the observation matrix with the rows of a matrix expressing
!  that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
        do i=kx2,nk1x
          ii = i-kx1
          do j=1,nk1y
!  initialize the new row
            do l=1,iband4
              h(l) = 0.
            end do
!  fill in the non-zero elements of the row. jrot records the column
!  number of the first non-zero element in the row.
            j1 = 1
            do l=1,kx2
              h(j1) = bx(ii,l)*pinv
              j1 = j1+nk1y
            end do
            zi = 0.
            jrot = (i-kx2)*nk1y+j
!  rotate the new row into triangle by givens transformations .
            do irot=jrot,ncof
              piv = h(1)
              i2 = min0(iband3,ncof-irot)
              if(piv/=0.) then
!  calculate the parameters of the givens transformation.
                call fpgivs(piv,q(irot,1),cs,sn)
!  apply that givens transformation to the right hand side.
                call fprota(cs,sn,zi,ff(irot))
                if(i2==0) exit
!  apply that givens transformation to the left hand side.
                do l=1,i2
                  l1 = l+1
                  call fprota(cs,sn,h(l1),q(irot,l1))
                end do
              else
                if(i2<=0) exit
              end if
              do l=1,i2
                h(l) = h(l+1)
              end do
              h(i2+1) = 0.
            end do
          end do
        end do
      end if
!  find dmax, the maximum value for the diagonal elements in the
!  reduced triangle.
      dmax = 0.
      do i=1,ncof
        if(q(i,1)<=dmax) cycle
        dmax = q(i,1)
      end do
!  check whether the matrix is rank deficient.
      sigma = eps*dmax
      do i=1,ncof
        if(q(i,1)<=sigma) go to 670
      end do
!  backward substitution in case of full rank.
      call fpback(q,ff,ncof,iband4,c,nc)
      rank = ncof
      go to 675
!  in case of rank deficiency, find the minimum norm solution.
      670 continue
      lwest = ncof*iband4+ncof+iband4
      if(lwrk<lwest) go to 780
      lf = 1
      lh = lf+ncof
      la = lh+iband4
      call fprank(q,ff,ncof,iband4,nc,sigma,c,sq,rank,wrk(la), &
       wrk(lf),wrk(lh))
      675 continue
      do i=1,ncof
        q(i,1) = q(i,1)/dmax
      end do
!  compute f(p).
      fp = 0.
      do num = 1,nreg
        num1 = num-1
        lx = num1/nyy
        ly = num1-lx*nyy
        jrot = lx*nk1y+ly
        in = idx(num)
        do
          if(in==0) exit
          store = 0.
          i1 = jrot
          do i=1,kx1
            hxi = spx(in,i)
            j1 = i1
            do j=1,ky1
             j1 = j1+1
              store = store+hxi*spy(in,j)*c(j1)
            end do
            i1 = i1+nk1y
          end do
          fp = fp+(w(in)*(z(in)-store))**2
          in = nummer(in)
        end do
      end do
!  test whether the approximation sp(x,y) is an acceptable solution.
      fpms = fp-s
      if(abs(fpms)<=acc) go to 820
!  test whether the maximum allowable number of iterations has been
!  reached.
      if(iter==maxit) go to 795
!  carry out one more step of the iteration process.
      p2 = p
      f2 = fpms
      if(ich3/=0) go to 740
      if((f2-f3)>acc) go to 730
!  our initial choice of p is too large.
      p3 = p2
      f3 = f2
      p = p*con4
      if(p<=p1) p = p1*con9 + p2*con1
      go to 770
   730    if(f2<0.) ich3 = 1
   740    if(ich1/=0) go to 760
      if((f1-f2)>acc) go to 750
!  our initial choice of p is too small
      p1 = p2
      f1 = f2
      p = p/con4
      if(p3<0.) go to 770
      if(p>=p3) p = p2*con1 + p3*con9
      go to 770
   750    if(f2>0.) ich1 = 1
!  test whether the iteration process proceeds as theoretically
!  expected.
   760    if(f2>=f1 .or. f2<=f3) go to 800
!  find the new value of p.
      p = fprati(p1,f1,p2,f2,p3,f3)
   770  continue
!  error codes and messages.
   780  ier = lwest
    go to 830
   785  ier = 5
    go to 830
   790  ier = 4
    go to 830
   795  ier = 3
    go to 830
   800  ier = 2
    go to 830
   810  ier = 1
    go to 830
   815  ier = -1
    fp = 0.
   820  if(ncof/=rank) ier = -rank
!  test whether x and y are in the original order.
   830  if(ichang<0) go to 930
!  if not, interchange x and y once more.
    l1 = 1
    do i=1,nk1x
      l2 = i
      do j=1,nk1y
        f(l2) = c(l1)
        l1 = l1+1
        l2 = l2+nk1x
      end do
    end do
    do i=1,ncof
      c(i) = f(i)
    end do
    do i=1,m
      store = x(i)
      x(i) = y(i)
      y(i) = store
    end do
    n = min0(nx,ny)
    do i=1,n
      store = tx(i)
      tx(i) = ty(i)
      ty(i) = store
    end do
    n1 = n+1
    if(nx<ny) then
      do i=n1,ny
        tx(i) = ty(i)
      end do
    else if(nx>ny) then
      do i=n1,nx
        ty(i) = tx(i)
      end do
    end if
    l = nx
    nx = ny
    ny = l
    930 continue
    if(iopt>=0) then
      nx0 = nx
      ny0 = ny
    end if
  end subroutine fpsurf
!-------------------------------------------------------
  subroutine fprank(a,f,n,m,na,tol,c,sq,rank,aa,ff,h)
!  subroutine fprank finds the minimum norm solution of a least-
!  squares problem in case of rank deficiency.
!
!  input parameters:
!    a : array, which contains the non-zero elements of the observation
!    matrix after triangularization by givens transformations.
!    f : array, which contains the transformed right hand side.
!    n : integer,wich contains the dimension of a.
!    m : integer, which denotes the bandwidth of a.
!  tol : real value, giving a threshold to determine the rank of a.
!
!  output parameters:
!    c : array, which contains the minimum norm solution.
!   sq : real value, giving the contribution of reducing the rank
!    to the sum of squared residuals.
! rank : integer, which contains the rank of matrix a.
!
    use kinds, only : r8
    implicit none
!  ..scalar arguments..
    integer n,m,na,rank
    real(kind=r8) tol,sq
!  ..array arguments..
    real(kind=r8) a(na,m),f(n),c(n),aa(n,m),ff(n),h(m)
!  ..local scalars..
    integer i,ii,ij,i1,i2,j,jj,j1,j2,j3,k,kk,m1,nl
    real(kind=r8) cs,fac,piv,sn,yi
    real(kind=r8) store,stor1,stor2,stor3
!  ..function references..
    integer min0
!  ..subroutine references..
!    fpgivs,fprota
  interface
    subroutine fpgivs(piv,ww,cs,sn)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: piv
      real(kind=r8), intent(out)   :: cs, sn
      real(kind=r8), intent(inout) :: ww
    end subroutine fpgivs
    subroutine fprota(cs,sn,a,b)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: cs, sn
      real(kind=r8), intent(inout) :: a, b
    end subroutine fprota
  end interface
!  ..
    m1 = m-1
!  the rank deficiency nl is considered to be the number of sufficient
!  small diagonal elements of a.
    nl = 0
    sq = 0.
    do i=1,n
      if(a(i,1)>tol) cycle
!  if a sufficient small diagonal element is found, we put it to
!  zero. the remainder of the row corresponding to that zero diagonal
!  element is then rotated into triangle by givens rotations .
!  the rank deficiency is increased by one.
      nl = nl+1
      if(i==n) cycle
      yi = f(i)
      do j=1,m1
        h(j) = a(i,j+1)
      end do
      h(m) = 0.
      i1 = i+1
      do ii=i1,n
        i2 = min0(n-ii,m1)
        piv = h(1)
        if(piv/=0.) then
          call fpgivs(piv,a(ii,1),cs,sn)
          call fprota(cs,sn,yi,f(ii))
          if(i2==0) exit
          do j=1,i2
            j1 = j+1
            call fprota(cs,sn,h(j1),a(ii,j1))
            h(j) = h(j1)
          end do
        else
          if(i2==0) exit
          do j=1,i2
            h(j) = h(j+1)
          end do
        end if
        h(i2+1) = 0.
      end do
!  add to the sum of squared residuals the contribution of deleting
!  the row with small diagonal element.
      sq = sq+yi**2
    end do
!  rank denotes the rank of a.
    rank = n-nl
!  let b denote the (rank*n) upper trapezoidal matrix which can be
!  obtained from the (n*n) upper triangular matrix a by deleting
!  the rows and interchanging the columns corresponding to a zero
!  diagonal element. if this matrix is factorized using givens
!  transformations as  b = (r) (u)  where
!    r is a (rank*rank) upper triangular matrix,
!    u is a (rank*n) orthonormal matrix
!  then the minimal least-squares solution c is given by c = b' v,
!  where v is the solution of the system  (r) (r)' v = g  and
!  g denotes the vector obtained from the old right hand side f, by
!  removing the elements corresponding to a zero diagonal element of a.
!  initialization.
    do i=1,rank
      do j=1,m
        aa(i,j) = 0.
      end do
    end do
!  form in aa the upper triangular matrix obtained from a by
!  removing rows and columns with zero diagonal elements. form in ff
!  the new right hand side by removing the elements of the old right
!  hand side corresponding to a deleted row.
    ii = 0
    do i=1,n
      if(a(i,1)<=tol) cycle
      ii = ii+1
      ff(ii) = f(i)
      aa(ii,1) = a(i,1)
      jj = ii
      kk = 1
      j = i
      j1 = min0(j-1,m1)
      do k=1,j1
        j = j-1
        if(a(j,1)<=tol) cycle
        kk = kk+1
        jj = jj-1
        aa(jj,kk) = a(j,k+1)
      end do
    end do
!  form successively in h the columns of a with a zero diagonal element.
    ii = 0
    do i=1,n
      ii = ii+1
      if(a(i,1)>tol) cycle
      ii = ii-1
      if(ii==0) cycle
      jj = 1
      j = i
      j1 = min0(j-1,m1)
      do k=1,j1
        j = j-1
        if(a(j,1)<=tol) cycle
        h(jj) = a(j,k+1)
        jj = jj+1
      end do
      do kk=jj,m
        h(kk) = 0.
      end do
!  rotate this column into aa by givens transformations.
      jj = ii
      do i1=1,ii
        j1 = min0(jj-1,m1)
        piv = h(1)
        if(piv==0.) then
          if(j1==0) exit
          do j2=1,j1
            j3 = j2+1
            h(j2) = h(j3)
          end do
        else
          call fpgivs(piv,aa(jj,1),cs,sn)
          if(j1==0) exit
          kk = jj
          do j2=1,j1
            j3 = j2+1
            kk = kk-1
            call fprota(cs,sn,h(j3),aa(kk,j3))
            h(j2) = h(j3)
          end do
        end if
        jj = jj-1
        h(j3) = 0.
      end do
    end do
!  solve the system (aa) (f1) = ff
    ff(rank) = ff(rank)/aa(rank,1)
    i = rank-1
    do j=2,rank
      store = ff(i)
      i1 = min0(j-1,m1)
      k = i
      do ii=1,i1
        k = k+1
        stor1 = ff(k)
        stor2 = aa(i,ii+1)
        store = store-stor1*stor2
      end do
      stor1 = aa(i,1)
      ff(i) = store/stor1
      i = i-1
    end do
!  solve the system  (aa)' (f2) = f1
    ff(1) = ff(1)/aa(1,1)
    do j=2,rank
      store = ff(j)
      i1 = min0(j-1,m1)
      k = j
      do ii=1,i1
        k = k-1
        stor1 = ff(k)
        stor2 = aa(k,ii+1)
        store = store-stor1*stor2
      end do
      stor1 = aa(j,1)
      ff(j) = store/stor1
    end do
!  premultiply f2 by the transpoze of a.
    k = 0
    do i=1,n
      store = 0.
      if(a(i,1)>tol) k = k+1
      j1 = min0(i,m)
      kk = k
      ij = i+1
      do j=1,j1
        ij = ij-1
        if(a(ij,1)<=tol) cycle
        stor1 = a(ij,j)
        stor2 = ff(kk)
        store = store+stor1*stor2
        kk = kk-1
      end do
      c(i) = store
    end do
!  add to the sum of squared residuals the contribution of putting
!  to zero the small diagonal elements of matrix (a).
    stor3 = 0.
    do i=1,n
      if(a(i,1)>tol) cycle
      store = f(i)
      i1 = min0(n-i,m1)
      do j=1,i1
        ij = i+j
        stor1 = c(ij)
        stor2 = a(i,j+1)
        store = store-stor1*stor2
      end do
      fac = a(i,1)*c(i)
      stor1 = a(i,1)
      stor2 = c(i)
      stor1 = stor1*stor2
      stor3 = stor3+stor1*(stor1-store-store)
    end do
    fac = stor3
    sq = sq+fac
  end subroutine fprank
!-------------------------------------------------------
  subroutine fporde(x,y,m,kx,ky,tx,nx,ty,ny,nummer,idx,nreg)
!  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
!  according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
!  to. for each panel a stack is constructed  containing the numbers
!  of data points lying inside; index(j),j=1,2,...,nreg points to the
!  first data point in the jth panel while nummer(i),i=1,2,...,m gives
!  the number of the next data point in the panel.
!  ..
    use kinds, only : r8
    implicit none
!  ..scalar arguments..
    integer m,kx,ky,nx,ny,nreg
!  ..array arguments..
    real(kind=r8) x(m),y(m),tx(nx),ty(ny)
    integer nummer(m),idx(nreg)
!  ..local scalars..
    real(kind=r8) xi,yi
    integer i,im,k,kx1,ky1,k1,l,l1,nk1x,nk1y,num,nyy
!  ..
    kx1 = kx+1
    ky1 = ky+1
    nk1x = nx-kx1
    nk1y = ny-ky1
    nyy = nk1y-ky
    do i=1,nreg
      idx(i) = 0
    end do
    do im=1,m
      xi = x(im)
      yi = y(im)
      l = kx1
      l1 = l+1
      do
        if(xi<tx(l1) .or. l==nk1x) exit
        l = l1
        l1 = l+1
      end do
      k = ky1
      k1 = k+1
      do
        if(yi<ty(k1) .or. k==nk1y) exit
        k = k1
        k1 = k+1
      end do
      num = (l-kx1)*nyy+k-ky
      nummer(im) = idx(num)
      idx(num) = im
    end do
  end subroutine fporde
!-------------------------------------------------------
subroutine regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,&
   nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
! given the set of values z(i,j) on the rectangular grid (x(i),y(j)),
! i=1,...,mx;j=1,...,my, subroutine regrid determines a smooth bivar-
! iate spline approximation s(x,y) of degrees kx and ky on the rect-
! angle xb <= x <= xe, yb <= y <= ye.
! if iopt = -1 regrid calculates the least-squares spline according
! to a given set of knots.
! if iopt >= 0 the total numbers nx and ny of these knots and their
! position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
! ally by the routine. the smoothness of s(x,y) is then achieved by
! minimalizing the discontinuity jumps in the derivatives of s(x,y)
! across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
! the amounth of smoothness is determined by the condition that f(p) =
! sum ((z(i,j)-s(x(i),y(j))))**2) be <= s, with s a given non-negative
! constant, called the smoothing factor.
! the fit is given in the b-spline representation (b-spline coefficients
! c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
! uated by means of subroutine bispev.
!
! calling sequence:
! call regrid(iopt,mx,x,my,y,z,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
!*  nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
!
! parameters:
!  iopt  : integer flag. on entry iopt must specify whether a least-
!      squares spline (iopt=-1) or a smoothing spline (iopt=0 or 1)
!      must be determined.
!      if iopt=0 the routine will start with an initial set of knots
!      tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
!      1,...,ky+1. if iopt=1 the routine will continue with the set
!      of knots found at the last call of the routine.
!      attention: a call with iopt=1 must always be immediately pre-
!                 ceded by another call with iopt=1 or iopt=0 and
!                 s/=0.
!      unchanged on exit.
!  mx    : integer. on entry mx must specify the number of grid points
!      along the x-axis. mx > kx . unchanged on exit.
!  x     : real array of dimension at least (mx). before entry, x(i)
!      must be set to the x-co-ordinate of the i-th grid point
!      along the x-axis, for i=1,2,...,mx. these values must be
!      supplied in strictly ascending order. unchanged on exit.
!  my    : integer. on entry my must specify the number of grid points
!      along the y-axis. my > ky . unchanged on exit.
!  y     : real array of dimension at least (my). before entry, y(j)
!      must be set to the y-co-ordinate of the j-th grid point
!      along the y-axis, for j=1,2,...,my. these values must be
!      supplied in strictly ascending order. unchanged on exit.
!  z     : real array of dimension at least (mx*my).
!      before entry, z(my*(i-1)+j) must be set to the data value at
!      the grid point (x(i),y(j)) for i=1,...,mx and j=1,...,my.
!      unchanged on exit.
!  xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
!  yb,ye   aries of the rectangular approximation domain.
!      xb<=x(i)<=xe,i=1,...,mx; yb<=y(j)<=ye,j=1,...,my.
!      unchanged on exit.
!  kx,ky : integer values. on entry kx and ky must specify the degrees
!      of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
!      (kx=ky=3) splines. unchanged on exit.
!  s     : real. on entry (in case iopt>=0) s must specify the smoothing
!      factor. s >=0. unchanged on exit.
!      for advice on the choice of s see further comments
!  nxest : integer. unchanged on exit.
!  nyest : integer. unchanged on exit.
!      on entry, nxest and nyest must specify an upper bound for the
!      number of knots required in the x- and y-directions respect.
!      these numbers will also determine the storage space needed by
!      the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
!      in most practical situation nxest = mx/2, nyest=my/2, will
!      be sufficient. always large enough are nxest=mx+kx+1, nyest=
!      my+ky+1, the number of knots needed for interpolation (s=0).
!      see also further comments.
!  nx    : integer.
!      unless ier=10 (in case iopt >=0), nx will contain the total
!      number of knots with respect to the x-variable, of the spline
!      approximation returned. if the computation mode iopt=1 is
!      used, the value of nx should be left unchanged between sub-
!      sequent calls.
!      in case iopt=-1, the value of nx should be specified on entry
!  tx    : real array of dimension nmax.
!      on succesful exit, this array will contain the knots of the
!      spline with respect to the x-variable, i.e. the position of
!      the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
!      position of the additional knots tx(1)=...=tx(kx+1)=xb and
!      tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
!      if the computation mode iopt=1 is used, the values of tx(1),
!      ...,tx(nx) should be left unchanged between subsequent calls.
!      if the computation mode iopt=-1 is used, the values tx(kx+2),
!      ...tx(nx-kx-1) must be supplied by the user, before entry.
!      see also the restrictions (ier=10).
!  ny    : integer.
!      unless ier=10 (in case iopt >=0), ny will contain the total
!      number of knots with respect to the y-variable, of the spline
!      approximation returned. if the computation mode iopt=1 is
!      used, the value of ny should be left unchanged between sub-
!      sequent calls.
!      in case iopt=-1, the value of ny should be specified on entry
!  ty    : real array of dimension nmax.
!      on succesful exit, this array will contain the knots of the
!      spline with respect to the y-variable, i.e. the position of
!      the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
!      position of the additional knots ty(1)=...=ty(ky+1)=yb and
!      ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
!      if the computation mode iopt=1 is used, the values of ty(1),
!      ...,ty(ny) should be left unchanged between subsequent calls.
!      if the computation mode iopt=-1 is used, the values ty(ky+2),
!      ...ty(ny-ky-1) must be supplied by the user, before entry.
!      see also the restrictions (ier=10).
!  c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
!      on succesful exit, c contains the coefficients of the spline
!      approximation s(x,y)
!  fp    : real. unless ier=10, fp contains the sum of squared
!      residuals of the spline approximation returned.
!  wrk   : real array of dimension (lwrk). used as workspace.
!      if the computation mode iopt=1 is used the values of wrk(1),
!      ...,wrk(4) should be left unchanged between subsequent calls.
!  lwrk  : integer. on entry lwrk must specify the actual dimension of
!      the array wrk as declared in the calling (sub)program.
!      lwrk must not be too small.
!       lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
!        my*(ky+1) +u
!       where u is the larger of my and nxest.
!  iwrk  : integer array of dimension (kwrk). used as workspace.
!      if the computation mode iopt=1 is used the values of iwrk(1),
!      ...,iwrk(3) should be left unchanged between subsequent calls
!  kwrk  : integer. on entry kwrk must specify the actual dimension of
!      the array iwrk as declared in the calling (sub)program.
!      kwrk >= 3+mx+my+nxest+nyest.
!  ier   : integer. unless the routine detects an error, ier contains a
!      non-positive value on exit, i.e.
!   ier=0  : normal return. the spline returned has a residual sum of
!        squares fp such that abs(fp-s)/s <= tol with tol a relat-
!        ive tolerance set to 0.001 by the program.
!   ier=-1 : normal return. the spline returned is an interpolating
!        spline (fp=0).
!   ier=-2 : normal return. the spline returned is the least-squares
!        polynomial of degrees kx and ky. in this extreme case fp
!        gives the upper bound for the smoothing factor s.
!   ier=1  : error. the required storage space exceeds the available
!        storage space, as specified by the parameters nxest and
!        nyest.
!        probably causes : nxest or nyest too small. if these param-
!        eters are already large, it may also indicate that s is
!        too small
!        the approximation returned is the least-squares spline
!        according to the current set of knots. the parameter fp
!        gives the corresponding sum of squared residuals (fp>s).
!   ier=2  : error. a theoretically impossible result was found during
!        the iteration proces for finding a smoothing spline with
!        fp = s. probably causes : s too small.
!        there is an approximation returned but the corresponding
!        sum of squared residuals does not satisfy the condition
!        abs(fp-s)/s < tol.
!   ier=3  : error. the maximal number of iterations maxit (set to 20
!        by the program) allowed for finding a smoothing spline
!        with fp=s has been reached. probably causes : s too small
!        there is an approximation returned but the corresponding
!        sum of squared residuals does not satisfy the condition
!        abs(fp-s)/s < tol.
!   ier=10 : error. on entry, the input data are controlled on validity
!        the following restrictions must be satisfied.
!        -1<=iopt<=1, 1<=kx,ky<=5, mx>kx, my>ky, nxest>=2*kx+2,
!        nyest>=2*ky+2, kwrk>=3+mx+my+nxest+nyest,
!        lwrk >= 4+nxest*(my+2*kx+5)+nyest*(2*ky+5)+mx*(kx+1)+
!         my*(ky+1) +max(my,nxest),
!        xb<=x(i-1)<x(i)<=xe,i=2,..,mx,yb<=y(j-1)<y(j)<=ye,j=2,..,my
!        if iopt=-1: 2*kx+2<=nx<=min(nxest,mx+kx+1)
!                    xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
!                    2*ky+2<=ny<=min(nyest,my+ky+1)
!                    yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
!                the schoenberg-whitney conditions, i.e. there must
!                be subset of grid co-ordinates xx(p) and yy(q) such
!                that   tx(p) < xx(p) < tx(p+kx+1) ,p=1,...,nx-kx-1
!                       ty(q) < yy(q) < ty(q+ky+1) ,q=1,...,ny-ky-1
!        if iopt>=0: s>=0
!                    if s=0 : nxest>=mx+kx+1, nyest>=my+ky+1
!        if one of these conditions is found to be violated,control
!        is immediately repassed to the calling program. in that
!        case there is no approximation returned.
!
! further comments:
!   regrid does not allow individual weighting of the data-values.
!   so, if these were determined to widely different accuracies, then
!   perhaps the general data set routine surfit should rather be used
!   in spite of efficiency.
!   by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating spline if
!   s=0 and the least-squares polynomial (degrees kx,ky) if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the accuracy of the data values.
!   if the user has an idea of the statistical errors on the data, he
!   can also find a proper estimate for s. for, by assuming that, if he
!   specifies the right s, regrid will return a spline s(x,y) which
!   exactly reproduces the function underlying the data he can evaluate
!   the sum((z(i,j)-s(x(i),y(j)))**2) to find a good estimate for this s
!   for example, if he knows that the statistical errors on his z(i,j)-
!   values is not greater than 0.1, he may expect that a good s should
!   have a value not larger than mx*my*(0.1)**2.
!   if nothing is known about the statistical error in z(i,j), s must
!   be determined by trial and error, taking account of the comments
!   above. the best is then to start with a very large value of s (to
!   determine the least-squares polynomial and the corresponding upper
!   bound fp0 for s) and then to progressively decrease the value of s
!   ( say by a factor 10 in the beginning, i.e. s=fp0/10,fp0/100,...
!   and more carefully as the approximation shows more detail) to
!   obtain closer fits.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if regrid is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. if the computation mode iopt=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   regrid once more with the selected value for s but now with iopt=0.
!   indeed, regrid may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!   the number of knots may also depend on the upper bounds nxest and
!   nyest. indeed, if at a certain stage in regrid the number of knots
!   in one direction (say nx) has reached the value of its upper bound
!   (nxest), then from that moment on all subsequent knots are added
!   in the other (y) direction. this may indicate that the value of
!   nxest is too small. on the other hand, it gives the user the option
!   of limiting the number of knots the routine locates in any direction
!   for example, by setting nxest=2*kx+2 (the lowest allowable value for
!   nxest), the user can indicate that he wants an approximation which
!   is a simple polynomial of degree kx in the variable x.
!
!  other subroutines required:
!fpback,fpbspl,fpregr,fpdisc,fpgivs,fpgrre,fprati,fprota,fpchec,
!fpknot
!
!  references:
!   dierckx p. : a fast algorithm for smoothing data on a rectangular
!            grid while using spline functions, siam j.numer.anal.
!            19 (1982) 1286-1304.
!   dierckx p. : a fast algorithm for smoothing data on a rectangular
!            grid while using spline functions, report tw53, dept.
!            computer science,k.u.leuven, 1980.
!   dierckx p. : curve and surface fitting with splines, monographs on
!            numerical analysis, oxford university press, 1993.
!
!  author:
!p.dierckx
!dept. computer science, k.u. leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : may 1979
!  latest update : march 1989
!
!  ..
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: iopt, mx, my, kx, ky, nxest, nyest, lwrk, kwrk
  integer, intent(out)   :: ier
  integer, intent(inout) :: nx, ny, iwrk(kwrk)
  real(kind=r8), intent(in)    :: xb, xe, yb, ye, x(mx), y(my), z(mx*my), s
  real(kind=r8), intent(out)   :: fp, c((nxest-kx-1)*(nyest-ky-1))
  real(kind=r8), intent(inout) :: tx(nxest), ty(nyest), wrk(lwrk)
! local variables
  integer :: i, j, jwrk, kndx, kndy, knrx, knry, kwest, kx1, kx2, &
             ky1, ky2, lfpx, lfpy, lwest, lww, nc, nminx, nminy, mz
! parameters
  integer      , parameter :: maxit = 20
  real(kind=r8), parameter :: tol = 0.1d-02
! intrinsic functions
!  integer, intrinsic :: max
! external subroutines
  interface
    subroutine fpchec(x,m,t,n,k,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)       :: m, n, k
      real(kind=r8), intent(in) :: x(m), t(n)
      integer, intent(out)      :: ier
    end subroutine fpchec
    subroutine fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,tol,&
       maxit,nc,nx,tx,ny,ty,c,fp,fp0,fpold,reducx,reducy,fpintx,fpinty,lastdi,&
       nplusx,nplusy,nrx,nry,nrdatx,nrdaty,wrk,lwrk,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: kx,ky,mx,my,mz,nxest,nyest,maxit,nc,lwrk,iopt
      integer, intent(inout) :: nx,ny,lastdi,nplusx,nplusy,                   &
                                nrdatx(nxest),nrdaty(nyest),nrx(mx),nry(my)
      integer, intent(out)   :: ier
      real(kind=r8), intent(in)    :: xb,xe,yb,ye,s,tol,x(mx),y(my),z(mz)
      real(kind=r8), intent(out)   :: fp,c(nc)
      real(kind=r8), intent(inout) :: fp0,fpold,reducx,reducy,tx(nxest),      &
                        ty(nyest),fpintx(nxest),fpinty(nyest),wrk(lwrk)
    end subroutine fpregr
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
  ier = 10
  if (kx<=0 .or. kx>5) return
  kx1 = kx+1
  kx2 = kx1+1
  if (ky<=0 .or. ky>5) return
  ky1 = ky+1
  ky2 = ky1+1
  if (iopt<(-1) .or. iopt>1) return
  nminx = 2*kx1
  if (mx<kx1 .or. nxest<nminx) return
  nminy = 2*ky1
  if (my<ky1 .or. nyest<nminy) return
  mz = mx*my
  nc = (nxest-kx1)*(nyest-ky1)
  lwest = 4+nxest*(my+2*kx2+1)+nyest*(2*ky2+1)+mx*kx1+&
   my*ky1+max(nxest,my)
  kwest = 3+mx+my+nxest+nyest
  if (lwrk<lwest .or. kwrk<kwest) return
  if (xb>x(1) .or. xe<x(mx)) return
  do i=2,mx
    if (x(i-1)>=x(i)) return
  end do
  if (yb>y(1) .or. ye<y(my)) return
  do i=2,my
    if (y(i-1)>=y(i)) return
  end do
  if (iopt>=0) then
    if (s<0.0d0) return
    if (s==0.0d0 .and. (nxest<(mx+kx1).or.nyest<(my+ky1)))&
   return
    ier = 0
  else
    if (nx<nminx .or. nx>nxest) return
    j = nx
    do i=1,kx1
      tx(i) = xb
      tx(j) = xe
      j = j-1
    end do
    call fpchec(x,mx,tx,nx,kx,ier)
    if (ier/=0) return
    if (ny<nminy .or. ny>nyest) return
    j = ny
    do i=1,ky1
      ty(i) = yb
      ty(j) = ye
      j = j-1
    end do
    call fpchec(y,my,ty,ny,ky,ier)
    if (ier/=0) return
  end if
!  we partition the working space and determine the spline approximation
  lfpx = 5
  lfpy = lfpx+nxest
  lww = lfpy+nyest
  jwrk = lwrk-4-nxest-nyest
  knrx = 4
  knry = knrx+mx
  kndx = knry+my
  kndy = kndx+nxest
  call fpregr(iopt,x,mx,y,my,z,mz,xb,xe,yb,ye,kx,ky,s,nxest,nyest,            &
   tol,maxit,nc,nx,tx,ny,ty,c,fp,wrk(1),wrk(2),wrk(3),wrk(4),                 &
   wrk(lfpx:lfpy-1),wrk(lfpy:lww-1),iwrk(1),iwrk(2),iwrk(3),iwrk(knrx:knry-1),&
   iwrk(knry:kndx-1),iwrk(kndx:kndy-1),iwrk(kndy:),wrk(lww:),jwrk,ier)
end subroutine regrid
!-------------------------------------------------------
subroutine parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,&
   wrk,lwrk,iwrk,kwrk,ier)
!  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
!  ,my the partial derivative ( order nux,nuy) of a bivariate spline
!  s(x,y) of degrees kx and ky, given in the b-spline representation.
!
!  calling sequence:
! call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk,
!* iwrk,kwrk,ier)
!
!  input parameters:
!   tx    : real array, length nx, which contains the position of the
!       knots in the x-direction.
!   nx    : integer, giving the total number of knots in the x-direction
!   ty    : real array, length ny, which contains the position of the
!       knots in the y-direction.
!   ny    : integer, giving the total number of knots in the y-direction
!   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!       b-spline coefficients.
!   kx,ky : integer values, giving the degrees of the spline.
!   nux   : integer values, specifying the order of the partial
!   nuy     derivative. 0<=nux<kx, 0<=nuy<ky.
!   x     : real array of dimension (mx).
!       before entry x(i) must be set to the x co-ordinate of the
!       i-th grid point along the x-axis.
!       tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
!   mx    : on entry mx must specify the number of grid points along
!       the x-axis. mx >=1.
!   y     : real array of dimension (my).
!       before entry y(j) must be set to the y co-ordinate of the
!       j-th grid point along the y-axis.
!       ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
!   my    : on entry my must specify the number of grid points along
!       the y-axis. my >=1.
!   wrk   : real array of dimension lwrk. used as workspace.
!   lwrk  : integer, specifying the dimension of wrk.
!       lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1)
!   iwrk  : integer array of dimension kwrk. used as workspace.
!   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.
!
!  output parameters:
!   z     : real array of dimension (mx*my).
!       on succesful exit z(my*(i-1)+j) contains the value of the
!       specified partial derivative of s(x,y) at the point
!       (x(i),y(j)),i=1,...,mx;j=1,...,my.
!   ier   : integer error flag
!ier=0 : normal return
!ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my
!   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1),
!   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
!   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my
!
!  other subroutines required:
!fpbisp,fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!            6 (1972) 50-62.
!   dierckx p. : curve and surface fitting with splines, monographs on
!            numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1989
!
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: nx, ny, kx, ky, nux, nuy, mx, my, lwrk, kwrk
  integer, intent(out)   :: ier
  integer, intent(inout) :: iwrk(kwrk)
  real(kind=r8), intent(in)    :: tx(nx), ty(ny), c((nx-kx-1)*(ny-ky-1)), &
                                  x(mx), y(my)
  real(kind=r8), intent(out)   :: z(mx*my)
  real(kind=r8), intent(inout) :: wrk(lwrk)
! local variables
  integer       :: i, iwx, iwy, j, kkx, kky, kx1, ky1, lx, ly, lwest,   &
                   l1, l2, m, m0, m1, nc, nkx1, nky1, nxx, nyy
  real(kind=r8) :: ak, fac
! external subroutines
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
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
  ier = 10
  kx1 = kx+1
  ky1 = ky+1
  nkx1 = nx-kx1
  nky1 = ny-ky1
  nc = nkx1*nky1
  if (nux<0 .or. nux>=kx) return
  if (nuy<0 .or. nuy>=ky) return
  lwest = nc +(kx1-nux)*mx+(ky1-nuy)*my
  if (lwrk<lwest) return
  if (kwrk<(mx+my)) return
  if (mx<1) return
  do i=2,mx
    if(x(i)<x(i-1)) return
  end do
  if (my<1) return
  do i=2,my
    if(y(i)<y(i-1)) return
  end do
  ier = 0
  nxx = nkx1
  nyy = nky1
  kkx = kx
  kky = ky
!  the partial derivative of order (nux,nuy) of a bivariate spline of
!  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
!  we calculate the b-spline coefficients of this spline
  do i=1,nc
    wrk(i) = c(i)
  end do
  if(nux/=0) then
    lx = 1
    do j=1,nux
      ak = kkx
      nxx = nxx-1
      l1 = lx
      m0 = 1
      do i=1,nxx
        l1 = l1+1
        l2 = l1+kkx
        fac = tx(l2)-tx(l1)
        if (fac>0.0d0) then
          do m=1,nyy
            m1 = m0+nyy
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+1
          end do
        end if
      end do
      lx = lx+1
      kkx = kkx-1
    end do
  end if
  if(nuy/=0) then
    ly = 1
    do j=1,nuy
      ak = kky
      nyy = nyy-1
      l1 = ly
      do i=1,nyy
        l1 = l1+1
        l2 = l1+kky
        fac = ty(l2)-ty(l1)
        if(fac>0.0d0) then
          m0 = i
          do m=1,nxx
            m1 = m0+1
            wrk(m0) = (wrk(m1)-wrk(m0))*ak/fac
            m0  = m0+nky1
          end do
        end if
      end do
      ly = ly+1
      kky = kky-1
    end do
    m0 = nyy
    m1 = nky1
    do m=2,nxx
      do i=1,nyy
        m0 = m0+1
        m1 = m1+1
        wrk(m0) = wrk(m1)
      end do
      m1 = m1+nuy
    end do
  end if
!  we partition the working space and evaluate the partial derivative
  iwx = 1+nxx*nyy
  iwy = iwx+mx*(kx1-nux)
  call fpbisp(tx(nux+1:),nx-2*nux,ty(nuy+1:),ny-2*nuy,wrk(:iwx-1),kkx,kky,&
   x,mx,y,my,z,wrk(iwx:iwy-1),wrk(iwy:),iwrk(:mx),iwrk(mx+1:))
end subroutine parder
!-------------------------------------------------------
subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,&
   n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)    :: iopt, m, k, nest, maxit, k1, k2
  integer, intent(out)   :: ier
  integer, intent(inout) :: n, nrdata(nest)
  real(kind=r8), intent(in)    :: xb, xe, s, tol, x(m), y(m), w(m)
  real(kind=r8), intent(out)   :: fp, c(nest)
  real(kind=r8), intent(inout) :: t(nest), fpint(nest), z(nest),        &
                    a(nest,k1), b(nest,k2), g(nest,k2), q(m,k1)
! local variables
  real(kind=r8) :: acc, cs, fpart, fpms, fpold, fp0, f1, f2, f3, p, pinv,&
                   piv, p1, p2, p3, rn, sn, store, term, wi, xi, yi, h(7)
  integer :: i, ich1, ich3, it, iter, i1, i2, i3, j, k3, l, l0, mk1, new,&
             nk1, nmax, nmin, nplus, npl1, nrint, n8
  logical :: rstart
! parameters
  real(kind=r8), parameter :: one=1.0d0, con1=0.1d0, con9=0.9d0,        &
                              con4=0.4d-01, half=0.5d0
! intrinsic functions
!  integer, intrinsic :: max, min
!  real, intrinsic :: abs
! external functions/subroutines
  interface
    subroutine fpback(a,z,n,k,c,nest)
      use kinds, only : r8
      implicit none
      integer, intent(in)          :: n, k, nest
      real(kind=r8), intent(in)    :: a(nest,k), z(n)
      real(kind=r8), intent(inout) :: c(n)
    end subroutine fpback
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
    subroutine fpdisc(t,n,k2,b,nest)
      use kinds, only : r8
      implicit none
      integer      , intent(in)  :: n, k2, nest
      real(kind=r8), intent(in)  :: t(n)
      real(kind=r8), intent(out) :: b(nest,k2)
    end subroutine fpdisc
    subroutine fpgivs(piv,ww,cs,sn)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: piv
      real(kind=r8), intent(out)   :: cs, sn
      real(kind=r8), intent(inout) :: ww
    end subroutine fpgivs
    subroutine fpknot(x,m,t,n,fpint,nrdata,nrint,nest,istart)
      use kinds, only : r8
      implicit none
      integer      , intent(in)    :: m, nest, istart
      integer      , intent(inout) :: n, nrint, nrdata(nest)
      real(kind=r8), intent(in)    :: x(m)
      real(kind=r8), intent(inout) :: t(nest), fpint(nest)
    end subroutine fpknot
    function fprati(p1,f1,p2,f2,p3,f3)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: fprati
      real(kind=r8), intent(in)    :: p2, f2
      real(kind=r8), intent(inout) :: p1, f1, p3, f3
    end function fprati
    subroutine fprota(cs,sn,a,b)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: cs, sn
      real(kind=r8), intent(inout) :: a, b
    end subroutine fprota
  end interface
! executable statements
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 1: determination of the number of knots and their position     c
!  **************************************************************      c
!  given a set of knots we compute the least-squares spline sinf(x),   c
!  and the corresponding sum of squared residuals fp=f(p=inf).         c
!  if iopt=-1 sinf(x) is the requested approximation.                  c
!  if iopt=0 or iopt=1 we check whether we can accept the knots:       c
!if fp <=s we will continue with the current set of knots.             c
!if fp > s we will increase the number of knots and compute the        c
!   corresponding least-squares spline until finally fp<=s.            c
!the initial choice of knots depends on the value of s and iopt.       c
!if s=0 we have spline interpolation; in that case the number of       c
!knots equals nmax = m+k+1.                                            c
!if s > 0 and                                                          c
!  iopt=0 we first compute the least-squares polynomial of             c
!  degree k; n = nmin = 2*k+2                                          c
!  iopt=1 we start with the set of knots found at the last             c
!  call of the routine, except for the case that s > fp0; then         c
!  we compute directly the least-squares polynomial of degree k.       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  determine nmin, the number of knots for polynomial approximation.
  nmin = 2*k1
  if (iopt>=0) then
!  calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol*s
!  determine nmax, the number of knots for spline interpolation.
    nmax = m+k1
    if (s>0.0d0) then
!  if s>0 our initial choice of knots depends on the value of iopt.
!  if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
!  polynomial of degree k which is a spline without interior knots.
!  if iopt=1 and fp0>s we start computing the least squares spline
!  according to the set of knots found at the last call of the routine.
      if (iopt==0 .or. n==nmin) then
        n = nmin
        fpold = 0.0d0
        nplus = 0
        nrdata(1) = m-2
      else
        fp0 = fpint(n)
        if (s>=fp0) then
          n = nmin
          fpold = 0.0d0
          nplus = 0
          nrdata(1) = m-2
        else
          fpold = fpint(n-1)
          nplus = nrdata(n)
        end if
      end if
    else
!  if s=0, s(x) is an interpolating spline.
!  test whether the required storage space exceeds the available one.
      n = nmax
      if (nmax>nest) then
        ier = 1
        return
      end if
!  find the position of the interior knots in case of interpolation.
      mk1 = m-k1
      if (mk1/=0) then
        k3 = k/2
        i = k2
        j = k3+2
        if (k3*2==k) then
          do l=1,mk1
            t(i) = (x(j)+x(j-1))*half
            i = i+1
            j = j+1
          end do
        else
          do l=1,mk1
            t(i) = x(j)
            i = i+1
            j = j+1
          end do
        end if
      end if
    end if
  end if
  do
    rstart=.false.
!  main loop for the different sets of knots. m is a save upper bound
!  for the number of trials.
    do iter = 1,m
      if(n==nmin) ier = -2
!  find nrint, tne number of knot intervals.
      nrint = n-nmin+1
!  find the position of the additional knots which are needed for
!  the b-spline representation of s(x).
      nk1 = n-k1
      i = n
      do j=1,k1
        t(j) = xb
        t(i) = xe
        i = i-1
      end do
!  compute the b-spline coefficients of the least-squares spline
!  sinf(x). the observation matrix a is built up row by row and
!  reduced to upper triangular form by givens transformations.
!  at the same time fp=f(p=inf) is computed.
      fp = 0.0d0
!  initialize the observation matrix a.
      do i=1,nk1
        z(i) = 0.0d0
        do j=1,k1
          a(i,j) = 0.0d0
        end do
      end do
      l = k1
      do it=1,m
!  fetch the current data point x(it),y(it).
        xi = x(it)
        wi = w(it)
        yi = y(it)*wi
!  search for knot interval t(l) <= xi < t(l+1).
        do while (xi>=t(l+1) .and. l/=nk1)
          l = l+1
        end do
!  evaluate the (k+1) non-zero b-splines at xi and store them in q.
        call fpbspl(t,n,k,xi,l,h(:6))
        do i=1,k1
          q(it,i) = h(i)
          h(i) = h(i)*wi
        end do
!  rotate the new row of the observation matrix into triangle.
        j = l-k1
        do i=1,k1
          j = j+1
          piv = h(i)
          if (piv/=0.0d0) then
!  calculate the parameters of the givens transformation.
            call fpgivs(piv,a(j,1),cs,sn)
!  transformations to right hand side.
            call fprota(cs,sn,yi,z(j))
            if(i==k1) exit
            i2 = 1
            i3 = i+1
            do i1 = i3,k1
              i2 = i2+1
!  transformations to left hand side.
              call fprota(cs,sn,h(i1),a(j,i2))
            end do
          end if
        end do
!  add contribution of this row to the sum of squares of residual
!  right hand sides.
        fp = fp+yi**2
      end do
      if(ier==(-2)) fp0 = fp
      fpint(n) = fp0
      fpint(n-1) = fpold
      nrdata(n) = nplus
!  backward substitution to obtain the b-spline coefficients.
      call fpback(a,z,nk1,k1,c,nest)
!  test whether the approximation sinf(x) is an acceptable solution.
      if(iopt<0) return
      fpms = fp-s
      if (abs(fpms)<acc) return
!  if f(p=inf) < s accept the choice of knots.
      if (fpms<0.0d0) exit
!  if n = nmax, sinf(x) is an interpolating spline.
      if (n==nmax) then
        ier = -1
        return
      end if
!  increase the number of knots.
!  if n=nest we cannot increase the number of knots because of
!  the storage capacity limitation.
      if (n==nest) then
        ier = 1
        return
      end if
!  determine the number of knots nplus we are going to add.
      if (ier==0) then
        npl1 = nplus*2
        rn = nplus
        if (fpold-fp>acc) npl1 = rn*fpms/(fpold-fp)
        nplus = min(nplus*2,max(npl1,nplus/2,1))
      else
        nplus = 1
        ier = 0
      end if
      fpold = fp
!  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
!  t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
      fpart = 0.0d0
      i = 1
      l = k2
      new = 0
      do it=1,m
        if (x(it)>=t(l) .and. l<=nk1) then
          new = 1
          l = l+1
        end if
        term = 0.0d0
        l0 = l-k2
        do j=1,k1
          l0 = l0+1
          term = term+c(l0)*q(it,j)
        end do
        term = (w(it)*(term-y(it)))**2
        fpart = fpart+term
        if (new/=0) then
          store = term*half
          fpint(i) = fpart-store
          i = i+1
          fpart = store
          new = 0
        end if
      end do
      fpint(nrint) = fpart
      do l=1,nplus
!  add a new knot.
        call fpknot(x,m,t,n,fpint,nrdata,nrint,nest,1)
!  if n=nmax we locate the knots as for interpolation.
        if (n==nmax) then
          rstart=.true.
          exit
        end if
!  test whether we cannot further increase the number of knots.
        if (n==nest) exit
      end do
      if (rstart) exit
!  restart the computations with the new set of knots.
    end do
    if (rstart) then
!  find the position of the interior knots in case of interpolation.
      mk1 = m-k1
      if (mk1/=0) then
        k3 = k/2
        i = k2
        j = k3+2
        if (k3*2==k) then
          do l=1,mk1
            t(i) = (x(j)+x(j-1))*half
            i = i+1
            j = j+1
          end do
        else
          do l=1,mk1
            t(i) = x(j)
            i = i+1
            j = j+1
          end do
        end if
      end if
    else
      exit
    end if
  end do
!  test whether the least-squares kth degree polynomial is a solution
!  of our approximation problem.
  if (ier==(-2)) return
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  part 2: determination of the smoothing spline sp(x).                c
!  ***************************************************                 c
!  we have determined the number of knots and their position.          c
!  we now compute the b-spline coefficients of the smoothing spline    c
!  sp(x). the observation matrix a is extended by the rows of matrix   c
!  b expressing that the kth derivative discontinuities of sp(x) at    c
!  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
!  ponding weights of these additional rows are set to 1/p.            c
!  iteratively we then have to determine the value of p such that      c
!  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
!  the least-squares kth degree polynomial corresponds to p=0, and     c
!  that the least-squares spline corresponds to p=infinity. the        c
!  iteration process which is proposed here, makes use of rational     c
!  interpolation. since f(p) is a convex and strictly decreasing       c
!  function of p, it can be approximated by a rational function        c
!  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
!  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
!  to calculate the new value of p such that r(p)=s. convergence is    c
!  guaranteed by taking f1>0 and f3<0.                                 c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  evaluate the discontinuity jump of the kth derivative of the
!  b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
  call fpdisc(t,n,k2,b,nest)
!  initial value for p.
  p1 = 0.0d0
  f1 = fp0-s
  p3 = -one
  f3 = fpms
  p = 0.0d0
  do i=1,nk1
    p = p+a(i,1)
  end do
  rn = nk1
  p = rn/p
  ich1 = 0
  ich3 = 0
  n8 = n-nmin
!  iteration process to find the root of f(p) = s.
  do iter=1,maxit
!  the rows of matrix b with weight 1/p are rotated into the
!  triangularised observation matrix a which is stored in g.
    pinv = one/p
    do i=1,nk1
      c(i) = z(i)
      g(i,k2) = 0.0d0
      do j=1,k1
        g(i,j) = a(i,j)
      end do
    end do
    do it=1,n8
!  the row of matrix b is rotated into triangle by givens transformation
      do i=1,k2
        h(i) = b(it,i)*pinv
      end do
      yi = 0.0d0
      do j=it,nk1
        piv = h(1)
!  calculate the parameters of the givens transformation.
        call fpgivs(piv,g(j,1),cs,sn)
!  transformations to right hand side.
        call fprota(cs,sn,yi,c(j))
        if (j==nk1) exit
        i2 = k1
        if (j>n8) i2 = nk1-j
        do i=1,i2
!  transformations to left hand side.
          i1 = i+1
          call fprota(cs,sn,h(i1),g(j,i1))
          h(i) = h(i1)
        end do
        h(i2+1) = 0.0d0
      end do
    end do
!  backward substitution to obtain the b-spline coefficients.
    call fpback(g,c,nk1,k2,c,nest)
!  computation of f(p).
    fp = 0.0d0
    l = k2
    do it=1,m
      if (x(it)>=t(l) .and. l<=nk1) l = l+1
      l0 = l-k2
      term = 0.0d0
      do j=1,k1
        l0 = l0+1
        term = term+c(l0)*q(it,j)
      end do
      fp = fp+(w(it)*(term-y(it)))**2
    end do
!  test whether the approximation sp(x) is an acceptable solution.
    fpms = fp-s
    if (abs(fpms)<acc) return
!  test whether the maximal number of iterations is reached.
    if (iter==maxit) exit
!  carry out one more step of the iteration process.
    p2 = p
    f2 = fpms
    if (ich3==0) then
      if ((f2-f3)<=acc) then
!  our initial choice of p is too large.
        p3 = p2
        f3 = f2
        p = p*con4
        if (p<=p1) p=p1*con9 + p2*con1
        cycle
      end if
      if (f2<0.0d0) ich3=1
    end if
    if (ich1==0) then
      if ((f1-f2)<=acc) then
!  our initial choice of p is too small
        p1 = p2
        f1 = f2
        p = p/con4
        if (p3<0.0d0) cycle
        if (p>=p3) p = p2*con1 + p3*con9
        cycle
      end if
      if (f2>0.0d0) ich1=1
    end if
!  test whether the iteration process proceeds as theoretically
!  expected.
    if (f2>=f1 .or. f2<=f3) then
      ier = 2
      return
    end if
!  find the new value for p.
    p = fprati(p1,f1,p2,f2,p3,f3)
  end do
  ier = 3
end subroutine fpcurf
!-------------------------------------------------------
subroutine curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,&
   wrk,lwrk,iwrk,ier)
!  given the set of data points (x(i),y(i)) and the set of positive
!  numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
!  approximation of degree k on the interval xb <= x <= xe.
!  if iopt=-1 curfit calculates the weighted least-squares spline
!  according to a given set of knots.
!  if iopt>=0 the number of knots of the spline s(x) and the position
!  t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
!  ness of s(x) is then achieved by minimalizing the discontinuity
!  jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
!  n-k-1. the amount of smoothness is determined by the condition that
!  f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
!  negative constant, called the smoothing factor.
!  the fit s(x) is given in the b-spline representation (b-spline coef-
!  ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
!subroutine splev.
!
!  calling sequence:
! call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
!* lwrk,iwrk,ier)
!
!  parameters:
!   iopt  : integer flag. on entry iopt must specify whether a weighted
!       least-squares spline (iopt=-1) or a smoothing spline (iopt=
!       0 or 1) must be determined. if iopt=0 the routine will start
!       with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
!       k+1. if iopt=1 the routine will continue with the knots
!       found at the last call of the routine.
!       attention: a call with iopt=1 must always be immediately
!       preceded by another call with iopt=1 or iopt=0.
!       unchanged on exit.
!   m     : integer. on entry m must specify the number of data points.
!       m > k. unchanged on exit.
!   x     : real array of dimension at least (m). before entry, x(i)
!       must be set to the i-th value of the independent variable x,
!       for i=1,2,...,m. these values must be supplied in strictly
!       ascending order. unchanged on exit.
!   y     : real array of dimension at least (m). before entry, y(i)
!       must be set to the i-th value of the dependent variable y,
!       for i=1,2,...,m. unchanged on exit.
!   w     : real array of dimension at least (m). before entry, w(i)
!       must be set to the i-th value in the set of weights. the
!       w(i) must be strictly positive. unchanged on exit.
!       see also further comments.
!   xb,xe : real values. on entry xb and xe must specify the boundaries
!       of the approximation interval. xb<=x(1), xe>=x(m).
!       unchanged on exit.
!   k     : integer. on entry k must specify the degree of the spline.
!       1<=k<=5. it is recommended to use cubic splines (k=3).
!       the user is strongly dissuaded from choosing k even,together
!       with a small s-value. unchanged on exit.
!   s     : real.on entry (in case iopt>=0) s must specify the smoothing
!       factor. s >=0. unchanged on exit.
!       for advice on the choice of s see further comments.
!   nest  : integer. on entry nest must contain an over-estimate of the
!       total number of knots of the spline returned, to indicate
!       the storage space available to the routine. nest >=2*k+2.
!       in most practical situation nest=m/2 will be sufficient.
!       always large enough is  nest=m+k+1, the number of knots
!       needed for interpolation (s=0). unchanged on exit.
!   n     : integer.
!       unless ier =10 (in case iopt >=0), n will contain the
!       total number of knots of the spline approximation returned.
!       if the computation mode iopt=1 is used this value of n
!       should be left unchanged between subsequent calls.
!       in case iopt=-1, the value of n must be specified on entry.
!   t     : real array of dimension at least (nest).
!       on succesful exit, this array will contain the knots of the
!       spline,i.e. the position of the interior knots t(k+2),t(k+3)
!       ...,t(n-k-1) as well as the position of the additional knots
!       t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
!       the b-spline representation.
!       if the computation mode iopt=1 is used, the values of t(1),
!       t(2),...,t(n) should be left unchanged between subsequent
!       calls. if the computation mode iopt=-1 is used, the values
!       t(k+2),...,t(n-k-1) must be supplied by the user, before
!       entry. see also the restrictions (ier=10).
!   c     : real array of dimension at least (nest).
!       on succesful exit, this array will contain the coefficients
!       c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
!   fp    : real. unless ier=10, fp contains the weighted sum of
!       squared residuals of the spline approximation returned.
!   wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
!       used as working space. if the computation mode iopt=1 is
!       used, the values wrk(1),...,wrk(n) should be left unchanged
!       between subsequent calls.
!   lwrk  : integer. on entry,lwrk must specify the actual dimension of
!       the array wrk as declared in the calling (sub)program.lwrk
!       must not be too small (see wrk). unchanged on exit.
!   iwrk  : integer array of dimension at least (nest).
!       used as working space. if the computation mode iopt=1 is
!       used,the values iwrk(1),...,iwrk(n) should be left unchanged
!       between subsequent calls.
!   ier   : integer. unless the routine detects an error, ier contains a
!       non-positive value on exit, i.e.
!ier=0  : normal return. the spline returned has a residual sum of
!         squares fp such that abs(fp-s)/s <= tol with tol a relat-
!         ive tolerance set to 0.001 by the program.
!ier=-1 : normal return. the spline returned is an interpolating
!         spline (fp=0).
!ier=-2 : normal return. the spline returned is the weighted least-
!         squares polynomial of degree k. in this extreme case fp
!         gives the upper bound fp0 for the smoothing factor s.
!ier=1  : error. the required storage space exceeds the available
!         storage space, as specified by the parameter nest.
!         probably causes : nest too small. if nest is already
!         large (say nest > m/2), it may also indicate that s is
!         too small
!         the approximation returned is the weighted least-squares
!         spline according to the knots t(1),t(2),...,t(n). (n=nest)
!         the parameter fp gives the corresponding weighted sum of
!         squared residuals (fp>s).
!ier=2  : error. a theoretically impossible result was found during
!         the iteration proces for finding a smoothing spline with
!         fp = s. probably causes : s too small.
!         there is an approximation returned but the corresponding
!         weighted sum of squared residuals does not satisfy the
!         condition abs(fp-s)/s < tol.
!ier=3  : error. the maximal number of iterations maxit (set to 20
!         by the program) allowed for finding a smoothing spline
!         with fp=s has been reached. probably causes : s too small
!         there is an approximation returned but the corresponding
!         weighted sum of squared residuals does not satisfy the
!         condition abs(fp-s)/s < tol.
!ier=10 : error. on entry, the input data are controlled on validity
!         the following restrictions must be satisfied.
!         -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
!         xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
!         if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
!                     xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
!                   the schoenberg-whitney conditions, i.e. there
!                   must be a subset of data points xx(j) such that
!                     t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
!         if iopt>=0: s>=0
!                     if s=0 : nest >= m+k+1
!         if one of these conditions is found to be violated,control
!         is immediately repassed to the calling program. in that
!         case there is no approximation returned.
!
!  further comments:
!   by means of the parameter s, the user can control the tradeoff
!   between closeness of fit and smoothness of fit of the approximation.
!   if s is too large, the spline will be too smooth and signal will be
!   lost ; if s is too small the spline will pick up too much noise. in
!   the extreme cases the program will return an interpolating spline if
!   s=0 and the weighted least-squares polynomial of degree k if s is
!   very large. between these extremes, a properly chosen s will result
!   in a good compromise between closeness of fit and smoothness of fit.
!   to decide whether an approximation, corresponding to a certain s is
!   satisfactory the user is highly recommended to inspect the fits
!   graphically.
!   recommended values for s depend on the weights w(i). if these are
!   taken as 1/d(i) with d(i) an estimate of the standard deviation of
!   y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
!   sqrt(2*m)). if nothing is known about the statistical error in y(i)
!   each w(i) can be set equal to one and s determined by trial and
!   error, taking account of the comments above. the best is then to
!   start with a very large value of s ( to determine the least-squares
!   polynomial and the corresponding upper bound fp0 for s) and then to
!   progressively decrease the value of s ( say by a factor 10 in the
!   beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
!   approximation shows more detail) to obtain closer fits.
!   to economize the search for a good s-value the program provides with
!   different modes of computation. at the first call of the routine, or
!   whenever he wants to restart with the initial set of knots the user
!   must set iopt=0.
!   if iopt=1 the program will continue with the set of knots found at
!   the last call of the routine. this will save a lot of computation
!   time if curfit is called repeatedly for different values of s.
!   the number of knots of the spline returned and their location will
!   depend on the value of s and on the complexity of the shape of the
!   function underlying the data. but, if the computation mode iopt=1
!   is used, the knots returned may also depend on the s-values at
!   previous calls (if these were smaller). therefore, if after a number
!   of trials with different s-values and iopt=1, the user can finally
!   accept a fit as satisfactory, it may be worthwhile for him to call
!   curfit once more with the selected value for s but now with iopt=0.
!   indeed, curfit may then return an approximation of the same quality
!   of fit but with fewer knots and therefore better if data reduction
!   is also an important objective for the user.
!
!  other subroutines required:
!fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota
!
!  references:
!   dierckx p. : an algorithm for smoothing, differentiation and integ-
!            ration of experimental data using spline functions,
!            j.comp.appl.maths 1 (1975) 165-184.
!   dierckx p. : a fast algorithm for smoothing data on a rectangular
!            grid while using spline functions, siam j.numer.anal.
!            19 (1982) 1286-1304.
!   dierckx p. : an improved algorithm for curve fitting with spline
!            functions, report tw54, dept. computer science,k.u.
!            leuven, 1981.
!   dierckx p. : curve and surface fitting with splines, monographs on
!            numerical analysis, oxford university press, 1993.
!
!  author:
!p.dierckx
!dept. computer science, k.u. leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  creation date : may 1979
!  latest update : march 1987
!
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)  :: iopt, m, k, nest, lwrk
  integer, intent(out) :: ier
  integer, intent(inout) :: n, iwrk(nest)
  real(kind=r8), intent(in)    :: xb, xe, s, x(m), y(m), w(m)
  real(kind=r8), intent(out)   :: fp, c(nest)
  real(kind=r8), intent(inout) :: t(nest), wrk(lwrk)
! local variables
  integer :: i, ia, ib, ifp, ig, iq, iz, j, k1, k2, lwest, nmin
! parameters
  integer, parameter :: maxit = 20
  real(kind=r8), parameter :: tol = 0.1d-02
! external subroutines
  interface
    subroutine fpchec(x,m,t,n,k,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)       :: m, n, k
      real(kind=r8), intent(in) :: x(m), t(n)
      integer, intent(out)      :: ier
    end subroutine fpchec
    subroutine fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,&
       n,t,c,fp,fpint,z,a,b,g,q,nrdata,ier)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: iopt, m, k, nest, maxit, k1, k2
      integer, intent(out)   :: ier
      integer, intent(inout) :: n, nrdata(nest)
      real(kind=r8), intent(in)    :: xb, xe, s, tol, x(m), y(m), w(m)
      real(kind=r8), intent(out)   :: fp, c(nest)
      real(kind=r8), intent(inout) :: t(nest), fpint(nest), z(nest),        &
                        a(nest,k1), b(nest,k2), g(nest,k2), q(m,k1)
    end subroutine fpcurf
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
  ier = 10
  if (k<=0 .or. k>5) return
  k1 = k+1
  k2 = k1+1
  if (iopt<(-1) .or. iopt>1) return
  nmin = 2*k1
  if (m<k1 .or. nest<nmin) return
  lwest = m*k1+nest*(7+3*k)
  if (lwrk<lwest) return
  if (xb>x(1) .or. xe<x(m) .or. w(1)<=0.0d0) return
  do i=2,m
    if (x(i-1)>=x(i) .or. w(i)<=0.0d0) return
  end do
  if (iopt>=0) then
    if (s<0.0d0) return
    if (s==0.0d0 .and. nest<(m+k1)) return
    ier = 0
  else
    if (n<nmin .or. n>nest) return
    j = n
    do i=1,k1
      t(i) = xb
      t(j) = xe
      j = j-1
    end do
    call fpchec(x,m,t,n,k,ier)
    if (ier/=0) return
  end if
! we partition the working space and determine the spline approximation.
  ifp = 1
  iz = ifp+nest
  ia = iz+nest
  ib = ia+nest*k1
  ig = ib+nest*k2
  iq = ig+nest*k2
  call fpcurf(iopt,x,y,w,m,xb,xe,k,s,nest,tol,maxit,k1,k2,n,t,c,fp,     &
   wrk(ifp:iz-1),wrk(iz:ia-1),wrk(ia:ib-1),wrk(ib:ig-1),wrk(ig:iq-1),   &
   wrk(iq:),iwrk,ier)
end subroutine curfit
!-------------------------------------------------------
subroutine splder(t,n,c,k,nu,x,y,m,wrk,ier)
!subroutine splder evaluates in a number of points x(i),i=1,2,...,m
!  the derivative of order nu of a spline s(x) of degree k,given in
!  its b-spline representation.
!
!  calling sequence:
! call splder(t,n,c,k,nu,x,y,m,wrk,ier)
!
!  input parameters:
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!nu   : integer, specifying the order of the derivative. 0<=nu<=k
!x    : array,length m, which contains the points where the deriv-
!       ative of s(x) must be evaluated.
!m    : integer, giving the number of points where the derivative
!       of s(x) must be evaluated
!wrk  : real array of dimension n. used as working space.
!
!  output parameters:
!y    : array,length m, giving the value of the derivative of s(x)
!       at the different points.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!0 <= nu <= k
!m >= 1
!t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl
!
!  references :
!de boor c : on calculating with b-splines, j. approximation theory
!            6 (1972) 50-62.
!cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
!            applics 10 (1972) 134-149.
!   dierckx p. : curve and surface fitting with splines, monographs on
!            numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)  :: n, k, nu, m
  integer, intent(out) :: ier
  real(kind=r8), intent(in)    :: t(n), c(n), x(m)
  real(kind=r8), intent(out)   :: y(m), wrk(n)
! local variables
  integer       :: i, j, kk, k1, k2, l, ll, l1, l2, nk1, nk2
  real(kind=r8) :: ak, arg, fac, sp, tb, te, h(6)
! external subroutines
  interface
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
  ier = 10
  if (nu<0 .or. nu>k) return
  if (m<1) return
  do i=2,m
    if(x(i)<x(i-1)) return
  end do
  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
  k1 = k+1
  nk1 = n-k1
  tb = t(k1)
  te = t(nk1+1)
!  the derivative of order nu of a spline of degree k is a spline of
!  degree k-nu,the b-spline coefficients wrk(i) of which can be found
!  using the recurrence scheme of de boor.
  l = 1
  kk = k
  do i=1,nk1
    wrk(i) = c(i)
  end do
  if (nu/=0) then
    nk2 = nk1
    do j=1,nu
      ak = kk
      nk2 = nk2-1
      l1 = l
      do i=1,nk2
        l1 = l1+1
        l2 = l1+kk
        fac = t(l2)-t(l1)
        if (fac>0.0d0) wrk(i) = ak*(wrk(i+1)-wrk(i))/fac
      end do
      l = l+1
      kk = kk-1
    end do
    if (kk==0) then
!  if nu=k the derivative is a piecewise constant function
      j = 1
      do i=1,m
        arg = x(i)
        do while (arg>=t(l+1) .and. l/=nk1)
          l = l+1
          j = j+1
        end do
        y(i) = wrk(j)
      end do
      return
    end if
  end if
  l = k1
  l1 = l+1
  k2 = k1-nu
!  main loop for the different points.
  do i=1,m
!  fetch a new x-value arg.
    arg = x(i)
    if(arg<tb) arg = tb
    if(arg>te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
    do while (arg>=t(l1) .and. l/=nk1)
      l = l1
      l1 = l+1
    end do
!  evaluate the non-zero b-splines of degree k-nu at arg.
    call fpbspl(t,n,kk,arg,l,h)
!  find the value of the derivative at x=arg.
    sp = 0.0d0
    ll = l-k1
    do j=1,k2
      ll = ll+1
      sp = sp+wrk(ll)*h(j)
    end do
    y(i) = sp
  end do
end subroutine splder
!-------------------------------------------------------
subroutine splev(t,n,c,k,x,y,m,ier)
!subroutine splev evaluates in a number of points x(i),i=1,2,...,m
!  a spline s(x) of degree k, given in its b-spline representation.
!
!  calling sequence:
! call splev(t,n,c,k,x,y,m,ier)
!
!  input parameters:
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!x    : array,length m, which contains the points where s(x) must
!       be evaluated.
!m    : integer, giving the number of points where s(x) must be
!       evaluated.
!
!  output parameter:
!y    : array,length m, giving the value of s(x) at the different
!       points.
!ier  : error flag
!  ier = 0 : normal return
!  ier =10 : invalid input data (see restrictions)
!
!  restrictions:
!m >= 1
!t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.
!
!  other subroutines required: fpbspl.
!
!  references :
!de boor c  : on calculating with b-splines, j. approximation theory
!             6 (1972) 50-62.
!cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!             applics 10 (1972) 134-149.
!dierckx p. : curve and surface fitting with splines, monographs on
!             numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)  :: n, k, m
  integer, intent(out) :: ier
  real(kind=r8), intent(in)  :: t(n), c(n), x(m)
  real(kind=r8), intent(out) :: y(m)
! local variables
  integer       :: i, j, k1, l, ll, l1, nk1
  real(kind=r8) :: arg, sp, tb, te, h(6)
! external subroutines
  interface
    subroutine fpbspl(t,n,k,x,l,h)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, k, l
      real(kind=r8), intent(in) :: x, t(n)
      real(kind=r8), intent(out) :: h(6)
    end subroutine fpbspl
  end interface
! executable statements
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
  ier = 10
  if (m<1) return
  do i=2,m
    if(x(i)<x(i-1)) return
  end do
  ier = 0
!  fetch tb and te, the boundaries of the approximation interval.
  k1 = k+1
  nk1 = n-k1
  tb = t(k1)
  te = t(nk1+1)
  l = k1
  l1 = l+1
!  main loop for the different points.
  do i=1,m
!  fetch a new x-value arg.
    arg = x(i)
    if(arg<tb) arg = tb
    if(arg>te) arg = te
!  search for knot interval t(l) <= arg < t(l+1)
    do while (arg>=t(l1) .and. l/=nk1)
      l = l1
      l1 = l+1
    end do
!  evaluate the non-zero b-splines at arg.
    call fpbspl(t,n,k,arg,l,h)
!  find the value of s(x) at x=arg.
    sp = 0.0d0
    ll = l-k1
    do j=1,k1
      ll = ll+1
      sp = sp+c(ll)*h(j)
    end do
    y(i) = sp
  end do
end subroutine splev
!-------------------------------------------------------
subroutine fpintb(t,n,bint,nk1,x,y)
!subroutine fpintb calculates integrals of the normalized b-splines
!  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
!  it makes use of the formulae of gaffney for the calculation of
!  indefinite integrals of b-splines.
!
!  calling sequence:
! call fpintb(t,n,bint,nk1,x,y)
!
!  input parameters:
!t    : real array,length n, containing the position of the knots.
!n    : integer value, giving the number of knots.
!nk1  : integer value, giving the number of b-splines of degree k,
!       defined on the set of knots ,i.e. nk1 = n-k-1.
!x,y  : real values, containing the end points of the integration
!       interval.
!  output parameter:
!bint : array,length nk1, containing the integrals of the b-splines.
!  ..
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in) :: n, nk1
  real(kind=r8), intent(in)  :: x, y, t(n)
  real(kind=r8), intent(out) :: bint(nk1)
! local variables
  integer :: i, ia, ib, it, j, j1, k, k1, l, li, lj, lk, l0, min
  real(kind=r8) :: a, ak, arg, b, f
  real(kind=r8), dimension(6) :: aaint, h, h1
! parameters
  real(kind=r8), parameter :: one = 1.0d0
! executable statements
  k1 = n-nk1
  ak = k1
  k = k1-1
  do i=1,nk1
    bint(i) = 0.0d0
  end do
!  the integration limits are arranged in increasing order.
  a = x
  b = y
  min = 0
  if (a==b) return
  if (a>b) then
    a = y
    b = x
    min = 1
  end if
  if (a<t(k1)) a = t(k1)
  if (b>t(nk1+1)) b = t(nk1+1)
!  using the expression of gaffney for the indefinite integral of a
!  b-spline we find that
!  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
!where for t(l) <= x < t(l+1)
!res(j,x) = 0, j=1,2,...,l-k-1
!         = 1, j=l+1,l+2,...,nk1
!         = aaint(j+k-l+1), j=l-k,l-k+1,...,l
!           = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
!             i=0,1,...,k
  l = k1
  l0 = l+1
!  set arg = a.
  arg = a
  do it=1,2
!  search for the knot interval t(l) <= arg < t(l+1).
    do while (arg>=t(l0) .and. l/=nk1)
      l = l0
      l0 = l+1
    end do
!  calculation of aaint(j), j=1,2,...,k+1.
!  initialization.
    do j=1,k1
      aaint(j) = 0.0d0
    end do
    aaint(1) = (arg-t(l))/(t(l+1)-t(l))
    h1(1) = one
    do j=1,k
!  evaluation of the non-zero b-splines of degree j at arg,i.e.
!h(i+1) = nl-j+i,j(arg), i=0,1,...,j.
      h(1) = 0.0d0
      do i=1,j
        li = l+i
        lj = li-j
        f = h1(i)/(t(li)-t(lj))
        h(i) = h(i)+f*(t(li)-arg)
        h(i+1) = f*(arg-t(lj))
      end do
!  updating of the integrals aaint.
      j1 = j+1
      do i=1,j1
        li = l+i
        lj = li-j1
        aaint(i) = aaint(i)+h(i)*(arg-t(lj))/(t(li)-t(lj))
        h1(i) = h(i)
      end do
    end do
    if (it==2) exit
!  updating of the integrals bint
    lk = l-k
    ia = lk
    do i=1,k1
      bint(lk) = -aaint(i)
      lk = lk+1
    end do
!  set arg = b.
    arg = b
  end do
!  updating of the integrals bint.
  lk = l-k
  ib = lk-1
  do i=1,k1
    bint(lk) = bint(lk)+aaint(i)
    lk = lk+1
  end do
  do i=ia,ib
    bint(i) = bint(i)+one
  end do
!  the scaling factors are taken into account.
  f = one/ak
  do i=1,nk1
    j = i+k1
    bint(i) = bint(i)*(t(j)-t(i))*f
  end do
!  the order of the integration limits is taken into account.
  if (min==0) return
  do i=1,nk1
    bint(i) = -bint(i)
  end do
end subroutine fpintb
!-------------------------------------------------------
function splint(t,n,c,k,a,b,wrk)
!  function splint calculates the integral of a spline function s(x)
!  of degree k, which is given in its normalized b-spline representation
!
!  calling sequence:
! aaint = splint(t,n,c,k,a,b,wrk)
!
!  input parameters:
!t    : array,length n,which contains the position of the knots
!       of s(x).
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, containing the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!a,b  : real values, containing the end points of the integration
!       interval. s(x) is considered to be identically zero outside
!       the interval (t(k+1),t(n-k)).
!
!  output parameter:
!aaint : real, containing the integral of s(x) between a and b.
!wrk  : real array, length n.  used as working space
!       on output, wrk will contain the integrals of the normalized
!       b-splines defined on the set of knots.
!
!  other subroutines required: fpintb.
!
!  references :
!gaffney p.w. : the calculation of indefinite integrals of b-splines
!               j. inst. maths applics 17 (1976) 37-41.
!dierckx p. : curve and surface fitting with splines, monographs on
!             numerical analysis, oxford university press, 1993.
!
!  author :
!p.dierckx
!dept. computer science, k.u.leuven
!celestijnenlaan 200a, b-3001 heverlee, belgium.
!e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
  use kinds, only : r8
  implicit none
  real(kind=r8) splint
! arguments
  integer, intent(in) :: n, k
  real(kind=r8), intent(in) :: a, b, t(n), c(n)
  real(kind=r8), intent(out) :: wrk(n)
! local variables
  integer :: i, nk1
  interface
    subroutine fpintb(t,n,bint,nk1,x,y)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, nk1
      real(kind=r8), intent(in)  :: x, y, t(n)
      real(kind=r8), intent(out) :: bint(nk1)
    end subroutine fpintb
  end interface
! executable statements
  nk1 = n-k-1
!  calculate the integrals wrk(i) of the normalized b-splines
!  ni,k+1(x), i=1,2,...nk1.
  call fpintb(t,n,wrk,nk1,a,b)
!  calculate the integral of s(x).
  splint = 0.0d0
  do i=1,nk1
    splint = splint+c(i)*wrk(i)
  end do
end function splint
!-------------------------------------------------------
subroutine profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
	use kinds, only : r8
	implicit none
!  if iopt=0 subroutine profil calculates the b-spline coefficients of
!  the univariate spline f(y) = s(u,y) with s(x,y) a bivariate spline of
!  degrees kx and ky, given in the b-spline representation.
!  if iopt = 1 it calculates the b-spline coefficients of the univariate
!  spline g(x) = s(x,u)
!
!  calling sequence:
!     call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier)
!
!  input parameters:
!   iopt  : integer flag, specifying whether the profile f(y) (iopt=0)
!           or the profile g(x) (iopt=1) must be determined.
!   tx    : real array, length nx, which contains the position of the
!           knots in the x-direction.
!   nx    : integer, giving the total number of knots in the x-direction
!   ty    : real array, length ny, which contains the position of the
!           knots in the y-direction.
!   ny    : integer, giving the total number of knots in the y-direction
!   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
!           b-spline coefficients.
!   kx,ky : integer values, giving the degrees of the spline.
!   u     : real value, specifying the requested profile.
!           tx(kx+1)<=u<=tx(nx-kx), if iopt=0.
!           ty(ky+1)<=u<=ty(ny-ky), if iopt=1.
!   nu    : on entry nu must specify the dimension of the array cu.
!           nu >= ny if iopt=0, nu >= nx if iopt=1.
!
!  output parameters:
!   cu    : real array of dimension (nu).
!           on succesful exit this array contains the b-spline
!   ier   : integer error flag
!    ier=0 : normal return
!    ier=10: invalid input data (see restrictions)
!
!  restrictions:
!   if iopt=0 : tx(kx+1) <= u <= tx(nx-kx), nu >=ny.
!   if iopt=1 : ty(ky+1) <= u <= ty(ny-ky), nu >=nx.
!
!  other subroutines required:
!    fpbspl
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
!  ..scalar arguments..
	integer,intent(in) :: iopt,nx,ny,kx,ky,nu
	integer,intent(out) :: ier
	real(kind=r8), intent(in) :: u
!  ..array arguments..
	real(kind=r8),intent(in) :: tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1))
	real(kind=r8),intent(out) :: cu(nu)
!  ..local scalars..
	integer :: i,j,kx1,ky1,l,l1,m,m0,nkx1,nky1
	real(kind=r8) :: sum
!  ..local array
	real(kind=r8) :: h(6)
!  ..
	interface
		subroutine fpbspl(t,n,k,x,l,h)
			use kinds, only : r8
			implicit none
			integer, intent(in) :: n, k, l
			real(kind=r8), intent(in) :: x, t(n)
			real(kind=r8), intent(out) :: h(6)
		end subroutine fpbspl
	end interface
!  before starting computations a data check is made. if the input data
!  are invalid control is immediately repassed to the calling program.
	kx1 = kx+1
	ky1 = ky+1
	nkx1 = nx-kx1
	nky1 = ny-ky1
	ier = 10
	if(iopt == 0) then
		if(nu.lt.ny) return
		if(u.lt.tx(kx1) .or. u.gt.tx(nkx1+1)) return
!  the b-splinecoefficients of f(y) = s(u,y).
		ier = 0
		l = kx1
		l1 = l+1
		do while (u >= tx(l1) .and. l /= nkx1)
			l = l1
			l1 = l+1
		end do
		call fpbspl(tx,nx,kx,u,l,h)
		m0 = (l-kx1)*nky1+1
		do i=1,nky1
			m = m0
			sum = 0.0d0
			do j=1,kx1
				sum = sum+h(j)*c(m)
				m = m+nky1
			end do
			cu(i) = sum
			m0 = m0+1
		end do
	else
		if(nu.lt.nx) return
		if(u.lt.ty(ky1) .or. u.gt.ty(nky1+1)) return
!  the b-splinecoefficients of g(x) = s(x,u).
		ier = 0
		l = ky1
		l1 = l+1
		do while (u >= ty(l1) .and. l /= nky1)
			l = l1
			l1 = l+1
		end do
		call fpbspl(ty,ny,ky,u,l,h)
		m0 = l-ky
		do i=1,nkx1
			m = m0
			sum = 0.0d0
			do j=1,ky1
				sum = sum+h(j)*c(m)
				m = m+1
			end do
			cu(i) = sum
			m0 = m0+nky1
		end do
	end if
end subroutine profil
!-------------------------------------------------------
integer function err_handle(ifail,ierror,srname,nrec,rec) ! P01ABF
! P01ABF either returns the value of IERROR through the routine
! name (soft failure), or terminates execution of the program
! (hard failure). Diagnostic messages may be output.
! If IERROR = 0 (successful exit from the calling routine),
! the value 0 is returned through the routine name, and no
! message is output
! If IERROR is non-zero (abnormal exit from the calling routine),
! the action taken depends on the value of IFAIL.
! IFAIL =  1: soft failure, silent exit (i.e. no messages are
!             output)
! IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
! IFAIL =-13: soft failure, noisy exit but standard messages from
!             P01ABF are suppressed
! IFAIL =  0: hard failure, noisy exit
! For compatibility with certain routines included before Mark 12
! P01ABF also allows an alternative specification of IFAIL in which
! it is regarded as a decimal integer with least significant digits
! cba. Then
! a = 0: hard failure  a = 1: soft failure
! b = 0: silent exit   b = 1: noisy exit
! except that hard failure now always implies a noisy exit.
!
  implicit none
! arguments
  integer, intent(in)          :: ierror, ifail, nrec
  character(len=*), intent(in) :: srname, rec(:)
! local variables
  integer           :: i
  character(len=72) :: mess
! intrinsic functions
!  integer, intrinsic :: mod
!  real   , intrinsic :: abs
! executable statements
  if (ierror/=0) then
!   abnormal exit from calling routine
    if (ifail==-1 .or. ifail==0 .or. ifail==-13 .or. &
       (ifail>0 .and. mod(ifail/10,10)/=0)) then
!     noisy exit
      do i = 1, nrec
        write(*,"(a)") rec(i)
      end do
      if (ifail/=-13) then
        write (mess, &
     fmt="(' ** Abnormal exit from ',a,': IFAIL =',i6)") srname, ierror
        write(*,*) mess
        if (abs(mod(ifail,10))/=1) then
!         hard failure
          write(*,*) ' ** Hard failure - execution terminated'
          stop
        else
          write(*,*) ' ** Soft failure - control returned'
        end if
      end if
    end if
  end if
  err_handle = ierror
end function err_handle
!-------------------------------------------------------
subroutine zeroin_rev(x,y,fx,tolx,ir,c,ind,ifail) ! C05AZF
! Find zero in given interval of continuous function by Bus and Dekker
! algorithm (reverse communication)
  use kinds, only : r8
  implicit none
! arguments
  real(kind=r8), intent(in)    :: tolx
  real(kind=r8), intent(inout) :: fx, x, y, c(17)
  integer, intent(in)          :: ir
  integer, intent(inout)       :: ifail, ind
! local variables
  real(kind=r8) :: ab, diff, diff1, diff2, rel, rmax, tol, tol1
  integer       :: i
  logical       :: t
! intrinsic functions
!  integer, intrinsic :: int, max
!  real   , intrinsic :: abs, sign, dble, sqrt, epsilon, tiny
! parameters
  character(len=*), parameter :: srname='ZEROIN_REV'
  character(len=1), dimension(1), parameter :: dumrec=(/' '/)
  real(kind=r8)            :: epsqrt
  real(kind=r8), parameter :: eps = epsilon(1.d0)/2.d0, &
                              small = tiny(1.d0)
! external functions
  interface
    integer function err_handle(ifail,ierror,srname,nrec,rec)
      implicit none
      integer, intent(in)          :: ierror, ifail, nrec
      character(len=*), intent(in) :: srname, rec(:)
    end function err_handle
  end interface
! executable statements
  epsqrt = sqrt(eps)
  i = 0
  if ((ind>0 .and. ind<=4) .or. ind==-1) then
    if (tolx>0.d0 .and. (ir==0 .or. ir==1 .or. ir==2)) then
      rel = 1.d0
      ab = 1.d0
      if (ir==1) then
        rel = 0.d0
      else if (ir==2) then
        ab = 0.d0
      end if
      if (ind==1) then
        c(3) = x
        ind = 2
        return
      end if
      if (ind==-1) c(3) = x
      if ((ind==-1 .or. ind==2) .and. fx/=0.d0) then
        c(4) = fx
        c(15) = abs(fx)
        c(16) = 0.d0
        x = y
        y = c(3)
        c(2) = c(4)
        c(5) = x
        if (ind==2) then
          ind = 3
          return
        end if
        fx = c(1)
        ind = 3
      end if
      if (fx/=0.d0) then
        if (ind==4 .or. sign(1.d0,fx)/=sign(1.d0,c(2))) then
          if (ind==4) then
            c(4) = fx
            rmax = abs(fx)
            if (c(13)*rmax<=c(15)) then
              c(16) = 0.d0
            else
              if (c(16)==1.d0) c(16) = -1.d0
              if (c(16)==0.d0) c(16) = 1.d0
            end if
            if (c(2)>=0.d0) then
              t = c(4) >= 0.d0
            else
              t = c(4) <= 0.d0
            end if
            if (.not.t) then
              i = int(c(17)+0.1d0)
              i = i + 1
              if (c(11)==c(12)) i = 0
              c(17) = dble(i)
            end if
          else
            c(6) = fx
            c(13) = epsqrt
            c(15) = max(c(15),abs(fx))
            c(14) = small
            c(16) = 0.0d0
          end if
          if (ind/=4 .or. t) then
            c(1) = c(5)
            c(2) = c(6)
            c(17) = 0.d0
          end if
          if (abs(c(2))<abs(c(4))) then
            if (c(1)/=c(5)) then
              c(7) = c(5)
              c(8) = c(6)
            end if
            c(5) = c(3)
            c(6) = c(4)
            x = c(1)
            c(3) = x
            c(4) = c(2)
            c(1) = c(5)
            c(2) = c(6)
          end if
          tol = 0.5d0*tolx*max(ab,rel*abs(c(3)))
          tol1 = 2.0d0*eps*max(ab,rel*abs(c(3)))
          diff2 = 0.5d0*(c(1)-c(3))
          c(12) = diff2
          diff2 = diff2 + c(3)
          if (abs(c(12))<=tol) then
            if (c(16)<0.d0) then
              i = 4
            else
              y = c(1)
              i = 0
            end if
            ind = 0
          else
            if (abs(c(12))>tol1) then
              if (c(17)<2.5d0) then
                tol = tol*sign(1.d0,c(12))
                diff1 = (c(3)-c(5))*c(4)
              end if
              if (.not.(c(17)<2.5d0 .and. c(17)>1.5d0 &
                  .and. (c(7)==c(3) .or. c(7)==c(5)))) then
                if (c(17)>=2.5d0) then
                  c(11) = c(12)
                else
                  if (c(17)>1.5d0) then
                    c(9) = (c(8)-c(4))/(c(7)-c(3))
                    c(10) = (c(8)-c(6))/(c(7)-c(5))
                    diff1 = c(10)*diff1
                    diff = c(9)*c(6) - c(10)*c(4)
                  else
                    diff = c(6) - c(4)
                  end if
                  if (diff1<0.d0) then
                    diff1 = -diff1
                    diff = -diff
                  end if
                  if (diff1>c(14) .and. diff1>diff*tol) then
                    if (diff1>=c(12)*diff) then
                      c(11) = c(12)
                    else
                      c(11) = diff1/diff
                    end if
                  else
                    c(11) = tol
                  end if
                end if
                c(7) = c(5)
                c(8) = c(6)
                c(5) = c(3)
                c(6) = c(4)
                c(3) = c(3) + c(11)
                x = c(3)
                y = c(1)
                ind = 4
                return
              end if
            end if
            ind = 0
            i = 5
          end if
        else
          ind = 0
          i = 1
        end if
      else
        y = x
        ind = 0
        i = 0
      end if
    else
      ind = 0
      i = 3
    end if
  else
!   user not checked ind or changed it
    ind = 0      
    i = 2
  end if
  ifail = err_handle(ifail,i,srname,0,dumrec)
end subroutine zeroin_rev
!-------------------------------------------------------
subroutine zeroin(a,b,eps,eta,f,x,ifail) ! C05ADF
! driver for ZEROIN_REV
  use kinds, only : r8
  implicit none
! arguments
  real(kind=r8), intent(in)    :: a, b, eps, eta
  real(kind=r8), intent(out)   :: x
  integer      , intent(inout) :: ifail
! local variables
  integer       :: ifail1, ind, ir
  real(kind=r8) :: fx, y
  real(kind=r8), dimension(17) :: c
! parameters
  character(len=1), dimension(1), parameter :: dumrec=(/' '/)
  character(len=*), parameter :: srname='ZEROIN'
! intrinsic functions
!  real, intrinsic :: abs, sign
! external functions/subroutines
  interface
    function f(x)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: f
      real(kind=r8), intent(in) :: x
    end function f
    subroutine zeroin_rev(x,y,fx,tolx,ir,c,ind,ifail)
      use kinds, only : r8
      implicit none
      real(kind=r8), intent(in)    :: tolx
      real(kind=r8), intent(inout) :: fx, x, y, c(17)
      integer, intent(in)          :: ir
      integer, intent(inout)       :: ifail, ind
    end subroutine zeroin_rev
    integer function err_handle(ifail,ierror,srname,nrec,rec)
      implicit none
      integer, intent(in)          :: ierror, ifail, nrec
      character(len=*), intent(in) :: srname, rec(:)
    end function err_handle
  end interface
! executable statements
  ifail1 = 1
! NO INPUT ERROR
  if (a/=b .and. eps>0.d0) then
    x = a
    fx = f(x)
    ifail1 = 0
! NO ZERO AT INITIAL POINT
    if (abs(fx)>=eta .and. fx/=0.d0) then
      y = x
      c(1) = fx
      x = b
      fx = f(x)
! NO ZERO AT INITIAL POINT
      if (abs(fx)>=eta .and. fx/=0.d0) then
        ifail1 = 1
! ROOT IN RANGE
        if (sign(1.d0,fx)/=sign(1.d0,c(1))) then
          ir = 1
          ind = -1
          do
            call zeroin_rev(x,y,fx,eps,ir,c,ind,ifail1)
            if (ind==0) then
! IMPOSSIBLE EXIT
              if (ifail1>=1.and.ifail1<=3) then
                ifail1 = 4
! PROBABLY A POLE
              else if (ifail1==4) then
                ifail1 = 3
! TOO MUCH ACCURACY REQUESTED
              else if (ifail1==5) then
                ifail1 = 2
              end if
              exit
            end if
            fx = f(x)
! ZERO HIT EXACTLY
            if (abs(fx)<eta .or. fx==0.d0) then
              ifail1 = 0
              exit
            end if
          end do
        end if
      end if
    end if
  end if
  ifail = err_handle(ifail,ifail1,srname,0,dumrec)
end subroutine zeroin
!-------------------------------------------------------
!-------------------------------------------------------
! hybrj1 from MINPACK
!
!-------------------------------------------------------
function enorm(n,x)
! **********
!
! function enorm
!
! given an n-vector x, this function calculates the
! euclidean norm of x.
!
! the euclidean norm is computed by accumulating the sum of
! squares in three different sums. the sums of squares for the
! small and large components are scaled so that no overflows
! occur. non-destructive underflows are permitted. underflows
! and overflows do not occur in the computation of the unscaled
! sum of squares for the intermediate components.
! the definitions of small, intermediate and large components
! depend on two constants, rdwarf and rgiant. the main
! restrictions on these constants are that rdwarf**2 not
! underflow and rgiant**2 not overflow. the constants
! given here are suitable for every known computer.
!
! the function statement is
!
!   double precision function enorm(n,x)
!
! where
!
!   n is a positive integer input variable.
!
!   x is an input array of length n.
!
! subprograms called
!
!   fortran-supplied ... abs,sqrt
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  implicit none
  real(kind=r8) enorm
! arguments
  integer      , intent(in)  :: n
  real(kind=r8), intent(in) :: x(n)
! local variables
  integer       :: i
  real(kind=r8) :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
! parameters
  real(kind=r8), parameter :: one = 1.0d0, zero = 0.0d0, &
                              rdwarf = 3.834d-20, rgiant = 1.304d19
! intrinsic functions
!  real, intrinsic :: abs, sqrt
! executable statements
  s1 = zero
  s2 = zero
  s3 = zero
  x1max = zero
  x3max = zero
  floatn = n
  agiant = rgiant/floatn
  do i = 1, n
    xabs = abs(x(i))
    if (xabs <= rdwarf) then
!
!     sum for small components.
!
      if (xabs <= x3max) then
        if (xabs /= zero) s3 = s3 + (xabs/x3max)**2
      else
        s3 = one + s3*(x3max/xabs)**2
        x3max = xabs
      end if
    else if (xabs >= agiant) then
!
!     sum for large components.
!
      if (xabs <= x1max) then
        s1 = s1 + (xabs/x1max)**2
      else
        s1 = one + s1*(x1max/xabs)**2
        x1max = xabs
      end if
    else
!
!     sum for intermediate components.
!
      s2 = s2 + xabs**2
    end if
  end do
!
! calculation of norm.
!
  if (s1 == zero) then
    if (s2 == zero) then
      enorm = x3max*sqrt(s3)
    else if (s2 >= x3max) then
      enorm = sqrt(s2*(one+(x3max/s2)*(x3max*s3)))
    else
      enorm = sqrt(x3max*((s2/x3max)+(x3max*s3)))
    end if
  else
    enorm = x1max*sqrt(s1+(s2/x1max)/x1max)
  end if
end function enorm
!-------------------------------------------------------
subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
! **********
!
! subroutine dogleg
!
! given an m by n matrix a, an n by n nonsingular diagonal
! matrix d, an m-vector b, and a positive number delta, the
! problem is to determine the convex combination x of the
! gauss-newton and scaled gradient directions that minimizes
! (a*x - b) in the least squares sense, subject to the
! restriction that the euclidean norm of d*x be at most delta.
!
! this subroutine completes the solution of the problem
! if it is provided with the necessary information from the
! qr factorization of a. that is, if a = q*r, where q has
! orthogonal columns and r is an upper triangular matrix,
! then dogleg expects the full upper triangle of r and
! the first n components of (q transpose)*b.
!
! the subroutine statement is
!
! subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
!
! where
!
!   n is a positive integer input variable set to the order of r.
!
!   r is an input array of length lr which must contain the upper
!     triangular matrix r stored by rows.
!
!   lr is a positive integer input variable not less than
!     (n*(n+1))/2.
!
!   diag is an input array of length n which must contain the
!     diagonal elements of the matrix d.
!
!   qtb is an input array of length n which must contain the first
!     n elements of the vector (q transpose)*b.
!
!   delta is a positive input variable which specifies an upper
!     bound on the euclidean norm of d*x.
!
!   x is an output array of length n which contains the desired
!     convex combination of the gauss-newton direction and the
!     scaled gradient direction.
!
!   wa1 and wa2 are work arrays of length n.
!
! subprograms called
!
!   minpack-supplied ... enorm
!
!   fortran-supplied ... abs,max,min,sqrt
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
! arguments
  integer, intent(in) :: n, lr
  real(kind=r8), intent(in)  :: delta, r(lr), diag(n), qtb(n)
  real(kind=r8), intent(out) :: x(n), wa1(n), wa2(n)
! local variables
  integer       :: i, j, jj, jp1, k, l
  real(kind=r8) :: alpha, bnorm, gnorm, qnorm, sgnorm, sum, temp
! parameters
  real(kind=r8), parameter :: one=1.0d0, zero=0.0d0
! intrinsic functions
!  real, intrinsic :: abs, max, min, sqrt
! external functions
  interface
    function enorm(n,x)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: enorm
      integer      , intent(in) :: n
      real(kind=r8), intent(in) :: x(n)
    end function enorm
  end interface
! executable statements
!
! first, calculate the gauss-newton direction.
!
  jj = (n*(n + 1))/2 + 1
  do k = 1, n
    j = n - k + 1
    jp1 = j + 1
    jj = jj - k
    l = jj + 1
    sum = zero
    do i = jp1, n
      sum = sum + r(l)*x(i)
      l = l + 1
    end do
    temp = r(jj)
    if (temp == zero) then
      l = j
      do i = 1, j
        temp = max(temp,abs(r(l)))
        l = l + n - i
      end do
      temp = eps_r8*temp
      if (temp == zero) temp = eps_r8
    end if
    x(j) = (qtb(j) - sum)/temp
  end do
!
! test whether the gauss-newton direction is acceptable.
!
  do j = 1, n
    wa1(j) = zero
    wa2(j) = diag(j)*x(j)
  end do
  qnorm = enorm(n,wa2)
  if (qnorm <= delta) return
!
! the gauss-newton direction is not acceptable.
! next, calculate the scaled gradient direction.
!
  l = 1
  do j = 1, n
    temp = qtb(j)
    do i = j, n
      wa1(i) = wa1(i) + r(l)*temp
      l = l + 1
    end do
    wa1(j) = wa1(j)/diag(j)
  end do
!
! calculate the norm of the scaled gradient and test for
! the special case in which the scaled gradient is zero.
!
  gnorm = enorm(n,wa1)
  sgnorm = zero
  alpha = delta/qnorm
  if (gnorm /= zero) then
!
!   calculate the point along the scaled gradient
!   at which the quadratic is minimized.
!
    do j = 1, n
      wa1(j) = (wa1(j)/gnorm)/diag(j)
    end do
    l = 1
    do j = 1, n
      sum = zero
      do i = j, n
        sum = sum + r(l)*wa1(i)
        l = l + 1
      end do
      wa2(j) = sum
    end do
    temp = enorm(n,wa2)
    sgnorm = (gnorm/temp)/temp
!
!   test whether the scaled gradient direction is acceptable.
!
    alpha = zero
    if (sgnorm < delta) then
!
!     the scaled gradient direction is not acceptable.
!     finally, calculate the point along the dogleg
!     at which the quadratic is minimized.
!
      bnorm = enorm(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2&
           + sqrt((temp-(delta/qnorm))**2&
                   +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
    end if
  end if
!
! form appropriate convex combination of the gauss-newton
! direction and the scaled gradient direction.
!
  temp = (one - alpha)*min(sgnorm,delta)
  do j = 1, n
    x(j) = temp*wa1(j) + alpha*x(j)
  end do
  continue
end subroutine dogleg
!-------------------------------------------------------
subroutine qform(m,n,q,ldq,wa)
! **********
!
! subroutine qform
!
! this subroutine proceeds from the computed qr factorization of
! an m by n matrix a to accumulate the m by m orthogonal matrix
! q from its factored form.
!
! the subroutine statement is
!
! subroutine qform(m,n,q,ldq,wa)
!
! where
!
!   m is a positive integer input variable set to the number
!     of rows of a and the order of q.
!
!   n is a positive integer input variable set to the number
!     of columns of a.
!
!   q is an m by m array. on input the full lower trapezoid in
!     the first min(m,n) columns of q contains the factored form.
!     on output q has been accumulated into a square matrix.
!
!   ldq is a positive integer input variable not less than m
!     which specifies the leading dimension of the array q.
!
!   wa is a work array of length m.
!
! subprograms called
!
!   fortran-supplied ... min
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in) :: m,n,ldq
  real(kind=r8), intent(out) :: wa(m)
  real(kind=r8), intent(inout) :: q(ldq,m)
! local variables
  integer       :: i, j, jm1, k, l, minmn, np1
  real(kind=r8) :: sum, temp
! parameters
  real(kind=r8), parameter :: one=1.0d0, zero=0.0d0
! intrinsic functions
!  integer, intrinsic :: min
! executable statements
!
! zero out upper triangle of q in the first min(m,n) columns.
!
  minmn = min(m,n)
  do j = 2, minmn
    jm1 = j - 1
    do i = 1, jm1
      q(i,j) = zero
    end do
  end do
!
! initialize remaining columns to those of the identity matrix.
!
  np1 = n + 1
  do j = np1, m
    do i = 1, m
      q(i,j) = zero
    end do
    q(j,j) = one
  end do
!
! accumulate q from its factored form.
!
  do l = 1, minmn
    k = minmn - l + 1
    do i = k, m
      wa(i) = q(i,k)
      q(i,k) = zero
    end do
    q(k,k) = one
    if (wa(k) /= zero) then
      do j = k, m
        sum = zero
        do i = k, m
          sum = sum + q(i,j)*wa(i)
        end do
        temp = sum/wa(k)
        do i = k, m
          q(i,j) = q(i,j) - temp*wa(i)
        end do
      end do
    end if
  end do
end subroutine qform
!-------------------------------------------------------
subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
! **********
!
! subroutine qrfac
!
! this subroutine uses householder transformations with column
! pivoting (optional) to compute a qr factorization of the
! m by n matrix a. that is, qrfac determines an orthogonal
! matrix q, a permutation matrix p, and an upper trapezoidal
! matrix r with diagonal elements of nonincreasing magnitude,
! such that a*p = q*r. the householder transformation for
! column k, k = 1,2,...,min(m,n), is of the form
!
!                       t
!       i - (1/u(k))*u*u
!
! where u has zeros in the first k-1 positions. the form of
! this transformation and the method of pivoting first
! appeared in the corresponding linpack subroutine.
!
! the subroutine statement is
!
! subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
! where
!
!   m is a positive integer input variable set to the number
!     of rows of a.
!
!   n is a positive integer input variable set to the number
!     of columns of a.
!
!   a is an m by n array. on input a contains the matrix for
!     which the qr factorization is to be computed. on output
!     the strict upper trapezoidal part of a contains the strict
!     upper trapezoidal part of r, and the lower trapezoidal
!     part of a contains a factored form of q (the non-trivial
!     elements of the u vectors described above).
!
!   lda is a positive integer input variable not less than m
!     which specifies the leading dimension of the array a.
!
!   pivot is a logical input variable. if pivot is set true,
!     then column pivoting is enforced. if pivot is set false,
!     then no column pivoting is done.
!
!   ipvt is an integer output array of length lipvt. ipvt
!     defines the permutation matrix p such that a*p = q*r.
!     column j of p is column ipvt(j) of the identity matrix.
!     if pivot is false, ipvt is not referenced.
!
!   lipvt is a positive integer input variable. if pivot is false,
!     then lipvt may be as small as 1. if pivot is true, then
!     lipvt must be at least n.
!
!   rdiag is an output array of length n which contains the
!     diagonal elements of r.
!
!   acnorm is an output array of length n which contains the
!     norms of the corresponding columns of the input matrix a.
!     if this information is not needed, then acnorm can coincide
!     with rdiag.
!
!   wa is a work array of length n. if pivot is false, then wa
!     can coincide with rdiag.
!
! subprograms called
!
!   minpack-supplied ... enorm
!
!   fortran-supplied ... max,sqrt,min
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
! arguments
  integer, intent(in)  :: m, n, lda, lipvt
  integer, intent(out) :: ipvt(lipvt)
  logical, intent(in)  :: pivot
  real(kind=r8), intent(out)   :: rdiag(n), acnorm(n), wa(n)
  real(kind=r8), intent(inout) :: a(lda,n)
! local variables
  integer       :: i, j, jp1, k, kmax, minmn
  real(kind=r8) :: ajnorm, sum, temp
! parameters
  real(kind=r8), parameter :: one=1.0d0, p05=5.0d-2, zero=0.0d0
! intrinsic functions
!  integer, intrinsic :: min
!  real   , intrinsic :: max, sqrt
! external functions
  interface
    function enorm(n,x)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: enorm
      integer      , intent(in) :: n
      real(kind=r8), intent(in) :: x(n)
    end function enorm
  end interface
! executable statements
!
! compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm(m,a(1,j))
    rdiag(j) = acnorm(j)
    wa(j) = rdiag(j)
    if (pivot) ipvt(j) = j
  end do
!
! reduce a to r with householder transformations.
!
  minmn = min(m,n)
  do j = 1, minmn
    if (pivot) then
!
!     bring the column of largest norm into the pivot position.
!
      kmax = j
      do k = j, n
        if (rdiag(k) > rdiag(kmax)) kmax = k
      end do
      if (kmax /= j) then
        do i = 1, m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
        end do
        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)
        k = ipvt(j)
        ipvt(j) = ipvt(kmax)
        ipvt(kmax) = k
      end if
    end if
!
!   compute the householder transformation to reduce the
!   j-th column of a to a multiple of the j-th unit vector.
!
    ajnorm = enorm(m-j+1,a(j,j))
    if (ajnorm /= zero) then
      if (a(j,j) < zero) ajnorm = -ajnorm
      do i = j, m
        a(i,j) = a(i,j)/ajnorm
      end do
      a(j,j) = a(j,j) + one
!
!     apply the transformation to the remaining columns
!     and update the norms.
!
      jp1 = j + 1
      if (n >= jp1) then
        do k = jp1, n
          sum = zero
          do i = j, m
            sum = sum + a(i,j)*a(i,k)
          end do
          temp = sum/a(j,j)
          do i = j, m
            a(i,k) = a(i,k) - temp*a(i,j)
          end do
          if (pivot .and. rdiag(k) /= zero) then
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*sqrt(max(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 <= eps_r8) then
              rdiag(k) = enorm(m-j,a(jp1,k))
              wa(k) = rdiag(k)
            end if
          end if
        end do
      end if
    end if
    rdiag(j) = -ajnorm
  end do
end subroutine qrfac
!-------------------------------------------------------
subroutine r1mpyq(m,n,a,lda,v,w)
! **********
!
! subroutine r1mpyq
!
! given an m by n matrix a, this subroutine computes a*q where
! q is the product of 2*(n - 1) transformations
!
!       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
! and gv(i), gw(i) are givens rotations in the (i,n) plane which
! eliminate elements in the i-th and n-th planes, respectively.
! q itself is not given, rather the information to recover the
! gv, gw rotations is supplied.
!
! the subroutine statement is
!
! subroutine r1mpyq(m,n,a,lda,v,w)
!
! where
!
!   m is a positive integer input variable set to the number
!     of rows of a.
!
!   n is a positive integer input variable set to the number
!     of columns of a.
!
!   a is an m by n array. on input a must contain the matrix
!     to be postmultiplied by the orthogonal matrix q
!     described above. on output a*q has replaced a.
!
!   lda is a positive integer input variable not less than m
!     which specifies the leading dimension of the array a.
!
!   v is an input array of length n. v(i) must contain the
!     information necessary to recover the givens rotation gv(i)
!     described above.
!
!   w is an input array of length n. w(i) must contain the
!     information necessary to recover the givens rotation gw(i)
!     described above.
!
! subroutines called
!
!   fortran-supplied ... abs,sqrt
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)          :: m, n, lda
  real(kind=r8), intent(in)    :: v(n),w(n)
  real(kind=r8), intent(inout) :: a(lda,n)
! local variables
  integer       :: i, j, nmj, nm1
  real(kind=r8) :: cos, sin, temp
! parameters
  real(kind=r8), parameter :: one=1.0d0
! intrinsic functions
!  real, intrinsic :: abs, sqrt
! executable statements
!
! apply the first set of givens rotations to a.
!
  nm1 = n - 1
  if (nm1 < 1) return
  do nmj = 1, nm1
    j = n - nmj
    if (abs(v(j)) > one) cos = one/v(j)
    if (abs(v(j)) > one) sin = sqrt(one-cos**2)
    if (abs(v(j)) <= one) sin = v(j)
    if (abs(v(j)) <= one) cos = sqrt(one-sin**2)
    do i = 1, m
      temp = cos*a(i,j) - sin*a(i,n)
      a(i,n) = sin*a(i,j) + cos*a(i,n)
      a(i,j) = temp
    end do
  end do
!
! apply the second set of givens rotations to a.
!
  do j = 1, nm1
    if (abs(w(j)) > one) cos = one/w(j)
    if (abs(w(j)) > one) sin = sqrt(one-cos**2)
    if (abs(w(j)) <= one) sin = w(j)
    if (abs(w(j)) <= one) cos = sqrt(one-sin**2)
    do i = 1, m
      temp = cos*a(i,j) + sin*a(i,n)
      a(i,n) = -sin*a(i,j) + cos*a(i,n)
      a(i,j) = temp
    end do
  end do
end subroutine r1mpyq
!-------------------------------------------------------
subroutine r1updt(m,n,s,ls,u,v,w,sing)
! **********
!
! subroutine r1updt
!
! given an m by n lower trapezoidal matrix s, an m-vector u,
! and an n-vector v, the problem is to determine an
! orthogonal matrix q such that
!
!               t
!       (s + u*v )*q
!
! is again lower trapezoidal.
!
! this subroutine determines q as the product of 2*(n - 1)
! transformations
!
!       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
!
! where gv(i), gw(i) are givens rotations in the (i,n) plane
! which eliminate elements in the i-th and n-th planes,
! respectively. q itself is not accumulated, rather the
! information to recover the gv, gw rotations is returned.
!
! the subroutine statement is
!
! subroutine r1updt(m,n,s,ls,u,v,w,sing)
!
! where
!
!   m is a positive integer input variable set to the number
!     of rows of s.
!
!   n is a positive integer input variable set to the number
!     of columns of s. n must not exceed m.
!
!   s is an array of length ls. on input s must contain the lower
!     trapezoidal matrix s stored by columns. on output s contains
!     the lower trapezoidal matrix produced as described above.
!
!   ls is a positive integer input variable not less than
!     (n*(2*m-n+1))/2.
!
!   u is an input array of length m which must contain the
!     vector u.
!
!   v is an array of length n. on input v must contain the vector
!     v. on output v(i) contains the information necessary to
!     recover the givens rotation gv(i) described above.
!
!   w is an output array of length m. w(i) contains information
!     necessary to recover the givens rotation gw(i) described
!     above.
!
!   sing is a logical output variable. sing is set true if any
!     of the diagonal elements of the output s are zero. otherwise
!     sing is set false.
!
! subprograms called
!
!   fortran-supplied ... abs,sqrt
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more,
! john l. nazareth
!
! **********
  use kinds, only : r8
  use constants, only : huge_r8
  implicit none
! arguments
  integer, intent(in)  :: m, n, ls
  logical, intent(out) :: sing
  real(kind=r8), intent(in)    :: u(m)
  real(kind=r8), intent(out)   :: w(m)
  real(kind=r8), intent(inout) :: s(ls), v(n)
! local variables
  integer       :: i, j, jj, l, nmj, nm1
  real(kind=r8) :: cos, cotan, sin, tan, tau, temp
! parameters
  real(kind=r8), parameter :: one=1.0d0, p5=5.0d-1, p25=2.5d-1, zero=0.0d0
! intrinsic functions
!  real, intrinsic :: abs, sqrt
! executable statements
!
! initialize the diagonal element pointer.
!
  jj = (n*(2*m - n + 1))/2 - (m - n)
!
! move the nontrivial part of the last column of s into w.
!
  l = jj
  do i = n, m
    w(i) = s(l)
    l = l + 1
  end do
!
! rotate the vector v into a multiple of the n-th unit vector
! in such a way that a spike is introduced into w.
!
  nm1 = n - 1
  do nmj = 1, nm1
    j = n - nmj
    jj = jj - (m - j + 1)
    w(j) = zero
    if (v(j) /= zero) then
!
!   determine a givens rotation which eliminates the
!   j-th element of v.
!
      if (abs(v(n)) >= abs(v(j))) then
        tan = v(j)/v(n)
        cos = p5/sqrt(p25+p25*tan**2)
        sin = cos*tan
        tau = sin
      else
        cotan = v(n)/v(j)
        sin = p5/sqrt(p25+p25*cotan**2)
        cos = sin*cotan
        tau = one
        if (abs(cos)*huge_r8 > one) tau = one/cos
      end if
!
!    apply the transformation to v and store the information
!    necessary to recover the givens rotation.
!
      v(n) = sin*v(j) + cos*v(n)
      v(j) = tau
!
!    apply the transformation to s and extend the spike in w.
!
      l = jj
      do i = j, m
        temp = cos*s(l) - sin*w(i)
        w(i) = sin*s(l) + cos*w(i)
        s(l) = temp
        l = l + 1
      end do
    end if
  end do
!
! add the spike from the rank 1 update to w.
!
  do i = 1, m
    w(i) = w(i) + v(n)*u(i)
  end do
!
! eliminate the spike.
!
  sing = .false.
  do j = 1, nm1
    if (w(j) /= zero) then
!
!     determine a givens rotation which eliminates the
!     j-th element of the spike.
!
      if (abs(s(jj)) >= abs(w(j))) then
        tan = w(j)/s(jj)
        cos = p5/sqrt(p25+p25*tan**2)
        sin = cos*tan
        tau = sin
      else
        cotan = s(jj)/w(j)
        sin = p5/sqrt(p25+p25*cotan**2)
        cos = sin*cotan
        tau = one
        if (abs(cos)*huge_r8 > one) tau = one/cos
      end if
!
!     apply the transformation to s and reduce the spike in w.
!
      l = jj
      do i = j, m
        temp = cos*s(l) + sin*w(i)
        w(i) = -sin*s(l) + cos*w(i)
        s(l) = temp
        l = l + 1
      end do
!
!     store the information necessary to recover the
!     givens rotation.
!
      w(j) = tau
    end if
!
!    test for zero diagonal elements in the output s.
!
    if (s(jj) == zero) sing = .true.
    jj = jj + (m - j + 1)
  end do
!
! move w back into the last column of the output s.
!
  l = jj
  do i = n, m
    s(l) = w(i)
    l = l + 1
  end do
  if (s(jj) == zero) sing = .true.
end subroutine r1updt
!-------------------------------------------------------
subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,factor, &
                   nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,wa3,wa4)
! **********
!
! subroutine hybrj
!
! the purpose of hybrj is to find a zero of a system of
! n nonlinear functions in n variables by a modification
! of the powell hybrid method. the user must provide a
! subroutine which calculates the functions and the jacobian.
!
! the subroutine statement is
!
! subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
!                    mode,factor,nprint,info,nfev,njev,r,lr,qtf,
!                    wa1,wa2,wa3,wa4)
!
! where
!
!   fcn is the name of the user-supplied subroutine which
!     calculates the functions and the jacobian. fcn must
!     be declared in an external statement in the user
!     calling program, and should be written as follows.
!
!   subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!     integer n,ldfjac,iflag
!     double precision x(n),fvec(n),fjac(ldfjac,n)
!     ----------
!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.
!     ---------
!     return
!     end
!
!     the value of iflag should not be changed by fcn unless
!     the user wants to terminate execution of hybrj.
!     in this case set iflag to a negative integer.
!
!   n is a positive integer input variable set to the number
!     of functions and variables.
!
!   x is an array of length n. on input x must contain
!     an initial estimate of the solution vector. on output x
!     contains the final estimate of the solution vector.
!
!   fvec is an output array of length n which contains
!     the functions evaluated at the output x.
!
!   fjac is an output n by n array which contains the
!     orthogonal matrix q produced by the qr factorization
!     of the final approximate jacobian.
!
!   ldfjac is a positive integer input variable not less than n
!     which specifies the leading dimension of the array fjac.
!
!   xtol is a nonnegative input variable. termination
!     occurs when the relative error between two consecutive
!     iterates is at most xtol.
!
!   maxfev is a positive integer input variable. termination
!     occurs when the number of calls to fcn with iflag = 1
!     has reached maxfev.
!
!   diag is an array of length n. if mode = 1 (see
!     below), diag is internally set. if mode = 2, diag
!     must contain positive entries that serve as
!     multiplicative scale factors for the variables.
!
!   mode is an integer input variable. if mode = 1, the
!     variables will be scaled internally. if mode = 2,
!     the scaling is specified by the input diag. other
!     values of mode are equivalent to mode = 1.
!
!   factor is a positive input variable used in determining the
!     initial step bound. this bound is set to the product of
!     factor and the euclidean norm of diag*x if nonzero, or else
!     to factor itself. in most cases factor should lie in the
!     interval (.1,100.). 100. is a generally recommended value.
!
!   nprint is an integer input variable that enables controlled
!     printing of iterates if it is positive. in this case,
!     fcn is called with iflag = 0 at the beginning of the first
!     iteration and every nprint iterations thereafter and
!     immediately prior to return, with x and fvec available
!     for printing. fvec and fjac should not be altered.
!     if nprint is not positive, no special calls of fcn
!     with iflag = 0 are made.
!
!   info is an integer output variable. if the user has
!     terminated execution, info is set to the (negative)
!     value of iflag. see description of fcn. otherwise,
!     info is set as follows.
!
!     info = 0   improper input parameters.
!
!     info = 1   relative error between two consecutive iterates
!                is at most xtol.
!
!     info = 2   number of calls to fcn with iflag = 1 has
!                reached maxfev.
!
!     info = 3   xtol is too small. no further improvement in
!                the approximate solution x is possible.
!
!     info = 4   iteration is not making good progress, as
!                measured by the improvement from the last
!                five jacobian evaluations.
!
!     info = 5   iteration is not making good progress, as
!                measured by the improvement from the last
!                ten iterations.
!
!   nfev is an integer output variable set to the number of
!     calls to fcn with iflag = 1.
!
!   njev is an integer output variable set to the number of
!     calls to fcn with iflag = 2.
!
!   r is an output array of length lr which contains the
!     upper triangular matrix produced by the qr factorization
!     of the final approximate jacobian, stored rowwise.
!
!   lr is a positive integer input variable not less than
!     (n*(n+1))/2.
!
!   qtf is an output array of length n which contains
!     the vector (q transpose)*fvec.
!
!   wa1, wa2, wa3, and wa4 are work arrays of length n.
!
! subprograms called
!
!   user-supplied ...... fcn
!
!   minpack-supplied ... dogleg,enorm,
!                        qform,qrfac,r1mpyq,r1updt
!
!   fortran-supplied ... abs,max,min,mod
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  use constants, only : eps_r8
  implicit none
! arguments
  integer, intent(in)  :: n, ldfjac, maxfev, mode, nprint, lr
  integer, intent(out) :: info, nfev, njev
  real(kind=r8), intent(in)    :: xtol, factor
  real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), r(lr), qtf(n), &
                                  wa1(n), wa2(n), wa3(n), wa4(n)
  real(kind=r8), intent(inout) :: x(n), diag(n)
! local variables
  integer       :: i, iflag, iter, j, jm1, l, ncfail, ncsuc, nslow1, nslow2
  integer, dimension(1) :: iwa
  logical       :: jeval, sing, term
  real(kind=r8) :: actred, delta, fnorm, fnorm1, pnorm, prered, &
          ratio, sum, temp, xnorm
! parameters
  real(kind=r8), parameter :: one = 1.0d0, zero = 0.0d0, p1 = 1.0d-1,   &
             p5 = 5.0d-1, p001 = 1.0d-3, p0001 = 1.0d-4
! intrinsic functions
!  real, intrinsic :: abs, max, min
!  integer, intrinsic :: mod
! external functions/subroutines
  interface
    subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: n, ldfjac
      integer, intent(inout) :: iflag
      real(kind=r8), intent(in)    :: x(n)
      real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fcn
    subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: n, lr
      real(kind=r8), intent(in)  :: delta, r(lr), diag(n), qtb(n)
      real(kind=r8), intent(out) :: x(n), wa1(n), wa2(n)
    end subroutine dogleg
    function enorm(n,x)
      use kinds, only : r8
      implicit none
      real(kind=r8) :: enorm
      integer      , intent(in)  :: n
      real(kind=r8), intent(in) :: x(n)
    end function enorm
    subroutine qform(m,n,q,ldq,wa)
      use kinds, only : r8
      implicit none
      integer, intent(in) :: m,n,ldq
      real(kind=r8), intent(out) :: wa(m)
      real(kind=r8), intent(inout) :: q(ldq,m)
    end subroutine qform
    subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: m, n, lda, lipvt
      integer, intent(out) :: ipvt(lipvt)
      logical, intent(in)  :: pivot
      real(kind=r8), intent(out)   :: rdiag(n), acnorm(n), wa(n)
      real(kind=r8), intent(inout) :: a(lda,n)
    end subroutine qrfac
    subroutine r1mpyq(m,n,a,lda,v,w)
      use kinds, only : r8
      implicit none
      integer, intent(in)          :: m, n, lda
      real(kind=r8), intent(in)    :: v(n),w(n)
      real(kind=r8), intent(inout) :: a(lda,n)
    end subroutine r1mpyq
    subroutine r1updt(m,n,s,ls,u,v,w,sing)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: m, n, ls
      logical, intent(out) :: sing
      real(kind=r8), intent(in)    :: u(m)
      real(kind=r8), intent(out)   :: w(m)
      real(kind=r8), intent(inout) :: s(ls), v(n)
    end subroutine r1updt
  end interface
! executable statements
  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
! check the input parameters for errors.
!
  if (n>0 .and. ldfjac>=n .and. xtol>=zero &
      .and. maxfev>0 .and. factor>zero      &
      .and. lr>=(n*(n+1))/2) then
    term=.false.
    if (mode==2) then
      do j = 1, n
        if (diag(j) <= zero) then
          term=.true.
          exit
        end if
      end do
    end if
    if (.not.term) then
!
! evaluate the function at the starting point
! and calculate its norm.
!
      iflag = 1
      call fcn(n,x,fvec,fjac,ldfjac,iflag)
      nfev = 1
      if (iflag>=0) then
        fnorm = enorm(n,fvec)
!
!       initialize iteration counter and monitors.
!
        iter = 1
        ncsuc = 0
        ncfail = 0
        nslow1 = 0
        nslow2 = 0
!
!       beginning of the outer loop.
!
        do
          jeval = .true.
!
!         calculate the jacobian matrix.
!
          iflag = 2
          call fcn(n,x,fvec,fjac,ldfjac,iflag)
          njev = njev + 1
          if (iflag < 0) exit
!
!         compute the qr factorization of the jacobian.
!
          call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!
!         on the first iteration and if mode is 1, scale according
!         to the norms of the columns of the initial jacobian.
!
          if (iter == 1) then
            if (mode /= 2) then
              do j = 1, n
                diag(j) = wa2(j)
                if (wa2(j) == zero) diag(j) = one
              end do
            end if
!
!         on the first iteration, calculate the norm of the scaled x
!         and initialize the step bound delta.
!
            do j = 1, n
              wa3(j) = diag(j)*x(j)
            end do
            xnorm = enorm(n,wa3)
            delta = factor*xnorm
            if (delta == zero) delta = factor
          end if
!
!         form (q transpose)*fvec and store in qtf.
!
          do i = 1, n
            qtf(i) = fvec(i)
          end do
          do j = 1, n
            if (fjac(j,j) /= zero) then
              sum = zero
              do i = j, n
                sum = sum + fjac(i,j)*qtf(i)
              end do
              temp = -sum/fjac(j,j)
              do i = j, n
                qtf(i) = qtf(i) + fjac(i,j)*temp
              end do
            end if
          end do
!
!         copy the triangular factor of the qr factorization into r.
!
          sing = .false.
          do j = 1, n
            l = j
            jm1 = j - 1
            do i = 1, jm1
              r(l) = fjac(i,j)
              l = l + n - i
            end do
            r(l) = wa1(j)
            if (wa1(j) == zero) sing = .true.
          end do
!
!         accumulate the orthogonal factor in fjac.
!
          call qform(n,n,fjac,ldfjac,wa1)
!
!         rescale if necessary.
!
          if (mode /= 2) then
            do j = 1, n
              diag(j) = max(diag(j),wa2(j))
            end do
          end if
!
!         beginning of the inner loop.
!
          do
!
!           if requested, call fcn to enable printing of iterates.
!
            if (nprint > 0) then
              iflag = 0
              if (mod(iter-1,nprint) == 0) &
                  call fcn(n,x,fvec,fjac,ldfjac,iflag)
              if (iflag < 0) then
                term=.true.
                exit
              end if
            end if
!
!           determine the direction p.
!
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!
!           store the direction p and x+p. calculate the norm of p.
!
            do j = 1, n
              wa1(j) = -wa1(j)
              wa2(j) = x(j) + wa1(j)
              wa3(j) = diag(j)*wa1(j)
            end do
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter == 1) delta = min(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
            nfev = nfev + 1
            if (iflag < 0) then
              term=.true.
              exit
            end if
            fnorm1 = enorm(n,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction.
!
            l = 1
            do i = 1, n
              sum = zero
              do j = i, n
                sum = sum + r(l)*wa1(j)
                l = l + 1
              end do
              wa3(i) = qtf(i) + sum
            end do
            temp = enorm(n,wa3)
            prered = zero
            if (temp < fnorm) prered = one - (temp/fnorm)**2
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered > zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio >= p1) then
              ncfail = 0
              ncsuc = ncsuc + 1
              if (ratio >= p5 .or. ncsuc > 1)&
                  delta = max(delta,pnorm/p5)
              if (abs(ratio-one) <= p1) delta = pnorm/p5
            else
              ncsuc = 0
              ncfail = ncfail + 1
              delta = p5*delta
            end if
!
!           test for successful iteration.
!
            if (ratio >= p0001) then
!
!             successful iteration. update x, fvec, and their norms.
!
              do j = 1, n
                x(j) = wa2(j)
                wa2(j) = diag(j)*x(j)
                fvec(j) = wa4(j)
              end do
              xnorm = enorm(n,wa2)
              fnorm = fnorm1
              iter = iter + 1
            end if
!
!           determine the progress of the iteration.
!
            nslow1 = nslow1 + 1
            if (actred >= p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred >= p1) nslow2 = 0
!
!           test for convergence.
!
            if (delta <= xtol*xnorm .or. fnorm == zero) info=1
            if (info /= 0) then
              term=.true.
              exit
            end if
!
!           tests for termination and stringent tolerances.
!
            if (nfev >= maxfev) info = 2
            if (p1*max(p1*delta,pnorm) <= eps_r8*xnorm) info=3
            if (nslow2 == 5) info = 4
            if (nslow1 == 10) info = 5
            if (info /= 0) then
              term=.true.
              exit
            end if
!
!           criterion for recalculating jacobian.
!
            if (ncfail == 2) exit
!
!           calculate the rank one modification to the jacobian
!           and update qtf if necessary.
!
            do j = 1, n
              sum = zero
              do i = 1, n
                sum = sum + fjac(i,j)*wa4(i)
              end do
              wa2(j) = (sum - wa3(j))/pnorm
              wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
              if (ratio >= p0001) qtf(j) = sum
            end do
!
!           compute the qr factorization of the updated jacobian.
!
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
!
!           end of the inner loop.
!
            jeval = .false.
          end do
          if (term) exit
!
!   end of the outer loop.
!
        end do
      end if
    end if
  end if
!
! termination, either normal or user imposed.
!
  if (iflag < 0) info = iflag
  iflag = 0
  if (nprint > 0) call fcn(n,x,fvec,fjac,ldfjac,iflag)
end subroutine hybrj
!-------------------------------------------------------
subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
! **********
!
! subroutine hybrj1
!
! the purpose of hybrj1 is to find a zero of a system of
! n nonlinear functions in n variables by a modification
! of the powell hybrid method. this is done by using the
! more general nonlinear equation solver hybrj. the user
! must provide a subroutine which calculates the functions
! and the jacobian.
!
! the subroutine statement is
!
! subroutine hybrj1(fcn,n,x,fvec,fjac,ldfjac,tol,info,wa,lwa)
!
! where
!
!   fcn is the name of the user-supplied subroutine which
!     calculates the functions and the jacobian. fcn must
!     be declared in an external statement in the user
!     calling program, and should be written as follows.
!
!   subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
!     integer n,ldfjac,iflag
!     double precision x(n),fvec(n),fjac(ldfjac,n)
!     ----------
!     if iflag = 1 calculate the functions at x and
!     return this vector in fvec. do not alter fjac.
!     if iflag = 2 calculate the jacobian at x and
!     return this matrix in fjac. do not alter fvec.
!     ---------
!     return
!     end
!
!     the value of iflag should not be changed by fcn unless
!     the user wants to terminate execution of hybrj1.
!     in this case set iflag to a negative integer.
!
!   n is a positive integer input variable set to the number
!     of functions and variables.
!
!   x is an array of length n. on input x must contain
!     an initial estimate of the solution vector. on output x
!     contains the final estimate of the solution vector.
!
!   fvec is an output array of length n which contains
!     the functions evaluated at the output x.
!
!   fjac is an output n by n array which contains the
!     orthogonal matrix q produced by the qr factorization
!     of the final approximate jacobian.
!
!   ldfjac is a positive integer input variable not less than n
!     which specifies the leading dimension of the array fjac.
!
!   tol is a nonnegative input variable. termination occurs
!     when the algorithm estimates that the relative error
!     between x and the solution is at most tol.
!
!   info is an integer output variable. if the user has
!     terminated execution, info is set to the (negative)
!     value of iflag. see description of fcn. otherwise,
!     info is set as follows.
!
!     info = 0   improper input parameters.
!
!     info = 1   algorithm estimates that the relative error
!                between x and the solution is at most tol.
!
!     info = 2   number of calls to fcn with iflag = 1 has
!                reached 100*(n+1).
!
!     info = 3   tol is too small. no further improvement in
!                the approximate solution x is possible.
!
!     info = 4   iteration is not making good progress.
!
!   wa is a work array of length lwa.
!
!   lwa is a positive integer input variable not less than
!     (n*(n+13))/2.
!
! subprograms called
!
!   user-supplied ...... fcn
!
!   minpack-supplied ... hybrj
!
! argonne national laboratory. minpack project. march 1980.
! burton s. garbow, kenneth e. hillstrom, jorge j. more
!
! **********
  use kinds, only : r8
  implicit none
! arguments
  integer, intent(in)          :: n, ldfjac, lwa
  integer, intent(out)         :: info
  real(kind=r8), intent(in)    :: tol
  real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), wa(lwa)
  real(kind=r8), intent(inout) :: x(n)
! local variables
  integer       :: j, lr, maxfev, mode, nfev, njev, nprint
  integer       :: n2,n3,n4,n5,n6
  real(kind=r8) :: xtol
! parameters
  real(kind=r8), parameter :: factor=1.0d2, one=1.0d0, zero=0.0d0
  interface
    subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
      use kinds, only : r8
      implicit none
      integer, intent(in)    :: n, ldfjac
      integer, intent(inout) :: iflag
      real(kind=r8), intent(in)    :: x(n)
      real(kind=r8), intent(inout) :: fvec(n), fjac(ldfjac,n)
    end subroutine fcn
    subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,factor, &
                       nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,wa3,wa4)
      use kinds, only : r8
      implicit none
      integer, intent(in)  :: n, ldfjac, maxfev, mode, nprint, lr
      integer, intent(out) :: info, nfev, njev
      real(kind=r8), intent(in)    :: xtol, factor
      real(kind=r8), intent(out)   :: fvec(n), fjac(ldfjac,n), r(lr), qtf(n), &
                                      wa1(n), wa2(n), wa3(n), wa4(n)
      real(kind=r8), intent(inout) :: x(n), diag(n)
      external fcn
    end subroutine hybrj
  end interface
! executable statements
  info = 0
!
! check the input parameters for errors.
!
  if (n>0 .and. ldfjac>=n .and. tol>=zero &
      .and. lwa>=(n*(n+13))/2) then
!
! call hybrj.
!
    n2=2*n; n3=3*n; n4=4*n; n5=5*n; n6=6*n
    maxfev = 100*(n + 1)
    xtol = tol
    mode = 2
    do j = 1, n
      wa(j) = one
    end do
    nprint = 0
    lr = (n*(n + 1))/2
    call hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,wa(1:n),mode,    &
               factor,nprint,info,nfev,njev,wa(n6+1:),lr,wa(n+1:n2),  &
               wa(n2+1:n3),wa(n3+1:n4),wa(n4+1:n5),wa(n5+1:n6))
    if (info == 5) info = 4
  end if
end subroutine hybrj1
!-------------------------------------------------------
subroutine fpcuro(a,b,c,d,x,n)
	use kinds, only : r8
	implicit none
!  subroutine fpcuro finds the real zeros of a cubic polynomial
!  p(x) = a*x**3+b*x**2+c*x+d.
!
!  calling sequence:
!     call fpcuro(a,b,c,d,x,n)
!
!  input parameters:
!    a,b,c,d: real values, containing the coefficients of p(x).
!
!  output parameters:
!    x      : real array,length 3, which contains the real zeros of p(x)
!    n      : integer, giving the number of real zeros of p(x).
!  ..
!  ..scalar arguments..
	real(kind=r8), intent(in) :: a,b,c,d
	integer, intent(out) :: n
!  ..array argument..
	real(kind=r8), intent(out) :: x(3)
!  ..local scalars..
	integer i
	real(kind=r8) a1,b1,c1,df,disc,d1,e3,f,four,half,ovfl,pi3,p3,q,r, &
		step,tent,three,two,u,u1,u2,y
!  ..function references..
	real(kind=r8) abs,max,atan,atan2,cos,sign,sqrt
!  set constants
	two = 0.2d+01
	three = 0.3d+01
	four = 0.4d+01
	ovfl =0.1d+05
	half = 0.5d+0
	tent = 0.1d+0
	e3 = tent/0.3d0
	pi3 = atan(0.1d+01)/0.75d0
	a1 = abs(a)
	b1 = abs(b)
	c1 = abs(c)
	d1 = abs(d)
	if(max(b1,c1,d1).lt.a1*ovfl) then
!  p(x) is a third degree polynomial.
		b1 = b/a*e3
		c1 = c/a
		d1 = d/a
		q = c1*e3-b1*b1
		r = b1*b1*b1+(d1-b1*c1)*half
		disc = q*q*q+r*r
		if(disc.gt.0.0d0) then
			u = sqrt(disc)
			u1 = -r+u
			u2 = -r-u
			n = 1
			x(1) = sign(abs(u1)**e3,u1)+sign(abs(u2)**e3,u2)-b1
		else
			u = sqrt(abs(q))
			if(r.lt.0.0d0) u = -u
			p3 = atan2(sqrt(-disc),abs(r))*e3
			u2 = u+u
			n = 3
			x(1) = -u2*cos(p3)-b1
			x(2) = u2*cos(pi3-p3)-b1
			x(3) = u2*cos(pi3+p3)-b1
		end if
	else if(max(c1,d1).lt.b1*ovfl) then
!  p(x) is a second degree polynomial.
		disc = c*c-four*b*d
		n = 0
		if(disc.lt.0.0d0) return
		n = 2
		u = sqrt(disc)
		b1 = b+b
		x(1) = (-c+u)/b1
		x(2) = (-c-u)/b1
	else if(d1.lt.c1*ovfl) then
!  p(x) is a first degree polynomial.
		n = 1
		x(1) = -d/c
	else
!  p(x) is a constant function.
		n = 0
		return
	end if
!  apply a newton iteration to improve the accuracy of the roots.
	do i=1,n
		y = x(i)
		f = ((a*y+b)*y+c)*y+d
		df = (three*a*y+two*b)*y+c
		step = 0.0d0
		if(abs(f).lt.abs(df)*tent) step = f/df
		x(i) = y-step
	end do
end subroutine fpcuro
!-------------------------------------------------------
subroutine sproota(val,t,n,c,zero,mest,m,ier)
	use kinds, only : r8
	implicit none
!  subroutine sproot finds the zeros of a cubic spline s(x),which is
!  given in its normalized b-spline representation.
!
!  calling sequence:
!     call sproot(t,n,c,zero,mest,m,ier)
!
!  input parameters:
!    t    : real array,length n, containing the knots of s(x).
!    n    : integer, containing the number of knots.  n>=8
!    c    : real array,length n, containing the b-spline coefficients.
!    mest : integer, specifying the dimension of array zero.
!
!  output parameters:
!    zero : real array,lenth mest, containing the zeros of s(x).
!    m    : integer,giving the number of zeros.
!    ier  : error flag:
!      ier = 0: normal return.
!      ier = 1: the number of zeros exceeds mest.
!      ier =10: invalid input data (see restrictions).
!
!  other subroutines required: fpcuro
!
!  restrictions:
!    1) n>= 8.
!    2) t(4) < t(5) < ... < t(n-4) < t(n-3).
!       t(1) <= t(2) <= t(3) <= t(4)
!       t(n-3) <= t(n-2) <= t(n-1) <= t(n)
!
!  author :
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be
!
!  latest update : march 1987
!
! ..
! ..scalar arguments..
	integer, intent(in) :: n,mest
	integer, intent(out) :: m,ier
	real(kind=r8), intent(in) :: val
!  ..array arguments..
	real(kind=r8), intent(in) :: t(n),c(n)
	real(kind=r8), intent(out) :: zero(mest)
!  ..local scalars..
	integer :: i,j,j1,l,n4
	real(kind=r8) :: ah,a0,a1,a2,a3,bh,b0,b1,c1,c2,c3,c4,c5,d4,d5,h1,h2, &
		three,two,t1,t2,t3,t4,t5,zz
	logical :: z0,z1,z2,z3,z4,nz0,nz1,nz2,nz3,nz4
!  ..local array..
	real(kind=r8) :: y(3)
!  ..
	interface
		subroutine fpcuro(a,b,c,d,x,n)
			use kinds, only : r8
			implicit none
			real(kind=r8), intent(in) :: a,b,c,d
			integer, intent(out) :: n
			real(kind=r8), intent(out) :: x(3)
		end subroutine fpcuro
	end interface
!  set some constants
	two = 0.2d+01
	three = 0.3d+01
!  before starting computations a data check is made. if the input data
!  are invalid, control is immediately repassed to the calling program.
	n4 = n-4
	ier = 10
	if(n<8) return
	j = n
	do i=1,3
		if(t(i) > t(i+1)) return
		if(t(j) < t(j-1)) return
		j = j-1
	end do
	do i=4,n4
		if(t(i) >= t(i+1)) return
	end do
!  the problem considered reduces to finding the zeros of the cubic
!  polynomials pl(x) which define the cubic spline in each knot
!  interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
!  the condition that it belongs to the knot interval.
!  the cubic polynomial pl(x) is determined by computing s(t(l)),
!  s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
!  s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
!  splines and their derivatives, the value of s(t(l)) and s'(t(l))
!  is already known from the foregoing knot interval.
	ier = 0
!  evaluate some constants for the first knot interval
	h1 = t(4)-t(3)
	h2 = t(5)-t(4)
	t1 = t(4)-t(2)
	t2 = t(5)-t(3)
	t3 = t(6)-t(4)
	t4 = t(5)-t(2)
	t5 = t(6)-t(3)
!  calculate a0 = s(t(4)) and ah = s'(t(4)).
	c1 = c(1)
	c2 = c(2)
	c3 = c(3)
	c4 = (c2-c1)/t4
	c5 = (c3-c2)/t5
	d4 = (h2*c1+t1*c2)/t4
	d5 = (t3*c2+h1*c3)/t5
	a0 = (h2*d4+h1*d5)/t2 - val
	ah = three*(h2*c4+h1*c5)/t2
	z1 = .true.
	if(ah.lt.0.0d0) z1 = .false.
	nz1 = .not.z1
	m = 0
!  main loop for the different knot intervals.
	do l=4,n4
!  evaluate some constants for the knot interval t(l) <= x <= t(l+1).
		h1 = h2
		h2 = t(l+2)-t(l+1)
		t1 = t2
		t2 = t3
		t3 = t(l+3)-t(l+1)
		t4 = t5
		t5 = t(l+3)-t(l)
!  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
		c1 = c2
		c2 = c3
		c3 = c(l)
		c4 = c5
		c5 = (c3-c2)/t5
		d4 = (h2*c1+t1*c2)/t4
		d5 = (h1*c3+t3*c2)/t5
		b0 = (h2*d4+h1*d5)/t2 - val
		bh = three*(h2*c4+h1*c5)/t2
!  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
!  pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
		a1 = ah*h1
		b1 = bh*h1
		a2 = three*(b0-a0)-b1-two*a1
		a3 = two*(a0-b0)+b1+a1
!  test whether or not pl(x) could have a zero in the range
!  t(l) <= x <= t(l+1).
		z3 = .true.
		if(b1 < 0.0d0) z3 = .false.
		nz3 = .not.z3
		if(a0*b0 > 0.0d0) then
			z0 = .true.
			if(a0 < 0.0d0) z0 = .false.
			nz0 = .not.z0
			z2 = .true.
			if(a2 < 0.0d0) z2 = .false.
			nz2 = .not.z2
			z4 = .true.
			if(3.0d0*a3+a2 < 0.0d0) z4 = .false.
			nz4 = .not.z4
		end if
		if( (a0*b0<=0.0d0) .or. &
				(z0.and.(nz1.and.(z3.or.z2.and.nz4).or.nz2.and.z3.and.z4).or. &
				nz0.and.(z1.and.(nz3.or.nz2.and.z4).or.z2.and.nz3.and.nz4))) then
!  find the zeros of ql(y).
			call fpcuro(a3,a2,a1,a0,y,j)
			if(j/=0) then
!  find which zeros of pl(x) are zeros of s(x).
				do i=1,j
					if(y(i) >= 0.0d0 .and. y(i) <= 1.0d0) then
!  test whether the number of zeros of s(x) exceeds mest.
						if(m >= mest) then
							ier = 1
							return
						end if
						m = m+1
						zero(m) = t(l)+h1*y(i)
					end if
	 			end do
 			end if
 		end if
		a0 = b0
		ah = bh
		z1 = z3
		nz1 = nz3
	end do
!  the zeros of s(x) are arranged in increasing order.
	if(m < 2) return
	do i=2,m
		j = i
		do
			j1 = j-1
			if(j1==0) exit
			if(zero(j) >= zero(j1)) exit
			zz = zero(j)
			zero(j) = zero(j1)
			zero(j1) = zz
			j = j1
		end do
	end do
	j = m
	m = 1
	do i=2,j
		if(zero(i) /= zero(m)) then
			m = m+1
			zero(m) = zero(i)
		end if
	end do
end subroutine sproota
