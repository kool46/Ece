module cubic
  implicit none
contains
!-----------------------------------------------------------------------
  subroutine aord(a, n)
!-----------------------------------------------------------------------
! THE AORD SORTING PROCEDURE IS USED TO REORDER THE ELEMENTS OF A SO THAT
! ABS(A(I)) .LE. ABS(A(I+1)) FOR I = 1,...,N-1.  IT IS ASSUMED THAT N > 1.
!-----------------------------------------------------------------------
    use kinds, only : r8
    real(r8), intent(inout) :: a(:)
    integer, intent(in)     :: n
! local variables
    integer  :: k(10) = (/ 1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 28524 /), &
                imax, i, ii, ki, jmax, j, l, ll
    real(r8) :: s
!             SELECTION OF THE INCREMENTS K(I) = (3**I-1)/2
    if (n > 1) then
      imax = 1
      do i = 3,10
        if (n <= k(i)) exit
        imax = imax + 1
      end do
!            STEPPING THROUGH THE INCREMENTS K(IMAX),...,K(1)
      i = imax
      do ii = 1,imax
        ki = k(i)
!             SORTING ELEMENTS THAT ARE KI POSITIONS APART
!                 SO THAT ABS(A(J)) .LE. ABS(A(J+KI))
        jmax = n - ki
        do j = 1,jmax
          l = j
          ll = j + ki
          s = a(ll)
          do
            if (abs(s) >= abs(a(l))) exit
            a(ll) = a(l)
            ll = l
            l = l - ki
            if (l <= 0) exit
          end do
          a(ll) = s
        end do
        i = i - 1
      end do
    end if
  end subroutine aord
!-----------------------------------------------------------------------
  subroutine ord(a, n)
!-----------------------------------------------------------------------
! THE ORD SORTING PROCEDURE IS USED TO REORDER THE ELEMENTS OF A SO THAT
! A(I) .LE. A(I+1) FOR I = 1,...,N-1.  IT IS ASSUMED THAT N > 1.
!-----------------------------------------------------------------------
    use kinds, only : r8
    real(r8), intent(inout) :: a(:)
    integer, intent(in)     :: n
! local variables
    integer  :: k(10) = (/ 1, 4, 13, 40, 121, 364, 1093, 3280, 9841, 28524 /), &
                imax, i, ii, ki, jmax, j, l, ll
    real(r8) :: s
!             SELECTION OF THE INCREMENTS K(I) = (3**I-1)/2
    if (n > 1) then
      imax = 1
      do i = 3,10
        if (n <= k(i)) exit
        imax = imax + 1
      end do
!            STEPPING THROUGH THE INCREMENTS K(IMAX),...,K(1)
      i = imax
      do ii = 1,imax
        ki = k(i)
!             SORTING ELEMENTS THAT ARE KI POSITIONS APART
!                 SO THAT ABS(A(J)) .LE. ABS(A(J+KI))
        jmax = n - ki
        do j = 1,jmax
          l = j
          ll = j + ki
          s = a(ll)
          do
            if (s >= a(l)) exit
            a(ll) = a(l)
            ll = l
            l = l - ki
            if (l <= 0) exit
          end do
          a(ll) = s
        end do
        i = i - 1
      end do
    end if
  end subroutine ord
!-----------------------------------------------------------------------
  function cbrt(x) result(fn_val)
!-----------------------------------------------------------------------
!                   CUBE ROOT OF A REAL NUMBER
!-----------------------------------------------------------------------
    use kinds, only : r8
    use constants, only : zero
    real(r8), intent(in) :: x
    real(r8)             :: fn_val
!     local variable
    real(r8) :: r
    real(r8), parameter :: three = 3.0_r8
    if (x < zero) then
      r = log(-x) / three
      fn_val = -exp(r)
    else if (x > zero) then
      r = log(x) / three
      fn_val = exp(r)
    else
      fn_val = zero
    end if
  end function cbrt
!-----------------------------------------------------------------------
  subroutine qdcrt (a, z, nreal, re)
!-----------------------------------------------------------------------
!        QDCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2
!     AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(3) IS NONZERO.
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
    use kinds, only : r8
    use constants, only : zero, half, two, eps_r8
    real(r8), intent(in)     :: a(:)
    complex(r8), intent(out) :: z(:)
    integer,  optional, intent(out) :: nreal
    real(r8), optional, intent(out) :: re(:)
! local variables
    real(r8)            :: d, r, w, x, y
    real(r8), parameter :: four = 4.0_r8
    if (a(1) == zero) then
!                 CASE WHEN A(1) = 0
      z(1) = cmplx(zero, zero)
      z(2) = cmplx(-a(2)/a(3), zero)
      if (present(nreal)) nreal = 2
      if (present(re)) then
        re(:2) = z(:2)
        call ord(re,2)
      end if
    else
      d = a(2)*a(2) - four*a(1)*a(3)
      if (abs(d) <= two*eps_r8*a(2)*a(2)) then
!                  EQUAL REAL ROOTS
        z(1) = cmplx(-half*a(2)/a(3), zero)
        z(2) = z(1)
        if (present(nreal)) nreal = 2
        if (present(re)) then
          re(:2) = z(:2)
          call ord(re,2)
        end if
      else
        r = sqrt(abs(d))
        if (d < zero) then
!                   COMPLEX ROOTS
          x = -half*a(2)/a(3)
          y = abs(half*r/a(3))
          z(1) = cmplx(x, y)
          z(2) = cmplx(x,-y)
          if (present(nreal)) nreal = 0
        else
!                 DISTINCT REAL ROOTS
          if (a(2) == zero) then
            x = abs(half*r/a(3))
            z(1) = cmplx(x, zero)
            z(2) = cmplx(-x, zero)
          else
            w = -(a(2) + sign(r,a(2)))
            z(1) = cmplx(two*a(1)/w, zero)
            z(2) = cmplx(half*w/a(3), zero)
          end if
          if (present(nreal)) nreal = 2
          if (present(re)) then
            re(:2) = z(:2)
            call ord(re,2)
          end if
        end if
      end if
    end if
  end subroutine qdcrt
!-----------------------------------------------------------------------
  subroutine cbcrt(a, z, nreal, re)
!-----------------------------------------------------------------------
!        CBCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!              A(1) + A(2)*Z + A(3)*Z**2 + A(4)*Z**3
!     AND STORES THE RESULTS IN Z. IT IS ASSUMED THAT A(4) IS NONZERO.
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
    use kinds, only : r8
    use constants, only : zero, half, one, two, eps_r8
    real(r8), intent(in)     :: a(:)
    complex(r8), intent(out) :: z(:)
    integer,  optional, intent(out) :: nreal
    real(r8), optional, intent(out) :: re(:)
    real(r8) :: aq(3), arg, c, cf, d, p, p1, q, q1, r, ra, rb, rq, rt, &
                r1, s, sf, sq, sum, t, tol, t1, w, &
                w1, w2, x, x1, x2, x3, y, y1, y2, y3
    real(r8), parameter :: three = 3.0_r8, rt3 = 1.7320508075688772935274463415059_r8
    logical  :: bigw
    if (a(1) == zero) then
!                  CASE WHEN A(1) = 0
      z(1) = (zero, zero)
      call qdcrt(a(2:), z(2:), nreal, re(2:))
      if (present(nreal)) nreal = nreal + 1
      if (present(re)) then
        re(1) = z(1)
        call ord(re,nreal)
      end if
    else
      p = a(3)/(three*a(4))
      q = a(2)/a(4)
      r = a(1)/a(4)
      tol = 4._r8*eps_r8
      c = zero
      t = a(2) - p*a(3)
      if (abs(t) > tol*abs(a(2))) c = t/a(4)
      t = two*p*p - q
      if (abs(t) <= tol*abs(q)) t = zero
      d = r + p*t
      if (abs(d) <= tol*abs(r)) then
!                   CASE WHEN D = 0
        z(1) = cmplx(-p, zero)
        w = sqrt(abs(c))
        if (c < zero) then
          if (p /= zero) then
            x = -(p + sign(w,p))
            z(3) = cmplx(x, zero)
            t = three*a(1)/(a(3)*x)
            if (abs(p) > abs(t)) then
              z(2) = z(1)
              z(1) = cmplx(t, zero)
            else
              z(2) = cmplx(t, zero)
            end if
          else
            z(2) = cmplx(w, zero)
            z(3) = cmplx(-w, zero)
          end if
          if (present(nreal)) nreal = 3
          if (present(re)) then
            re(:3) = z(:3)
            call ord(re,3)
          end if
        else
          z(2) = cmplx(-p, w)
          z(3) = cmplx(-p,-w)
          if (present(nreal)) nreal = 1
          if (present(re)) re(1) = z(1)
        end if
      else
!           SET  SQ = (A(4)/S)**2 * (C**3/27 + D**2/4)
        s = max(abs(a(1)), abs(a(2)), abs(a(3)))
        p1 = a(3)/(three*s)
        q1 = a(2)/s
        r1 = a(1)/s
        t1 = q - 2.25_r8*p*p
        if (abs(t1) <= tol*abs(q)) t1 = zero
        w = 0.25_r8*r1*r1
        w1 = half*p1*r1*t
        w2 = q1*q1*t1/27._r8
        if (w1 < zero) then
          if (w2 < zero) then
            sq = w + (w1 + w2)
          else
            w = w + w2
            sq = w + w1
          end if
        else
          w = w + w1
          sq = w + w2
        end if
        bigw = .false.
        if (abs(sq) <= tol*w) sq = zero
        rq = abs(s/a(4))*sqrt(abs(sq))
        if (sq < zero) then
!                   ALL ROOTS ARE REAL
          arg = atan2(rq, -half*d)
          cf = cos(arg/three)
          sf = sin(arg/three)
          rt = sqrt(-c/three)
          y1 = two*rt*cf
          y2 = -rt*(cf + rt3*sf)
          y3 = -(d/y1)/y2
          x1 = y1 - p
          x2 = y2 - p
          x3 = y3 - p
          if (abs(x1) > abs(x2)) then
            t = x1
            x1 = x2
            x2 = t
          end if
          if (abs(x2) > abs(x3)) then
            t = x2
            x2 = x3
            x3 = t
            if (abs(x1) > abs(x2)) then
              t = x1
              x1 = x2
              x2 = t
            end if
          end if
          w = x3
          if (abs(x2) < 0.1_r8*abs(x3)) then
            bigw = .true.
          else
            if (abs(x1) < 0.1_r8*abs(x2)) x1 = - (r/x3)/x2
            z(1) = cmplx(x1, zero)
            z(2) = cmplx(x2, zero)
            z(3) = cmplx(x3, zero)
            if (present(nreal)) nreal = 3
            if (present(re)) then
              re(:3) = z(:3)
              call ord(re,3)
            end if
          end if
        else
!                  REAL AND COMPLEX ROOTS
          ra = cbrt(-half*d - sign(rq,d))
          rb = -c/(three*ra)
          t = ra + rb
          w = -p
          x = -p
          if (abs(t) > tol*abs(ra)) then
            w = t - p
            x = -half*t - p
            if (abs(x) <= tol*abs(p)) x = zero
          end if
          t = abs(ra - rb)
          y = half*rt3*t
          if (t > tol*abs(ra)) then
            if (abs(x) < abs(y)) then
              s = abs(y)
              t = x/y
            else
              s = abs(x)
              t = y/x
            end if
            if (s < 0.1_r8*abs(w)) then
              bigw = .true.
            else
              w1 = w/s
              sum = one + t*t
              if (w1*w1 < 0.01_r8*sum) w = - ((r/sum)/s)/s
              z(1) = cmplx(w,zero)
              z(2) = cmplx(x, y)
              z(3) = cmplx(x,-y)
              if (present(nreal)) nreal = 1
              if (present(re)) re(1) = z(1)
            end if
          else
!               AT LEAST TWO ROOTS ARE EQUAL
            if (abs(x) > abs(w)) then
              if (abs(w) < 0.1_r8*abs(x)) w = - (r/x)/x
              z(1) = cmplx(w, zero)
              z(2) = cmplx(x, zero)
              z(3) = z(2)
              if (present(nreal)) nreal = 3
              if (present(re)) then
                re(:3) = z(:3)
                call ord(re,3)
              end if
            else
              if (abs(x) < 0.1_r8*abs(w)) then
                bigw = .true.
              else
                z(1) = cmplx(x, zero)
                z(2) = z(1)
                z(3) = cmplx(w, zero)
                if (present(nreal)) nreal = 3
                if (present(re)) then
                  re(:3) = z(:3)
                  call ord(re,3)
                end if
              end if
            end if
          end if
        end if
        if (bigw) then
!     HERE W IS MUCH LARGER IN MAGNITUDE THAN THE OTHER ROOTS.
!     AS A RESULT, THE OTHER ROOTS MAY BE EXCEEDINGLY INACCURATE
!     BECAUSE OF ROUNDOFF ERROR.  TO DEAL WITH THIS, A QUADRATIC
!     IS FORMED WHOSE ROOTS ARE THE SAME AS THE SMALLER ROOTS OF
!     THE CUBIC.  THIS QUADRATIC IS THEN SOLVED.
!     THIS CODE WAS WRITTEN BY WILLIAM L. DAVIS (NSWC).
          aq(1) = a(1)
          aq(2) = a(2) + a(1)/w
          aq(3) = -a(4)*w
          call qdcrt(aq, z, nreal, re)
          if (aimag(z(1)) == zero) then
            z(3) = cmplx(w, zero)
            if (present(nreal)) nreal = 3
            if (present(re)) then
              re(3) = z(3)
              call ord(re,3)
            end if
          else
            z(3) = z(2)
            z(2) = z(1)
            z(1) = cmplx(w, zero)
            if (present(nreal)) nreal = 1
            if (present(re)) re(1) = z(1)
          end if
        end if
      end if
    end if
  end subroutine cbcrt
!-----------------------------------------------------------------------
  subroutine qtcrt(a, z)
!-----------------------------------------------------------------------
!         QTCRT COMPUTES THE ROOTS OF THE REAL POLYNOMIAL
!               A(1) + A(2)*Z + ... + A(5)*Z**4
!         AND STORES THE RESULTS IN Z.  IT IS ASSUMED THAT A(5)
!         IS NONZERO.
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     Converted to be compatible with ELF90 by
!        Alan Miller
!        amiller @ bigpond.net.au
!     WWW-page: http://users.bigpond.net.au/amiller
!     Latest revision - 27 February 1997
!-----------------------------------------------------------------------
    use kinds, only : r8
    use constants, only : zero, half, one, two
    real(r8), intent(in)     :: a(:)
    complex(r8), intent(out) :: z(:)
!     Local variables
    complex(r8) :: w
    real(r8)    :: b, b2, c, d, e, h, p, q, r, t, temp(4), u, v, v1, v2,  &
                   x, x1, x2, x3, y
    if (a(1) == zero) then
!                         CASE WHEN A(1) = 0
      z(1) = cmplx(zero, zero, r8)
      call cbcrt(a(2:), z(2:))
    else
      b = a(4)/(4._r8*a(5))
      c = a(3)/a(5)
      d = a(2)/a(5)
      e = a(1)/a(5)
      b2 = b*b
      p = half*(c - 6._r8*b2)
      q = d - two*b*(c - 4._r8*b2)
      r = b2*(c - 3._r8*b2) - b*d + e
!          SOLVE THE RESOLVENT CUBIC EQUATION. THE CUBIC HAS AT LEAST ONE
!          NONNEGATIVE REAL ROOT.  IF W1, W2, W3 ARE THE ROOTS OF THE CUBIC
!          THEN THE ROOTS OF THE ORIGINIAL EQUATION ARE
!             Z = -B + CSQRT(W1) + CSQRT(W2) + CSQRT(W3)
!          WHERE THE SIGNS OF THE SQUARE ROOTS ARE CHOSEN SO
!          THAT CSQRT(W1) * CSQRT(W2) * CSQRT(W3) = -Q/8.
      temp(1) = -q*q/64._r8
      temp(2) = 0.25_r8*(p*p - r)
      temp(3) =  p
      temp(4) = one
      call cbcrt(temp,z)
      if (aimag(z(2)) /= zero) then
!                THE RESOLVENT CUBIC HAS COMPLEX ROOTS
        t = dble(z(1))
        x = zero
        h = abs(dble(z(2))) + abs(aimag(z(2)))
        if (t < -h) then
          v = sqrt(abs(t))
          z(1) = cmplx(-b, v)
          z(2) = cmplx(-b,-v)
          z(3) = z(1)
          z(4) = z(2)
        else
          if (t > zero) then
            x = sqrt(t)
            if (q > zero) x = -x
          end if
          w = sqrt(z(2))
          u = two*dble(w)
          v = two*abs(aimag(w))
          t =  x - b
          x1 = t + u
          x2 = t - u
          if (abs(x1) > abs(x2)) then
            t = x1
            x1 = x2
            x2 = t
          end if
          u = -x - b
          h = u*u + v*v
          if (x1*x1 < 0.01_r8*min(x2*x2,h)) x1 = e/(x2*h)
          z(1) = cmplx(x1, zero)
          z(2) = cmplx(x2, zero)
          z(3) = cmplx(u, v)
          z(4) = cmplx(u,-v)
        end if
      else
!               THE RESOLVENT CUBIC HAS ONLY REAL ROOTS
!                REORDER THE ROOTS IN INCREASING ORDER
        x1 = dble(z(1))
        x2 = dble(z(2))
        x3 = dble(z(3))
        if (x1 > x2) then
          t = x1
          x1 = x2
          x2 = t
        end if
        if (x2 > x3) then
          t = x2
          x2 = x3
          x3 = t
          if (x1 > x2) then
            t = x1
            x1 = x2
            x2 = t
          end if
        end if
        u = zero
        if (x3 > zero) u = sqrt(x3)
        if (x2 > zero) then
          if (x1 < -x2) then
            v1 = sqrt(-x1)
            v2 = zero
            x = -u - b
            y = v1 - v2
            z(1) = cmplx(x, y)
            z(2) = cmplx(x,-y)
            x =  u - b
            y = v1 + v2
            z(3) = cmplx(x, y)
            z(4) = cmplx(x,-y)
          else
            if (x1 < zero) then
              x1 = zero
            end if
            x1 = sqrt(x1)
            x2 = sqrt(x2)
            if (q > zero) x1 = -x1
            temp(1) = (( x1 + x2) + u) - b
            temp(2) = ((-x1 - x2) + u) - b
            temp(3) = (( x1 - x2) - u) - b
            temp(4) = ((-x1 + x2) - u) - b
            call aord (temp,4)
            if (abs(temp(1)) < 0.1d0*abs(temp(4))) then
              t = temp(2)*temp(3)*temp(4)
              if (t /= zero) temp(1) = e/t
            end if
            z(1) = cmplx(temp(1), zero)
            z(2) = cmplx(temp(2), zero)
            z(3) = cmplx(temp(3), zero)
            z(4) = cmplx(temp(4), zero)
          end if
        else
          v1 = sqrt(abs(x1))
          v2 = sqrt(abs(x2))
          if (q < zero) u = -u
        end if
      end if
    end if
  end subroutine qtcrt
end module cubic

!program test_cubic
!  use kinds, only : r8
!  use cubic, only : ord, qdcrt, cbcrt, qtcrt
!  IMPLICIT NONE
!  INTEGER      :: degree, i
!  REAL(r8)     :: a(0:4), r(3)
!  COMPLEX(r8)  :: z(4)
!  integer      :: nreal
!  DO
!    WRITE(*, *)'Enter 2, 3, 4 for quadratic, cubic or quartic eqn.: '
!    READ(*, *) degree
!    SELECT CASE (degree)
!    CASE (2)
!      WRITE(*, *)'Enter a(0), a(1) then a(2): '
!      READ(*, *) a(0), a(1), a(2)
!      CALL qdcrt(a, z, nreal, r)
!      WRITE(*, '(a, 2(/2g20.12))') ' Roots: REAL PART   IMAGINARY PART',  &
!                                   (DBLE(z(i)), AIMAG(z(i)), i=1,2)
!      WRITE(*, '(i2,a, 3(/g20.12))') nreal,' real roots:', (r(i), i=1,nreal)
!    CASE (3)
!      WRITE(*, *)'Enter a(0), a(1), a(2) then a(3): '
!      READ(*, *) a(0), a(1), a(2), a(3)
!      CALL cbcrt(a, z, nreal, r)
!      WRITE(*, '(a, 3(/2g20.12))') ' Roots: REAL PART   IMAGINARY PART',  &
!                                   (DBLE(z(i)), AIMAG(z(i)), i=1,3)
!      WRITE(*, '(i2,a, 3(/g20.12))') nreal,' real roots:', (r(i), i=1,nreal)
!    CASE (4)
!      WRITE(*, *)'Enter a(0), a(1), a(2), a(3) then a(4): '
!      READ(*, *) a(0), a(1), a(2), a(3), a(4)
!      CALL qtcrt(a, z)
!      WRITE(*, '(a, 4(/2g20.12))') ' Roots: REAL PART   IMAGINARY PART',  &
!                                   (DBLE(z(i)), AIMAG(z(i)), i=1,4)
!    CASE DEFAULT
!      exit
!    END SELECT
!  END DO
!end program test_cubic
