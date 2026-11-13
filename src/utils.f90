module utils
  use kinds, only : wp_=>r8
  implicit none

contains

  function locatef(a,n,x) result(j)
! Given an array a(n), and a value x, with a(n) monotonic, either
! increasing or decreasing, returns a value j such that
! a(j) < x <= a(j+1) for a increasing, and such that
! a(j+1) < x <= a(j) for a decreasing.
! j=0 or j=n  indicate that x is out of range  (Numerical Recipes)
    implicit none
    integer, intent(in) :: n
    real(wp_), dimension(n), intent(in) :: a
    real(wp_), intent(in) :: x
    integer :: j
    integer :: jl,ju,jm
    logical :: incr
    jl=0
    ju=n+1
    incr=a(n)>a(1)
    do while ((ju-jl)>1)
      jm=(ju+jl)/2
      if(incr.eqv.(x>a(jm))) then
        jl=jm
      else
        ju=jm
      endif
    end do
    j=jl
  end function locatef

  subroutine locate_unord(a,n,x,j,m,nj)
    implicit none
    integer, intent(in) :: n,m
    integer, intent(out) :: nj
    real(wp_), dimension(n), intent(in) :: a
    real(wp_), intent(in) :: x
    integer, dimension(m), intent(inout) :: j
    integer :: i
    nj=0
    do i=1,n-1
      if (x>a(i).neqv.x>a(i+1)) then
        nj=nj+1
        if (nj<=m) j(nj)=i
      end if
    end do
  end subroutine locate_unord

  function intlinf(x1,y1,x2,y2,x) result(y)
    !linear interpolation
    !must be x1 != x2
    implicit none
    real(wp_),intent(in) :: x1,y1,x2,y2,x
    real(wp_) :: y
    real(wp_) :: a
    a=(x2-x)/(x2-x1)
    y=a*y1+(1.0_wp_-a)*y2
  end function intlinf

  subroutine order(p,q)
! returns p,q in ascending order
    implicit none
    real(wp_), intent(inout) :: p,q
    real(wp_) :: temp
    if (p>q) then
      temp=p
      p=q
      q=temp
    end if
  end subroutine order

  subroutine bubble(a,n)
! bubble sorting of array a
    implicit none
    integer, intent(in) :: n
    real(wp_), dimension(n), intent(inout) :: a
    integer :: i, j
    do i=1,n
      do j=n,i+1,-1
        call order(a(j-1), a(j))
      end do
    end do
  end subroutine bubble

function inside(xc,yc,n,x,y)
  implicit none
  integer, intent(in) :: n
  real(wp_), dimension(n), intent(in) :: xc,yc
  real(wp_), intent(in) :: x,y
  logical :: inside
  integer, dimension(n) :: jint
  real(wp_), dimension(n) :: xint
  real(wp_), dimension(n+1) :: xclosed,yclosed
  integer :: i,nj
  xclosed(1:n)=xc(1:n)
  yclosed(1:n)=yc(1:n)
  xclosed(n+1)=xc(1)
  yclosed(n+1)=yc(1)
  call locate_unord(yclosed,n+1,y,jint,n,nj)
  inside=.false.
  if (nj==0) return
  do i=1,nj
    xint(i)=intlinf(yclosed(jint(i)),xclosed(jint(i)), &
                   yclosed(jint(i)+1),xclosed(jint(i)+1),y)
  end do
  call bubble(xint,nj)
  inside=(mod(locatef(xint,nj,x),2)==1)
end function inside

end module utils

