module intervals
  use kinds, only : r8
  use constants, only : huge_r8
  implicit none
  private :: r8, huge_r8
  type interval
    real(kind=r8) :: lo, hi
  end type interval
  private :: interval2real, real2interval
  interface assignment(=)
    module procedure interval2real, real2interval
  end interface
  private :: intersect_r8
  interface intersect
    module procedure intersect_r8
  end interface
contains
  function intersect_r8(a,b)
    implicit none
    type(interval), intent(in) :: a, b
    type(interval) :: at, bt, intersect_r8
    if (a%lo > a%hi) then
      at%lo = a%hi
      at%hi = a%lo
    else
      at = a
    end if
    if (b%lo > b%hi) then
      bt%lo = b%hi
      bt%hi = b%lo
    else
      bt = b
    end if
    if (at%lo > bt%lo) then
      intersect_r8%lo = at%lo
    else
      intersect_r8%lo = bt%lo
    end if
    if (at%hi < bt%hi) then
      intersect_r8%hi = at%hi
    else
      intersect_r8%hi = bt%hi
    end if
    if (intersect_r8%hi < intersect_r8%lo) then
      intersect_r8%lo = -huge_r8
      intersect_r8%hi = -huge_r8
    end if
  end function intersect_r8
  subroutine interval2real(a,b)
    implicit none
    real(kind=r8), intent(out) :: a
    type(interval), intent(in) :: b
    a = (b%lo + b%hi)/2._r8
  end subroutine interval2real
  subroutine real2interval(a,b)
    implicit none
    type(interval), intent(out) :: a
    real(kind=r8), intent(in) :: b
    a%lo = b
    a%hi = b
  end subroutine real2interval
end module intervals
