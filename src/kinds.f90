module kinds
  implicit none
  integer, parameter :: r4=selected_real_kind(6,37), r8=selected_real_kind(15,307)
  integer, parameter :: i1=selected_int_kind(2)    , i2=selected_int_kind(4)
  integer, parameter :: i4=selected_int_kind(9)    , i8=selected_int_kind(18)
end module kinds
