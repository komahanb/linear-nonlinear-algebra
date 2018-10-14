!=====================================================================!
! Run gradient-based descent algorithms to solve linear algebra
!=====================================================================!

program test

  use direct_linear_solve, only : mgs_qrfactor, cgs_qrfactor  

  implicit none
  
  real(8) :: A(3,3), B(3,3)
  real(8) :: Q(3,3), R(3,3)
  integer :: i, info
  
  call random_number(A)
  
  do i = 1, 3
     print *, A(i,:)
  end do
  
  call mgs_qrfactor(A, Q, R, info)

  B = matmul(Q,R)

  print *, ''
  do i = 1, 3
     print *, B(i,:)
  end do
  
  
end program test
