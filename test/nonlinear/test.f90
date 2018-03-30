!=====================================================================!
! Module containing test problems for Newton's method
!=====================================================================!

module test_problems

  implicit none

contains
  
  ! f = 2x^2 -5
  real(8) function f1(x)    
    real(8), intent(in) :: x
    f1 = 2.0d0*x*x - 5.0d0    
  end function f1

  ! dfdx = 4x
  real(8) function dfdx1(x)    
    real(8), intent(in) :: x
    dfdx1 = 4.0d0*x
  end function dfdx1

  ! f = sin(x) + x
  real(8) function f2(x)    
    real(8), intent(in) ::  x   
    f2 = sin(x) + x
  end function f2

  ! dfdx = cos(x) + 1
  real(8) function dfdx2(x)    
    real(8), intent(in) :: x
    dfdx2 = cos(x) + 1.0d0
  end function dfdx2

  ! f = cos(x)
  real(8) function f3(x)    
    real(8), intent(in) ::  x  
    f3 = cos(x)
  end function f3

  ! dfdx = sin(x)
  real(8) function dfdx3(x)
    real(8), intent(in) :: x
    dfdx3 = sin(x)
  end function dfdx3

end module test_problems

!=====================================================================!
! Main program to solve the test problems using different methods
!=====================================================================!

program test_nonlinear

  use test_problems, only : f1, f2, f3
  ! use nonlinear, only : newton

  implicit none

  test_newton: block

    ! Test newton solver

  end block test_newton

end program test_nonlinear
