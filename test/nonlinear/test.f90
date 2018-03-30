!=====================================================================!
! Module containing test problems for Newton's method
!=====================================================================!

module test_problems

  implicit none

contains
  
  ! f = 2x^2 -5
  pure function f1(x)    
    real(8), intent(in) :: x(:)
    real(8) :: f1(size(x))
    f1 = 2.0d0*x*x - 5.0d0    
  end function f1

  ! dfdx = 4x
  pure function dfdx1(x)
    real(8), intent(in) :: x(:)
    real(8) :: dfdx1(size(x),size(x))
    dfdx1(1,1) = 4.0d0*x(1)
  end function dfdx1

  ! f = sin(x) + x
  pure function f2(x)    
    real(8), intent(in) :: x(:)
    real(8) :: f2(size(x))
    f2 = sin(x) + x
  end function f2

  ! dfdx = cos(x) + 1
  pure function dfdx2(x)
    real(8), intent(in) :: x(:)
    real(8) :: dfdx2(size(x),size(x))
    dfdx2(1,1) = cos(x(1)) + 1.0d0
  end function dfdx2

  ! f = cos(x)
  pure function f3(x)    
    real(8), intent(in) ::  x(:)
    real(8) :: f3(size(x))
    f3 = cos(x)
  end function f3

  ! dfdx = sin(x)
  pure function dfdx3(x)
    real(8), intent(in) :: x(:)
    real(8) :: dfdx3(size(x),size(x))
    dfdx3(1,1) = sin(x(1))
  end function dfdx3
  
  function res(x)
    real(8), intent(in) :: x(:)
    real(8) :: res(size(x))
  end function res
  
  function jac(x)
    real(8), intent(in) :: x(:)
    real(8) :: jac(size(x), size(x))
  end function jac

end module test_problems

!=====================================================================!
! Main program to solve the test problems using different methods
!=====================================================================!

program test_nonlinear

  use test_problems, only : f1, f2, f3, dfdx1, dfdx2, dfdx3
  ! use nonlinear, only : newton

  implicit none

  test_newton: block

    ! Test newton solver

  end block test_newton

end program test_nonlinear
