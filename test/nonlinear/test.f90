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
    dfdx3(1,1) = -sin(x(1))
  end function dfdx3

  ! Chandrasekhar H equations
  function chandra_res(x)

    real(8), intent(in) :: x(:)
    real(8)  :: chandra_res(size(x))    
    integer, parameter  :: N = 200
    real(8), parameter  :: c = 0.9d0

    if ( n .ne. 200) stop "variable-mismatch"
    
    assemble: block
      
      real(8) :: mui, muj, alpha, beta
      integer :: i, j

      do i = 1, N
         
         mui = (dble(i) - 0.5d0)/dble(N)

         ! Find the term within integral
         alpha = 0.0d0
         do j = 1, N
            muj = (dble(j) - 0.5d0)/dble(N)
            alpha = mui*x(j)/(mui+muj)
         end do

         ! Evaluate the inverse term
         beta = 1.0d0 - c*alpha/(2.0d0*dble(N))

         ! Assmble residual
         chandra_res(i) = x(i) - 1.0d0/beta
         
      end do

    end block assemble

  end function chandra_res

  !===================================================================!
  ! Finite difference approximation for jacobian
  !===================================================================!
  
  function chandra_jac(x)

    use nonlinear_algebra, only : diffjac
    
    real(8), intent(in) :: x(:)
    real(8) :: chandra_jac(size(x),size(x))
    integer, parameter :: N = 200
    
    call diffjac(x, chandra_res, chandra_jac)
    
  end function chandra_jac

end module test_problems

!=====================================================================!
! Main program to solve the test problems using different methods
!=====================================================================!

program test_nonlinear
  
  use test_problems, only : f1, f2, f3, dfdx1, dfdx2, dfdx3, &
       & chandra_res, chandra_jac
  
  use nonlinear_algebra, only : newton, secant, chord, fixed_point, &
       & shamanskii

  implicit none
  
  real(8), parameter :: tau_r = 1.0d-6
  real(8), parameter :: tau_a = 1.0d-6
  integer, parameter :: maxit = 100
  
  integer :: iter, flag
  real(8) :: tol

  ! Test nonlinear solvers on chandrasekhar method
  test_chandra: block

    integer, parameter :: npts = 200
    real(8) :: x(npts,4)

    print *, "Chandrasekhar Equation using Newton method"
  
    x(:,1) = 2.0d0
    call newton(chandra_res, chandra_jac, tau_r, tau_a, maxit, x(:,1))
    !print *, x(:,1)

    print *, "Chandrasekhar Equation using chord method"
    x(:,2) = 1.0d0
    call chord(chandra_res, chandra_jac, tau_r, tau_a, maxit, x(:,2))
    !print *, x(:,2)

    print *, "Chandrasekhar Equation using fixed point method"
    x(:,3) = 1.0d0
    call fixed_point(chandra_res, tau_r, tau_a, maxit, x(:,3))
    !print *, x(:,3)

    print *, "Chandrasekhar Equation using shamanskii method"
    x(:,4) = 1.0d0
    call shamanskii(chandra_res, chandra_jac, 2, tau_r, tau_a, maxit, x(:,4))
    !print *, x(:,4)

  end block test_chandra

  test_f1 : block

    integer, parameter :: npts = 1
    real(8) :: x(npts), x0(npts)

    print *, "function 1 using newton method"
    x = 10.0d0
    call newton(f1, dfdx1, tau_r, tau_a, maxit, x)
    print *, x

    print *, "function 1 using secant method"
    x = 10.0d0
    x0 = 0.99d0*x
    call secant(f1, tau_r, tau_a, maxit, x0, x)
    print *, x

    print *, "function 1 using chord method"
    x = 10.0d0
    call chord(f1, dfdx1, tau_r, tau_a, maxit, x)
    print *, x

  end block test_f1

  test_f2: block

    integer, parameter :: npts = 1
    real(8) :: x(npts), x0(npts)

    print *, "function 2 using newton method"
    x = 0.5d0
    call newton(f2, dfdx2, tau_r, tau_a, maxit, x)
    print *, x

    print *, "function 2 using secant method"
    x = 0.5d0
    x0 = 0.99d0*x
    call secant(f2, tau_r, tau_a, maxit, x0, x)
    print *, x

    print *, "function 2 using chord method"
    x = 0.5d0
    call chord(f2, dfdx2, tau_r, tau_a, maxit, x)
    print *, x

  end block test_f2
  
  test_f3: block

    integer, parameter :: npts = 1
    real(8) :: x(npts), x0(npts)

    print *, "function 3 using newton method"
    x = 3.0d0
    call newton(f3, dfdx3, tau_r, tau_a, maxit, x)
    print *, x

    print *, "function 3 using secant method"
    x = 3.0d0
    x0 = 0.99d0*x
    call secant(f3, tau_r, tau_a, maxit, x0, x)
    print *, x

    print *, "function 3 using chord method"
    x = 3.0d0
    call chord(f3, dfdx3, tau_r, tau_a, maxit, x)
    print *, x

  end block test_f3

end program test_nonlinear
