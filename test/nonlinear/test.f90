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

    real(8), intent(in) :: x(200)
    real(8) :: chandra_res(size(x))

    assemble: block
      
      real(8), parameter :: c = 0.9d0
      integer, parameter :: N = 200
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

  function chandra_jac(x)
    real(8), intent(in) :: x(:)
    real(8) :: chandra_jac(size(x), size(x))
  end function chandra_jac

end module test_problems

!=====================================================================!
! Main program to solve the test problems using different methods
!=====================================================================!

program test_nonlinear
  
  use test_problems, only : f1, f2, f3, dfdx1, dfdx2, dfdx3, &
       & chandra_res, chandra_jac
  
  use nonlinear_algebra, only : newton, secant, chord

  implicit none
  
  integer, parameter :: npts = 1
  real(8), parameter :: tau_r = 1.0d-6
  real(8), parameter :: tau_a = 1.0d-6
  integer, parameter :: maxit = 100
  
  integer :: iter, flag
  real(8) :: x(npts)
  real(8) :: tol

  ! Uses proper jacobian (analytic)
  test_newton: block

    x = 10.0d0
    call newton(f1, dfdx1, tau_r, tau_a, maxit, x)
    print *, x

    x = 0.5d0
    call newton(f2, dfdx2, tau_r, tau_a, maxit, x)
    print *, x

    x = 3.0d0
    call newton(f3, dfdx3, tau_r, tau_a, maxit, x)
    print *, x

  end block test_newton

  ! approximates jacobian with finite difference
  test_secant: block

    x = 10.0d0
    call secant(f1, tau_r, tau_a, maxit, x)
    print *, x

    x = 0.5d0
    call secant(f2, tau_r, tau_a, maxit, x)
    print *, x

    x = 3.0d0
    call secant(f3, tau_r, tau_a, maxit, x)
    print *, x

  end block test_secant

  ! Uses jacobian at first iteration
  test_chord: block

    x = 10.0d0
    call chord(f1, dfdx1, tau_r, tau_a, maxit, x)
    print *, x

    x = 0.5d0
    call chord(f2, dfdx2, tau_r, tau_a, maxit, x)
    print *, x

    x = 3.0d0
    call chord(f3, dfdx3, tau_r, tau_a, maxit, x)
    print *, x

  end block test_chord

end program test_nonlinear
