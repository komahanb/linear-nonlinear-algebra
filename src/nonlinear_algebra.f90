!=====================================================================!
! Library of nonlinear solvers
!=====================================================================!

module nonlinear_algebra

  ! Import dependencies
  use iso_fortran_env, only : dp => REAL64

  ! Import linear solver
  use linear_algebra

  ! disable implicit datatypes
  implicit none

  ! all members are private by default
  private

  public :: newton, chord, secant, fixed_point, shamanskii
  public :: diffjac

contains

  !===================================================================!
  ! Form the jacobian using finite differences
  !===================================================================!
  
  subroutine diffjac(x, F, jac)
    
    real(8), intent(in) :: x(:)
    
    interface
        function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
    end interface

    real(8), intent(inout) :: jac(size(x),size(x))

    integer             :: nvars, i
    real(8), parameter  :: h = 1.0e-8
    real(8)             :: xhat(size(x))

    nvars = size(x)

    xhat = x
    
    ! Compute Jacobian with first order approx
    do i = 1, nvars

       ! Perturb x
       xhat(i) = x(i) + h

       ! Evaluate column
       jac(:,i) = (F(xhat) - F(x))/h

       ! Restore x
       xhat(i) = x(i)

    end do

  end subroutine diffjac

  !===================================================================!
  ! A method that features quadratic convergence
  !===================================================================!
  
  subroutine newton(F, FPRIME, tau_r, tau_a, max_it, x)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x(:)
   
    interface
       function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       function FPRIME(x)
         real(8), intent(in) :: x(:)
         real(8) :: FPRIME(size(x),size(x))
       end function FPRIME
    end interface

    ! local variables
    real(8) :: r0 
    real(8) :: jac(size(x),size(x))
    real(8) :: s(size(x))
    integer :: iter

    ! Initial residual
    r0 = norm2(F(x))

    iter = 0
    do while (norm2(F(x)) > tau_r*r0 + tau_a .and. iter .le. max_it)

       ! Increment the iteration count
       iter = iter + 1

       ! Compute Jacobian
       jac = FPRIME(x)

       ! Solve the linear system
       s = solve(jac, F(x))

       ! Apply the update
       x = x - s

       ! print details
       print *, iter, norm2(F(x))

    end do

  end subroutine newton

  !===================================================================!
  ! Use use supplied two points to form Jacobian (super-linear
  ! convergence)
  !===================================================================!
  
  subroutine secant(F, tau_r, tau_a,  max_it, x0, x1)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x0(:), x1(:)
   
    interface
       function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       function FPRIME(x)
         real(8), intent(in) :: x(:)
         real(8) :: FPRIME(size(x),size(x))
       end function FPRIME
    end interface

    ! local variables
    real(8) :: r0 
    real(8) :: jac(size(x0),size(x0))
    real(8) :: s(size(x0))
    integer :: iter
    real(8) :: xtmp(size(x0))
    
    ! Initial residual
    r0 = norm2(F(x1))

    iter = 0
    
    do while (norm2(F(x1)) > tau_r*r0 + tau_a .and. iter .le. max_it)

       ! Increment the iteration count
       iter = iter + 1

       ! Apply the update
       xtmp = x1 
       x1 = x1 - F(x1)*(x1-x0)/(F(x1)-F(x0))
       x0 = xtmp
       
       ! print details
       print *, iter, norm2(F(x1))

    end do

  end subroutine secant

  !===================================================================!
  ! Use the same jacobian each iteration (linear convergence)
  !===================================================================!
  
  subroutine chord(F, FPRIME, tau_r, tau_a, max_it, x)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x(:)
    
    interface
       function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       function FPRIME(x)
         real(8), intent(in) :: x(:)
         real(8) :: FPRIME(size(x),size(x))
       end function FPRIME
    end interface

    ! local variables
    real(8) :: r0 
    real(8) :: jac(size(x),size(x))
    real(8) :: s(size(x))
    integer :: iter

    ! Initial residual
    r0 = norm2(F(x))

    ! Compute the jacobian once with the initial iterate
    jac = FPRIME(x)

    iter = 0
    do while (norm2(F(x)) > tau_r*r0 + tau_a .and. iter .le. max_it )

       ! Increment the iteration count
       iter = iter + 1
       
       ! Solve the linear system
       s = solve(jac, F(x))

       ! Apply the update
       x = x - s

       ! print details
       print *, iter, norm2(F(x))

    end do

  end subroutine chord

  !===================================================================!
  ! Eliminate linear solve by using identity jacobian
  !===================================================================!
  
  subroutine fixed_point(F, tau_r, tau_a, max_it, x)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x(:)

    interface
       function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
    end interface

    ! local variables
    real(8) :: r0 
    integer :: iter

    ! Initial residual
    r0 = norm2(F(x))

    iter = 0
    do while (norm2(F(x)) > tau_r*r0 + tau_a .and. iter .le. max_it )

       ! Increment the iteration count
       iter = iter + 1

       ! Apply the update
       x = x - F(x)

       ! print details
       print *, iter, norm2(F(x))

    end do

  end subroutine fixed_point

  !===================================================================!
  ! Combines Newton and chord methods. This is efficient with LU
  ! factorization than the way its implemented here.
  !===================================================================!
  
  subroutine shamanskii(F, FPRIME, chord_iterations, tau_r, tau_a, max_it, x)

    real(8), intent(in)    :: tau_r, tau_a
    integer, intent(in)    :: max_it
    real(8), intent(inout) :: x(:)
    integer, intent(in)    :: chord_iterations

    interface
       function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       function FPRIME(x)
         real(8), intent(in) :: x(:)
         real(8) :: FPRIME(size(x),size(x))
       end function FPRIME
    end interface

    ! Local variables
    real(8) :: r0 
    real(8) :: jac(size(x),size(x))
    real(8) :: s(size(x))
    integer :: newton_iter, chord_iter

    ! Initial residual
    r0 = norm2(F(x))
    
    newton_iter = 0
    quad_newton: do while (norm2(F(x)) > tau_r*r0 + tau_a .and. newton_iter .le. max_it)
       
       ! Increment the newton_iteration count
       newton_iter = newton_iter + 1

       ! Compute Jacobian
       jac = FPRIME(x)

       ! superlinear part 
       chord_iter = 0
       suplin_chord: do while ( chord_iter .le. chord_iterations )

          ! Increment the iteration count
          chord_iter = chord_iter + 1

          ! Solve the linear system
          s = solve(jac, F(x))

          ! Apply the update
          x = x - s

          ! print details
          print *, "sub-chord", chord_iter, norm2(F(x))

       end do suplin_chord

       ! Solve the linear system
       s = solve(jac, F(x))

       ! Apply the update
       x = x - s

       ! print details
       print *, "newton-", newton_iter, norm2(F(x))

    end do quad_newton
    
  end subroutine shamanskii

end module nonlinear_algebra
