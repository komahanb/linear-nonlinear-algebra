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

  public :: newton, chord, secant

contains
  
  subroutine newton(F, FPRIME, tau_r, tau_a, max_it, x)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x(:)
   
    interface
       pure function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       pure function FPRIME(x)
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

  subroutine secant(F, tau_r, tau_a,  max_it, x0, x1)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x0(:), x1(:)
   
    interface
       pure function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       pure function FPRIME(x)
         real(8), intent(in) :: x(:)
         real(8) :: FPRIME(size(x),size(x))
       end function FPRIME
    end interface

    ! local variables
    real(8) :: r0 
    real(8) :: jac(size(x0),size(x0))
    real(8) :: s(size(x0))
    integer :: iter
    integer :: nvars, i
    real(8), parameter :: h = 1.0e-8
    real(8) :: xtmp(size(x0))
    
    nvars = size(x0)
    
    ! Initial residual
    r0 = norm2(F(x1))

    iter = 0
    
    do while (norm2(F(x1)) > tau_r*r0 + tau_a .and. iter .le. max_it)

       ! Increment the iteration count
       iter = iter + 1
!!$
!!$       ! Compute Jacobian with first order approx
!!$       do i = 1, nvars
!!$          
!!$          ! Increment the iteration count
!!$          iter = iter + 1
!!$          
!!$          ! Perturb x
!!$          xhat(i) = x(i) + h
!!$
!!$          ! Evaluate column
!!$          jac(:,i) = (F(xhat) - F(x))/h
!!$
!!$          ! Restore x
!!$          xhat(i) = x(i)
!!$          
!!$       end do

       ! Solve the linear system
       ! s = solve(jac, F(x))

       ! Apply the update
       xtmp = x1 
       x1 = x1 - F(x1)*(x1-x0)/(F(x1)-F(x0))
       x0 = xtmp
       
       ! print details
       print *, iter, norm2(F(x1))

    end do

  end subroutine secant
  
  subroutine chord(F, FPRIME, tau_r, tau_a, max_it, x)

    real(8), intent(in) :: tau_r, tau_a
    integer, intent(in) :: max_it
    real(8), intent(inout) :: x(:)
    
    interface
       pure function F(x)
         real(8), intent(in) ::  x(:)
         real(8) :: F(size(x))
       end function F
       pure function FPRIME(x)
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

end module nonlinear_algebra
