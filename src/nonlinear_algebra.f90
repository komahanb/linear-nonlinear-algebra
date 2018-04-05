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

  public :: newton

contains
  
  subroutine newton(F, FPRIME, tau_r, tau_a, x)

    real(8), intent(in) :: tau_r, tau_a
    real(8), intent(inout) :: x(:)

    ! residual and jacobian
    procedure :: F
    procedure :: FRIME

    ! local variables
    real(8) :: r0 
    real(8) :: jac(size(x),size(x))
    real(8) :: s(size(x))
    
    ! Initial residual
    r0 = norm2(F(x))

    do while (norm2(F(x)) > tau_r*r0 + tau_a)

       ! Compute Jacobian
       jac = FPRIME(x)

       ! Solve the linear system
       s =  solve(jac, -F(x))

       ! Apply the update
       x = x + s

    end do

  end subroutine newton

end module nonlinear_algebra
