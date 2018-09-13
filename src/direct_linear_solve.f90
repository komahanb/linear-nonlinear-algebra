!=====================================================================!
! Module that solves linear systems using direct methods
!=====================================================================!

module direct_linear_solve

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64

  ! disable implicit datatypes
  implicit none

  ! restrict access to all functions
  private

  public :: dgetlu

contains

  !===================================================================!
  ! Plain vanilla LU factorization algorithm without pivoting
  !===================================================================!

  pure subroutine dgetlu(A, info)

    real(dp), intent(inout)  :: A(:,:)
    integer , intent(out)    :: info

  end subroutine dgetlu

end module direct_linear_solve
