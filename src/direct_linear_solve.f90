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

  public :: dlufactor, dluppfactor

contains
  
  !===================================================================!
  ! Plain vanilla LU factorization algorithm without pivoting. How to
  ! store L and U into the original matrix (like LAPACK)?
  !===================================================================!
  
  pure subroutine dlufactor(A, L, U, info)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)
    real(dp), intent(inout)  :: L(:,:)
    real(dp), intent(inout)  :: U(:,:)
    integer , intent(out)    :: info

    ! Local variables
    integer :: k, j, m

    ! Initialize
    m = size(A,1)
    U = A 
    do concurrent(k=1:m)
       L(k,k) = 1.0d0
    end do

    ! Algorithm
    do k = 1, m-1
       do j = k + 1, m
          L(j,k) = U(j,k)/U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
       end do
    end do

  end subroutine dlufactor

  !===================================================================!
  ! Plain vanilla LU factorization algorithm with partial pivoting
  !===================================================================!
  
  subroutine dluppfactor(A, L, U, P, info)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)
    real(dp), intent(inout)  :: L(:,:)
    real(dp), intent(inout)  :: U(:,:)
    real(dp), intent(inout)  :: P(:,:)
    integer , intent(out)    :: info

    ! Local variables
    integer :: k, j, m, pivot, i, jj

    ! Initialize
    m = size(A,1)
    U = A 
    do concurrent(k=1:m)
       L(k,k) = 1.0d0
       P(k,k) = 1.0d0
    end do
    
    ! Algorithm
    do k = 1, m - 1      

       ! Pivoting logic
       pivot = maxloc(abs(U(:,k)), dim = 1)

       ! Exchange rows
       U(k,k:m)   = U(pivot,k:m)
       L(k,1:k-1) = L(pivot,1:k-1)
       P(k,:)     = P(pivot,:)

       ! LU logic
       do j = k + 1, m
          L(j,k) = U(j,k)/U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
       end do

    end do

  end subroutine dluppfactor

!!$  !===================================================================!
!!$  ! Plain vanilla LU factorization algorithm without pivoting
!!$  !===================================================================!
!!$
!!$  pure subroutine dgetlu(A, info)
!!$
!!$    real(dp), intent(inout)  :: A(:,:)
!!$    integer , intent(out)    :: info
!!$
!!$  end subroutine dgetlu

end module direct_linear_solve
