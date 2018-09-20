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

  public :: dlufactor, dluppfactor, thomas

contains
 
  !============================================================
  ! Thomas algorithm for efficent solution of banded systems
  !===========================================================
  
   pure subroutine thomas(bandwidth, A, b)

    implicit none 

    ! Arguments
    integer , intent(in)    :: bandwidth
    real(dp), intent(inout) :: A(:,:), b(:)

    ! Local variables
    real(dp) :: coeff
    integer  :: i
    integer  :: n 

    ! Find the number of unknowns
    n = size(b)

    ! TODO generalize this for pentadiagonal matrix too

    ! Forward elimination
    do i = 2, n
       coeff  = A(i,1)/A(i-1,2)
       A(i,2) = A(i,2)-coeff*A(i-1,3)
       b(i)   = b(i)-A(i,1)*b(i-1)
    end do

    ! Back substitution
    b(n) = b(n)/A(n,2)
    do i = n-1, 1, -1
       b(i) = (b(i)- A(i,3)*b(i+1))/A(i,2)
    end do
    
  end subroutine thomas

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
    do k = 1, m - 1
       do j = k + 1, m
          L(j,k) = U(j,k)/U(k,k)
          U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m)
       end do
    end do

  end subroutine dlufactor

  !===================================================================!
  ! Swap two arrays
  !===================================================================!
  pure subroutine swap(a, b)
    real(dp), intent(inout), dimension(:) :: a, b
    real(dp), dimension(size(a)) :: work
    work = a
    a = b
    b = work
  end subroutine swap

  !===================================================================!
  ! LU factorization algorithm with partial pivoting
  !===================================================================!
  
  pure subroutine dluppfactor(A, L, U, P, info)

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
    do k = 1, m-1

       ! Find the entry in the left column with the largest absolute
       ! value. This entry is called the pivot.
       pivot = maxloc(abs(U(k:m,k)), dim = 1)
       i = pivot + k - 1
       
       ! Exchange rows
       call swap(U(k,k:m)   , U(i,k:m))
       call swap(L(k,1:k-1) , L(i,1:k-1))
       call swap(P(k,:)     , P(i,:))

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
