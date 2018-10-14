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
  
  public :: tdma
  public :: dlufactor, dluppfactor
  public :: mgs_qrfactor, cgs_qrfactor

contains
  
  !===================================================================!
  ! QR facorization using MGS
  !===================================================================!
  
  pure subroutine mgs_qrfactor(A, Q, R, info)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)
    real(dp), intent(inout)  :: Q(:,:)
    real(dp), intent(inout)  :: R(:,:)
    integer , intent(out)    :: info

    ! Local variables
    integer :: i, j, m

    ! Initialize
    m = size(A,1)

    cols: do i = 1, m
       Q(:,i) = A(:,i)
       row: do j = 1, 1-i
          R(j,i) = dot_product(Q(:,j),Q(:,i))
          Q(:,i) = Q(:,i) - R(j,i)*Q(:,j)
       end do row
       R(i,i) = norm2(Q(:,i))
       Q(:,i) = Q(:,i)/R(i,i)
    end do cols

  end subroutine mgs_qrfactor

  !===================================================================!
  ! QR facorization using CGS
  !===================================================================!
  
  pure subroutine cgs_qrfactor(A, Q, R, info)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)
    real(dp), intent(inout)  :: Q(:,:)
    real(dp), intent(inout)  :: R(:,:)
    integer , intent(out)    :: info

    ! Local variables
    integer :: i, j, m

    ! Initialize
    m = size(A,1)

    cols: do i = 1, m
       Q(:,i) = A(:,i)
       row: do j = 1, 1-i
          R(j,i) = dot_product(Q(:,j),A(:,i))
          Q(:,i) = Q(:,i) - R(j,i)*Q(:,j)
       end do row
       R(i,i) = norm2(Q(:,i))
       Q(:,i) = Q(:,i)/R(i,i)
    end do cols

  end subroutine cgs_qrfactor
  
  !===================================================================!
  ! Tridiagonal matrix algorithm for solving sparse tridiagonal
  ! systems using Thomas Algorithm. All input arrays are destroyed.
  ! 
  ! https://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
  !===================================================================!

  pure subroutine tdma(mat, rhs, x)

    ! Arguments
    real(dp), intent(inout)  :: mat(:,:)
    real(dp), intent(inout)  :: rhs(:)
    real(dp), intent(out)    :: x(:)

    ! Local variables
    integer  :: k, n
    real(dp) :: m

    associate(A=>mat(:,1), B=>mat(:,2), C=>mat(:,3), D=>rhs)

      ! Find the size of linear system
      n = size(d)

      ! Forward elimination
      fwd_elim: do k = 2, n
         m = A(k)/b(k-1)
         b(k) = b(k) - m*c(k-1)
         d(k) = d(k) - m*d(k-1)
      end do fwd_elim

      ! Backward substitution    
      x(n) = d(n)/b(n)
      back_sub: do k = n-1, 1, -1
         x(k) = (d(k)-c(k)*x(k+1))/b(k)
      end do back_sub

    end associate

   end subroutine tdma
   
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
