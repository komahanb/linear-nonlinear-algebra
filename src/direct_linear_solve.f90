!=====================================================================!
! Module that solves linear systems using direct methods
!=====================================================================!

module direct_linear_solve

  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  use matrix_utils, only : print 

  ! disable implicit datatypes
  implicit none

  ! restrict access to all functions
  private
  
  public :: tdma
  public :: dlufactor, dluppfactor
  public :: mgs_qrfactor, cgs_qrfactor, householder, banded_householder
  public :: householder_factorization, givens_factorization
  public :: givens, qriteration

contains
  
  !===================================================================!
  ! QR iteration to solve EVP
  !===================================================================!
  
  subroutine QRiteration(A, T, U, maxiters, tol)

    ! Arguments
    real(dp), intent(inout) :: A(:,:), U(:,:), T(:,:)
    integer , intent(in)    :: maxiters
    real(dp), intent(in)    :: tol

    ! Locals
    real(dp), allocatable   :: Q(:,:), R(:,:), D(:,:)
    real(dp) :: norm
    integer  :: k, i, n, m, iter
    real(dp) :: G(2,2)

    ! Initialize
    m = size(A,1)
    n = size(A,2)

    ! Copy A into T
    T = A

    ! Allocate space for QR factorization
    allocate(Q, R, D, mold=A)
    Q = 0
    R = 0
    D = 0
    
    ! Initialize Q with Identity
    U = 0
    do concurrent(i=1:n)
       U(i,i) = 1.0d0
    end do

    norm = huge(1.0d0)

    iter = 1
    do while (&
         & norm .gt. tol .or. &
         & iter .gt. maxiters)

       ! QR
       call givens_factorization(T, Q, R)
       call print(T)
       call print (Q)
       call print(R)


       stop
       ! A = RQ
       T = matmul(R,Q)

       ! Update transformation matrix
       U = matmul(U,Q)
       
       ! Elements below diagonal should have converged to zero to
       ! specified tolerance
       D = T
       do i = 1, n       
          D(i,i) = 0.0d0
       end do
       norm = norm2(D)

       print *, n, iter , norm

       iter = iter + 1
       
    end do

    deallocate(R,Q,D)

  end subroutine QRiteration

  !===================================================================!
  ! QR facorization using Householder algorithm for banded matrix. The
  ! upper triangular part of a contains R.
  !===================================================================!
  
   subroutine banded_householder(A, p, q)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)
    integer, intent(in) :: p, q
    
    ! Local variables
    real(dp), allocatable :: x(:), v(:)
    integer  :: k, j, n, m, nrows, ncols
    real(dp) :: scalar

    allocate(x, source = A(:,1))
    allocate(v, source = A(:,1))

    ! Initialize
    nrows = size(A,1)
    ncols = size(A,2)

    cols: do k = 1, ncols
       
       m = min(nrows, k + p)

       !print *, m, nrows, k+p, k, ncols
       ! Extract the column below the diagonal
       x(k:m) = A(k:m, k)

       ! print *, x(k:m)

       ! Fimd the reflection vector 
       v(k:m) = x(k:m)
       v(k) = v(k) + sign(norm2(x(k:m)),x(k))

       ! Normalize the reflection vector
       v(k:m) = v(k:m)/norm2(v(k:m))

       ! Perform householder transformation to zero out the lower
       ! entries except one for each column of the matrix
       do j = k , min(ncols,k+p+q) !limit to upper bandwidth
          !print *, 'tranforming col', j, 'rows ', k, m
          scalar = dot_product(v(k:m), A(k:m,j))
          A(k:m,j) = A(k:m,j) - 2.0_dp*scalar*v(k:m)
       end do

    end do cols

    deallocate(x,v)

  end subroutine banded_householder

  !===================================================================!
  ! QR facorization using Householder. The upper triangular part of a
  ! contains R.
  !===================================================================!
  
  pure subroutine householder(A)

    ! Arguments
    real(dp), intent(inout)  :: A(:,:)

    ! Local variables
    real(dp), allocatable :: x(:), v(:)
    integer  :: k, j, n, m
    real(dp) :: scalar

    allocate(x, source = A(:,1))
    allocate(v, source = A(:,1))

    ! Initialize
    m = size(A,1)
    n = size(A,2)

    cols: do k = 1, n

       x = 0
       v = 0
       
       ! Extract the column below the diagonal
       x(k:m) = A(k:m, k)

       ! Fimd the reflection vector 
       v(k:m) = x(k:m)
       v(k) = v(k) + sign(norm2(x(k:m)),x(k))

       ! Normalize the reflection vector
       v(k:m) = v(k:m)/norm2(v(k:m))

       ! Perform householder transformation to zero out the lower
       ! entries except one for each column of the matrix
       do j = k , n
          scalar = dot_product(v(k:m), A(k:m,j))
          A(k:m,j) = A(k:m,j) - 2.0_dp*scalar*v(k:m)
       end do

    end do cols

    deallocate(x,v)

  end subroutine householder
  
  !===================================================================!
  ! QR facorization using Householder. The upper triangular part of a
  ! contains R.
  !===================================================================!
  
  pure subroutine householder_factorization(A, Q, R)

    ! Arguments
    real(dp), intent(in)  :: A(:,:)
    real(dp), intent(inout) :: Q(:,:), R(:,:)

    ! Local variables
    real(dp), allocatable :: x(:), v(:), ek(:)
    integer  :: k, j, n, m
    real(dp) :: scalar

    allocate(x, source = A(:,1))
    allocate(v, source = A(:,1))

    R = A

    ! Initialize
    m = size(A,1)
    n = size(A,2)

    findR: do k = 1, n
       
       ! Extract the column below the diagonal
       x(k:m) = R(k:m, k)

       ! Fimd the reflection vector 
       v(k:m) = x(k:m)
       v(k) = v(k) + sign(norm2(x(k:m)),x(k))

       ! Normalize the reflection vector
       v(k:m) = v(k:m)/norm2(v(k:m))

       ! Perform householder transformation to zero out the lower
       ! entries except one for each column of the matrix
       do j = k , n
          scalar = dot_product(v(k:m), R(k:m,j))
          R(k:m,j) = R(k:m,j) - 2.0_dp*scalar*v(k:m)
       end do

    end do findR

    ! Form Q explicitly by dotting with identity vectors
    findQ: do k = n, 1, -1

       ! Extract the column below the diagonal
       x(k:m) = A(k:m, k)

       ! Fimd the reflection vector 
       v(k:m) = x(k:m)
       v(k) = v(k) + sign(norm2(x(k:m)),x(k))

       ! Normalize the reflection vector
       v(k:m) = v(k:m)/norm2(v(k:m))

       ! Perform householder transformation to zero out the lower
       ! entries except one for each column of the matrix
       Q(k,k) = 1.0d0
       do j = k, n
          scalar = dot_product(v(k:m), Q(k:m,j))
          Q(k:m,j) = Q(k:m,j) - 2.0_dp*scalar*v(k:m)
       end do

    end do findQ

    deallocate(x,v)

  end subroutine householder_factorization
  
  !===================================================================!
  ! Take a 2-entry vector, perform Givens rotation and return the new
  ! vector in place of the old vector
  !===================================================================!
  
  pure subroutine givens(x)

    ! Arguments
    real(dp), intent(inout)  :: x(2)

    ! Locals
    real(dp) :: c, s, d, h

    ! Prevent division by zero
    if (x(2) .ne. 0.0d0) then

       ! Find sin and cosine
       h = hypot(x(1),x(2))
       d = 1.0d0/h
       c = abs(x(1))*d
       s = sign(d,x(1))*x(2)

       ! Set the new first entry
       x(1) = sign(1.0_dp,x(1))*h

    end if

    ! Zero out the second entry
    x(2) = 0.0d0

  end subroutine givens

  pure subroutine givens_rotation(x,y, G)

    real(dp), intent(in)  :: x, y
    real(dp), intent(out) :: G(2,2)
    real(dp) :: s, c, h
    
    ! Prevent division by zero
    if (y .eq. 0.0d0) then

       g(1,1) = 1.0d0
       g(2,2) = 1.0d0
       g(1,2) = 0.0d0
       G(2,1) = 0.0d0

    else

       h = hypot(x,y)

       c =  x/h
       s = -y/h

       g(1,1) = c
       g(2,2) = c
       g(1,2) = -s
       G(2,1) = s

    end if

  end subroutine givens_rotation
  
  !===================================================================!
  ! QR facorization using Givens rotation.
  !===================================================================!
  
  subroutine givens_factorization(A, Q, R)

    ! Arguments
    real(dp), intent(in)    :: A(:,:)
    real(dp), intent(inout) :: Q(:,:), R(:,:)

    integer  :: k, i, n, m
    real(dp) :: G(2,2)

    ! Get matrix size
    m = size(A,1)
    n = size(A,2)

    ! Copy matrix
    R = A

    ! Initialize Q with Identity
    do concurrent(i=1:n)
       Q(i,i) = 1.0d0
    end do

    ! Left multiply A by G to get R
    ! Right multiply G to get Q^T
    cols: do k = 1, n
       rows: do i = m-1, k, -1

          ! Form 2 by 2 block of G matrix
          call givens_rotation(R(i,k),R(i+1,k), G) 

          ! Multiply only necessary rows and columns to form R and Q
          R(i:i+1,:) = matmul(G, R(i:i+1,:))
          Q(:,i:i+1) = matmul(Q(:,i:i+1),transpose(G))

       end do rows
    end do cols

  end subroutine givens_factorization

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
    integer :: i, j, n

    ! Initialize
    n = size(A,2)

    R = 0
    
    cols: do i = 1, n
       Q(:,i) = A(:,i)
       row: do j = 1, i-1
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
    integer :: i, j, n

    ! Initialize
    n = size(A,2)

    R = 0
    cols: do i = 1, n
       Q(:,i) = A(:,i)
       row: do j = 1, i-1
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
