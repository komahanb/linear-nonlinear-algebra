program main

  use iso_fortran_env     , only : dp => REAL64
  use matrix_utils        , only : toepliz, print
  use direct_linear_solve , only : mgs_qrfactor, cgs_qrfactor, &
       & householder, banded_householder, householder_factorization
  use direct_linear_solve , only : givens_factorization, givens, qriteration
  use linear_algebra      , only : eigvals, svdvals
  
  implicit none

  real(dp), parameter :: lambda = 2.0d0  
  integer , parameter :: npts = 5
  real(8) , allocatable, dimension(:,:) :: A, Q, R, D, U

  integer  :: i
  real(dp) :: h

  allocate(A(npts,npts))
  allocate(Q(npts,npts))
  allocate(R(npts,npts))
  allocate(D(npts,npts))
  allocate(U(npts,npts))

  call assemble_matrix(0.0d0, 1.0d0, npts, lambda, A)
  call qriteration(A, R, Q, 100, 1.0d-3)

  stop
  
  print *, 'Performing Givens Transformation'
  Q = 0
  R = 0
  call givens_factorization(A, Q, R)
  print *, 'matrix'
  call print(A)

  print *, 'Q'
  call print(Q)
  
  print *, 'R'
  call print(R)

  ! Check similarity with A
  print *, 'consistency A' 
  D = matmul(Q, R)
  print *, norm2(A - D)

  ! Check Orthogonality
  D = matmul(Q,transpose(Q))
  print *, 'orthogonality Q'
  do i = 1, npts       
     D(i,i) = 1.0d0 - D(i,i)
  end do
  print *, norm2(D)

  deallocate(A,Q,R,D)

contains
  
  ! Model problem to solve
  subroutine assemble_matrix(a, b, npts, lambda, V)

    real(8), intent(in)  :: a, b ! bound of domain
    integer              :: npts ! number of points
    real(8), intent(out) :: V(npts,npts)
    real(8), intent(in)  :: lambda
    real(8), parameter   :: PI = 3.141592653589793d0
    integer              :: m, n, i, j

    h = (b-a)/dble(npts+1)
    V = 0.0d0

    m = npts
    n = npts
    do i = 1, m
       do j = 1, n
          if (i .eq. j-1) then
             V(i,j) = -1.0d0
          else if (i .eq. j+1) then           
             V(i,j) = -1.0d0
          else if (i .eq. j) then           
             V(i,i) = 2.0d0 + lambda*h*h
          else
             ! skip zeros
          end if
       end do
    end do

  end subroutine assemble_matrix

end program main
