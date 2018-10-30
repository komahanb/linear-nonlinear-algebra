!=====================================================================!
! 
!=====================================================================!

program test

  use iso_fortran_env     , only : dp => REAL64
  use matrix_utils        , only : toepliz
  use direct_linear_solve , only : mgs_qrfactor, cgs_qrfactor, &
       & householder, banded_householder, householder_factorization
  !use direct_linear_solve , only : givens_factorization
  use linear_algebra      , only : eigvals, svdvals

  implicit none

  call test_givens_factorization()

contains
  
  subroutine test_givens_factorization

    real(dp), allocatable :: A(:,:), Q(:,:), R(:,:), D(:,:)
    integer , parameter :: matsize = 4
    integer :: i

    allocate(Q(matsize, matsize))
    allocate(R(matsize, matsize))
    allocate(D(matsize, matsize))

    allocate(A(matsize, matsize))
    call matrix(A)
    print *, 'matrix'
    do i = 1, matsize
       print *, A(i,:)
    end do

    print *, 'Performing Givens Transformation'
    Q = 0
    R = 0
    !call givens_factorization(A, Q, R)

    print *, 'Upper Triangular Matrix'
    do i = 1, matsize
       print *, R(i,:)
    end do
    
    ! Check Orthogonality
    D = matmul(Q,transpose(Q))
    print *, 'orthogonality'
    do i = 1, matsize
       D(i,i) = 1.0d0 - D(i,i)
       print *, D(i,i)
    end do

    deallocate(A,Q,R,D)

  end subroutine test_givens_factorization

  !===================================================================!
  ! Setup Toepliz matrix with supplied columns
  !===================================================================!
  
  subroutine matrix(A)

    real(8), allocatable :: A(:,:)
    real(8), allocatable :: diag(:)
    integer :: m, j

    A = 0

    call random_number(A(:,1))
    call random_number(A(1,:))

    m = size(A, 1)

    allocate(diag(m))
    call random_number(diag)    
    do concurrent(j=1:m)
       A(j,j) = diag(j)
    end do
    deallocate(diag)

  end subroutine matrix
  
  subroutine test_cgs_mgs_householder

    integer, parameter :: n  = 200

    ! integer, parameter :: sizes(10) = [100,200,300,400,500,600,700,800,900,1000]
    ! integer, parameter :: sizes(7) = [100,200,400,800,1600,3200,6400]

    real(8), allocatable :: A(:,:), B(:,:)
    real(8), allocatable :: Q(:,:), R(:,:)
    real(8), allocatable :: U(:,:), V(:,:), S(:,:), RC(:,:), RM(:,:)
    real(8), allocatable :: svda(:)

    real(8) :: walltime1, walltime2, time_end, time_Start
    integer :: i, info

    ! Allocte matrices
    allocate(A(n,n), B(n,n), Q(n,n), R(n,n))
    allocate(U(n,n), V(n,n), S(n,n), RC(n,n), RM(n,n))
    U = 0
    V = 0
    S = 0
    allocate(svda(n))

    ! Form a matrix with required singular values
    call random_number(A)
    call mgs_qrfactor(A, U, R, info)

    call random_number(A)
    call mgs_qrfactor(A, V, R, info)
    do i = 1, n
       S(i,i) = 3.0d0**dble(-i)
    end do
    A = matmul(U,matmul(S,V))

    print *, 'performing classical GS on A'
    call cgs_qrfactor(A, Q, RC, info)
    S = matmul(Q,transpose(Q))
    do i = 1, n
       S(i,i) = 1.0d0 - S(i,i)
       print *, S(i,i)
    end do

    ! Perform orthogonality check
    print *, n, norm2(S), norm2(A-matmul(Q,RC))

    print *, 'performing modified GS on A'
    call mgs_qrfactor(A, Q, RM, info)
    S = matmul(Q,transpose(Q))
    do i = 1, n
       S(i,i) = 1.0d0 - S(i,i)
       print *, S(i,i)
    end do

    ! Perform orthogonality check
    print *, n, norm2(S), norm2(A-matmul(Q,RM))

    print *, 'Householder'
    Q = 0
    R = 0
    !call householder(A)
    call householder_factorization(A, Q, R)
    S = matmul(Q,transpose(Q))
    do i = 1, n
       S(i,i) = 1.0d0 - S(i,i)
       print *, S(i,i)
    end do

    print *, "index ", "cgs ", "mgs ", "householder"
    do i = 1, n
       print *, i, RC(i,i), RM(i,i), abs(R(i,i))
    end do

    stop

    ! system 2
!!$
!!$  print *, "npts ", "walltime1 ", "walltime2"
!!$  
!!$  do i = 1, size(sizes)
!!$
!!$     n = sizes(i)
!!$
!!$     allocate(svda(n))
!!$     allocate(A(n,n), Q(n,n), R(n,n))
!!$     svda = 0
!!$     A = 0
!!$     Q = 0
!!$     R = 0
!!$
!!$     svda(1) = 10.0d0
!!$     svda(2) = 3.0d0
!!$     svda(3) = 2.0d0
!!$     svda(4) = 1.0d0
!!$
!!$     call toepliz(svda, A)
!!$
!!$     call cpu_time(time_start)              
!!$     call householder(A)
!!$     call cpu_time(time_end)
!!$     walltime1 = time_end - time_start
!!$
!!$     A = 0  
!!$     call toepliz(svda, A)
!!$
!!$     call cpu_time(time_start)              
!!$     call banded_householder(A,3,3)
!!$     call cpu_time(time_end)
!!$     walltime2 = time_end - time_start
!!$
!!$     print *, n, walltime1, walltime2
!!$
!!$     
!!$     deallocate(svda, A, Q, R)
!!$     
!!$  end do
!!$
!!$  !do i = 1, n
!!$  !   write(*,*) i, A(i,:)
!!$  !end do
!!$
!!$  stop
!!$  print *, 'performing classical GS on A'
!!$
!!$  !call mgs_qrfactor(A, Q, R, info)
!!$  !call banded_householder(A, 2, 2)
!!$
!!$  !print *, A- matmul(Q,R)
!!$  do i = 1, n
!!$     print *, i, A(i,i)
!!$  end do
!!$  

  end subroutine test_cgs_mgs_householder

end program test
