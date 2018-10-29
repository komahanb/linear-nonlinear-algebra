!=====================================================================!
! 
!=====================================================================!

program test

  use direct_linear_solve, only : mgs_qrfactor, cgs_qrfactor, &
       & householder, banded_householder, householder_factorization
  use linear_algebra, only : eigvals, svdvals

  implicit none

  integer, parameter :: n  = 20

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

contains
  
  subroutine toepliz(col, A)

    real(8), allocatable :: A(:,:)
    real(8) :: col(:)
    integer :: m, i, j

    m = size(A, 1)
    do j = 1, m
       do i = 1, m
          if (i.eq.j) then
             A(i,j) = col(1)
          else if (i+1.eq.j .or.i-1.eq.j) then
             A(i,j) = col(2)
             A(j,i) = col(2)
          else if (i+2.eq.j .or.i-2.eq.j) then
             A(i,j) = col(3)
             A(j,i) = col(3)
          else if (i+3.eq.j .or.i-3.eq.j) then
             A(i,j) = col(4)
             A(j,i) = col(4)
          else if (i+4.eq.j .or.i-4.eq.j) then
             A(i,j) = col(5)
             A(j,i) = col(5)
          end if
       end do
    end do

  end subroutine toepliz

end program test
