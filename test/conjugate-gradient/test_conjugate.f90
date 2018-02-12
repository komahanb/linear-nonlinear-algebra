!=====================================================================!
! Run gradient-based descent algorithms to solve linear algebra
!=====================================================================!

program test

  implicit none

  call check_conjugate

end program test

!---------------------------------------------------------------------!
! 
!---------------------------------------------------------------------!

subroutine check_conjugate

  use linear_algebra

  implicit none

  dirichlet : block

    integer, parameter :: npts = 64
    real(8), parameter :: max_tol = 1.0d-6
    integer, parameter :: max_it = 100000
    real(8) :: x(npts), xtmp(npts), b(npts), A(npts,npts), P(npts, npts)
    integer :: iter, flag, i, j
    real(8) :: tol

    ! Solve using LU factorization
    call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P)
    xtmp = solve(A,b)

    ! solve using CG
    call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P) 
    call dcg(A, b, max_it, max_tol, x, iter, tol, flag)
    print *, 'cg', tol, iter
    
    ! Solve using preconditioned CG
    call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P) 
    call dpcg(A, P,  b, max_it, max_tol, x, iter, tol, flag)
    print *, 'pcg', tol, iter 

    open(11, file='dirichlet.dat')
    do i = 1, npts
       write(11, *) dble(i)/dble(npts+1), x(i), xtmp(i)
    end do
    close(11)

  end block dirichlet
  
  mixed : block
    
    integer, parameter :: npts = 64
    real(8), parameter :: max_tol = 1.0d-6
    integer, parameter :: max_it = 100000
    real(8) :: x(npts+1), xtmp(npts+1), b(npts+1), A(npts+1,npts+1), P(npts+1, npts+1)
    integer :: iter, flag, i, j
    real(8) :: tol

    ! Solve using LU factorization
    call assemble_system_mixed(0.0d0, 1.0d0, npts, A, b, x, P)
    xtmp = solve(A,b)

    ! Solve using CG
    call assemble_system_mixed(0.0d0, 1.0d0, npts, A, b, x, P) 
    call dcg(A, b, max_it, max_tol, x, iter, tol, flag)
    print *, 'cg', tol, iter

    ! Solve using preconditioned CG
    call assemble_system_mixed(0.0d0, 1.0d0, npts, A, b, x, P) 
    call dpcg(A, P,  b, max_it, max_tol, x, iter, tol, flag)
    print *, 'pcg', tol, iter
    
    open(11, file='mixed.dat')
    do i = 1, npts + 1
       write(11, *) dble(i)/dble(npts+1), x(i), xtmp(i)
    end do
    close(11)

  end block mixed

end subroutine check_conjugate

!---------------------------------------------------------------------!
! Assemble -U_xx = 2x - 0.5, U(0) = U(1)= 0, x in [0,1]
!---------------------------------------------------------------------!

subroutine assemble_system_dirichlet(a, b, npts, V, rhs, u, P)

  implicit none
  
  real(8), intent(in)  :: a, b ! bounds of the domain
  integer, intent(in)  :: npts ! number of interior points
  real(8), intent(out) :: V(npts,npts) ! banded matrix
  real(8), intent(out) :: rhs(npts)
  real(8), intent(out) :: u(npts)
  real(8), intent(out) :: P(npts,npts)
  real(8) :: PTMP(npts,npts)
  
  real(8), parameter :: PI = 3.141592653589793d0
  real(8)            :: h, alpha
  integer            :: M, N
  integer            :: i, j, k
  
  ! h = width / num_interior_pts + 1
  h = (b-a)/dble(npts+1)
  V = 0.0d0

  ! Size of the linear system = unknowns (interior nodes)
  M = npts ! nrows
  N = npts ! ncols
  
  ! Set the inner block
  rows: do i = 1, M
     cols: do j = 1, N
        if (j .eq. i-1) then
           ! lower triangle
           V(j,i) = -1.0d0
        else if (j .eq. i+1) then           
           ! upper triangle
           V(j,i) = -1.0d0
        else if (j .eq. i) then           
           ! diagonal
           V(j,i) = 2.0d0
        else
           ! skip
        end if
     end do cols
  end do rows
  
  ! Assemble the RHS
  do i = 1, M
     rhs(i) = h*h*(2.0d0*dble(i)*h - 0.5d0)
  end do
  rhs(1) = rhs(1) + 1.0d0
  rhs(M) = rhs(M) + 1.0d0
  
  ! Initial solution profile use sin function as a first guess
  do i = 1, M
     u(i) =  sin(dble(i)*h*PI)
  end do

  ! Set a preconditioner
  !P = inv(V)
  alpha = sqrt(2.0d0/dble(npts+1))
  do j = 1, M
     do k = 1, N
        PTMP(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
     end do
  end do

  p = matmul(ptmp, ptmp)

end subroutine assemble_system_dirichlet

!---------------------------------------------------------------------!
! Assemble -U_xx = 2x - 0.5, U(0) = U_x(1)= 0, x in [0,1]
!---------------------------------------------------------------------!

subroutine assemble_system_mixed(a, b, npts, V, rhs, u, P)

  implicit none
  
  real(8), intent(in)  :: a, b ! bounds of the domain
  integer, intent(in)  :: npts ! number of interior points
  real(8), intent(out) :: V(npts+1,npts+1) ! banded matrix
  real(8), intent(out) :: rhs(npts+1)
  real(8), intent(out) :: u(npts+1)
  real(8), intent(out) :: P(npts+1,npts+1)
  real(8) :: ptmp(npts+1, npts+1)  
  real(8), parameter :: PI = 3.141592653589793d0
  real(8)            :: h, alpha
  integer            :: M, N
  integer            :: i, j, k
  
  ! h = width / num_interior_pts + 1
  h = (b-a)/dble(npts+1)
  V = 0.0d0

  ! Size of the linear system = unknowns
  M = npts + 1 ! nrows
  N = npts + 1 ! ncols
  
  ! Set the inner block
  rows: do i = 1, M
     cols: do j = 1, npts
        if (j .eq. i-1) then
           ! lower triangle
           V(j,i) = -1.0d0
        else if (j .eq. i+1) then           
           ! upper triangle
           V(j,i) = -1.0d0
        else if (j .eq. i) then           
           ! diagonal
           V(j,i) = 2.0d0
        else
           ! skip
        end if
     end do cols
  end do rows
  V(M,M) = 1.0d0
  V(M,N-1) = -1.0d0

  ! Assemble the RHS
  do i = 1, M
     rhs(i) = h*h*(2.0d0*dble(i)*h - 0.5d0)
  end do
  rhs(1) = rhs(1) + 1.0d0
  rhs(M) = 0.0d0

  ! Initial solution profile use sin function as a first guess
  do i = 1, M
     u(i) = sin(dble(i)*h*PI)
  end do

  alpha = sqrt(2.0d0/dble(npts+1))
  do j = 1, M
     do k = 1, N
        PTMP(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
     end do
  end do

  p = matmul(ptmp, ptmp)

  ! Set a preconditioner
  !P = inv(V)
!!$  alpha = sqrt(2.0d0/dble(npts+1))
!!$  do j = 1, npts
!!$     do k = 1, npts
!!$        P(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
!!$     end do
!!$  end do

end subroutine assemble_system_mixed
