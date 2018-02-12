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
  
  integer, parameter :: npts = 64
  real(8), parameter :: max_tol = 1.0d-6
  integer, parameter :: max_it = 100000
  real(8) :: x(npts), xtmp(npts), b(npts), A(npts,npts), P(npts, npts)
  integer :: iter, flag, i, j
  real(8) :: tol

  call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P)
  xtmp = solve(A,b)

  call assemble_system_dirichlet(0.0d0, 1.0d0, npts, A, b, x, P)
  call dsd(A, b, max_it, max_tol, x, iter, tol, flag)
  !call dcg(A, b, max_it, max_tol, x, iter, tol, flag)
  !call dpcg(A, P,  b, max_it, max_tol, x, iter, tol, flag)

  print *, 'sd', tol, iter 

  open(11, file='checkcg.dat')
  do i = 1, npts
     write(11, *) dble(i-1)/dble(npts), x(i), xtmp(i)
  end do
  close(11)

end subroutine check_conjugate

subroutine assemble_system_mixed(a, b, npts, V, rhs, x, P)

  implicit none
  
  real(8), intent(in)  :: a, b         ! bound of domain
  integer              :: npts         ! number of total points
  real(8), intent(out) :: V(npts,npts) ! banded matrix
  real(8), intent(out) :: rhs(npts)
  real(8), intent(out) :: x(npts)
  real(8), intent(out) :: P(npts,npts)
  integer              :: i, j, k
  real(8)              :: h, alpha
  real(8), parameter   :: PI = 3.141592653589793d0

  ! h = width / num_interior_pts + 1
  h = (b-a)/dble(npts)
  V = 0.0d0

  ! Set the inner block
  do j = 1, npts - 1
     do i = 1, npts - 1
        if (i .eq. j-1) then
           ! lower triangle
           V(i,j) = -1.0d0
        else if (i .eq. j+1) then           
           ! upper triangle
           V(i,j) = -1.0d0
        else if (i .eq. j) then           
           ! diagonal
           V(i,i) = 2.0d0
        else
           ! skip
        end if
     end do
  end do

  ! Set the first and last rows
  V(npts, npts) = 1.0d0
  V(npts, npts-1) = -1.0d0
  
  ! Assemble the RHS
  do i = 1, npts - 1
     rhs(i) = h*h*(2.0d0*dble(i)*h - 0.5d0)
  end do
  rhs(1) = rhs(1) + 1.0d0
  rhs(npts) = 0.0d0


  ! Initial solution profile
  do i = 1, npts
     !x(i) =  1.0d0 - 1.1d0*(dble(i)*h)
  end do
  
  ! Set a preconditioner
  !P = inv(V)
  alpha = sqrt(2.0d0/dble(npts+1))
  do j = 1, npts
     do k = 1, npts
        P(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
     end do
  end do

end subroutine assemble_system_mixed

!---------------------------------------------------------------------!
! Assemble -U_xx = 2x - 0.5, U(0) = U(1)= 0, x in [0,1]
!---------------------------------------------------------------------!

subroutine assemble_system_dirichlet(a, b, npts, V, rhs, u, P)

  implicit none
  
  real(8), intent(in)  :: a, b ! bounds of the domain
  integer              :: npts ! number of interior points
  real(8), intent(out) :: V(npts,npts) ! banded matrix
  real(8), intent(out) :: rhs(npts)
  real(8), intent(out) :: u(npts)
  real(8), intent(out) :: P(npts,npts)
  
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
     print *, dble(i)*h, u(i)
  end do

  ! Set a preconditioner
  !P = inv(V)
!!$  alpha = sqrt(2.0d0/dble(npts+1))
!!$  do j = 1, npts
!!$     do k = 1, npts
!!$        P(k,j) = alpha*sin(PI*dble(j*k)/dble(npts+1))
!!$     end do
!!$  end do

end subroutine assemble_system_dirichlet
