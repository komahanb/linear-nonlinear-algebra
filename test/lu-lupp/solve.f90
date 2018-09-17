program test

  ! Import dependencies
  use lapack, only : dgetrf
  use linear_algebra 
  use direct_linear_solve, only : dlufactor, dluppfactor

  implicit none

  integer, parameter :: npts = 4000
  real(8) :: x(npts), xtmp(npts), b(npts), y(npts)
  integer :: iter, flag, i, j
  real(8) :: tol, omega(20)
  real(8) :: scale = 1.0d0
  real(8) :: lambda = 1.0d0

  real(8), allocatable :: A(:,:)
  real(8), allocatable :: L(:,:)
  real(8), allocatable :: U(:,:)
  real(8), allocatable :: P(:,:)

  allocate(A(npts, npts))
  allocate(L, U, P, source = A)
  
  real :: start, finish
  call cpu_time(start)

  call assemble_system(0.0d0, 1.0d0, npts, A, b, x, lambda)
!!$  xtmp = solve(A,b)
  print *, "assembled system"
  
  ! Initialize matrices
  L = 0.0d0
  U = 0.0d0
  P = 0.0d0

  ! Ax = b  is now LUx=b
  ! call dlufactor(A, L, U, flag)
  
  call dluppfactor(A, L, U, P, flag)
  b = matmul(P,b)
  print *, "factorized matrix"
  
  ! Solve Ly = b
  call fwdsub(L, b, y, flag)

  ! Solve Ux = y
  call backsub(U, y, x, flag)

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start

  ! Write the time for each linear solve
  open(11, file='time.dat')
  write(11, *) "size", "time"
  write(11, *) npts, finish-start
  close(11)

  ! Write solution
  open(11, file='lu.dat')
  write(11, *) "x ", "lu ", "exact"
  do i = 1, npts
     write(11, *) dble(i)/dble(npts), x(i), xtmp(i)
  end do
  close(11)

contains


  
  ! Model problem to solve
  subroutine assemble_system(a, b, npts, V, rhs, u, lambda)

    use linear_algebra, only : inv

    real(8), intent(in)  :: a, b ! bound of domain
    integer              :: npts ! number of points
    real(8), intent(out) :: V(npts,npts)
    real(8), intent(out) :: rhs(npts)
    real(8), intent(out) :: u(npts)
    real(8)              :: h
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

    ! Assemble the RHS
    rhs = 0.0d0
    rhs(1)    = 1.0d0 + h*h*(3.0d0*(dble(i)*h)**2.0d0 - 0.5d0)
    rhs(npts) = 0.0d0 + h*h*(3.0d0*(dble(npts)*h)**2.0d0 - 0.5d0)

    ! Initial solution profile (straight line)
    do i = 1, npts
       u(i) = 1.0d0 + ((0.0d0+1.0d0)/(b-a))*(dble(i)*h)
    end do

  end subroutine assemble_system

end program test
