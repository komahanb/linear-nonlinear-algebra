program test
  
  ! Import dependencies
  use lapack, only : dgetrf
  use linear_algebra 
  use direct_linear_solve, only : dgetlu

  implicit none
  
  integer, parameter :: npts = 100
  real(8), parameter :: max_tol = 1.0d-6
  integer, parameter :: max_it  = 100000

  real(8) :: x(npts), xtmp(npts), b(npts), A(npts,npts)
  integer :: iter, flag, i, j
  real(8) :: tol, omega(20)
  real(8) :: scale = 1.0d0

  call assemble_system(0.0d0, 1.0d0, npts, A, b, x)
  xtmp = solve(A,b)

  print *, 'gauss seidel'    
  call assemble_system(0.0d0, 1.0d0, npts, A, b, x)
  call dseidel(A, b, 1.0d0, max_it, max_tol, x, iter, tol, flag)
  print *, 'seidel', tol, iter    
  open(11, file='seidel.dat')
  do i = 1, npts
     write(11, *) dble(i)/dble(npts), x(i), xtmp(i)
  end do
  close(11)

end program test

subroutine assemble_system(a, b, npts, V, rhs, x)

  use linear_algebra, only : inv

  real(8), intent(in)  :: a, b ! bound of domain
  integer              :: npts ! number of points
  real(8), intent(out) :: V(npts,npts)

  real(8), intent(out) :: rhs(npts)
  real(8), intent(out) :: x(npts)
  integer              :: m, n, i, j
  real(8) :: h
  real(8), parameter :: PI = 3.141592653589793d0

  h = (b-a)/dble(npts+1)
  V = 0.0d0

  m = npts
  n = npts
  do j = 1, m 
     do i = 1, n
        if (i .eq. j-1) then
           V(i,j) = -1.0d0
        else if (i .eq. j+1) then           
           V(i,j) = -1.0d0
        else if (i .eq. j) then           
           V(i,i) = 2.0d0 + 4.0d0*h*h
        else
           ! skip
        end if
     end do
  end do

  ! Assemble the RHS
  rhs = 0.0d0
  rhs(1) = -1.0d0
  rhs(npts) = 2.0d0

  ! Initial solution profile
  do i = 1, npts
     x(i) = - 1.0d0 + 2.75d0*(dble(i)*h)
  end do
  
end subroutine assemble_system
