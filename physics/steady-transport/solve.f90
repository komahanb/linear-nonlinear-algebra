program test

  ! Import dependencies
  use linear_algebra, only: solve

  implicit none

  ! Problem setup
  integer, parameter :: npts = 20

  ! Matrices and vectors
  real(8), allocatable, dimension(:,:) :: A
  real(8), allocatable, dimension(:)   :: phi, b

  ! Physics parameters
  integer :: flag, i, k

  ! Filename
  character(len=50) :: filename
  character(len=50) :: strnpts

  ! Create filename and open
  write(strnpts,*) int(npts)
  filename = 'case1-solution'//"-npts-" // trim(adjustl(strnpts)) //'.dat'
  open(11,file=filename)
  write(11,*) "phi ", "u"

  ! Allocate matrices
  allocate(A(npts, npts))
  A = 0.0d0

  ! Allocate vectors
  allocate(b(npts), phi(npts))
  b = 0.0d0
  phi = 0.0d0

  ! Assemble linear system
  call assemble_system(0.0d0, 1.0d0, npts, A, b)

  ! Solve the system
  phi = solve(A, b)

  ! Write output
  do i = 1, npts
     write(11, *) dble(i)/dble(npts), phi(i), exact1(dble(i)/dble(npts))
  end do

  close(11)

  deallocate(A,phi,b)

contains

  pure real(8) function exact1(x)

    real(8), intent(in) :: x
    real(8), parameter  :: L     = 1.0d0
    real(8), parameter  :: gamma = 0.1d0
    real(8), parameter  :: U     = 1.0d0

    exact1 = 1.0d0 - 1.0d0*(exp(x*U/gamma)-1.0d0)/(exp(L*U/gamma)-1.0d0)
    
  end function exact1
  
  ! Model problem to solve
  subroutine assemble_system(a, b, npts, V, rhs)

    real(8), intent(in)  :: a, b ! bound of domain
    real(8), intent(out) :: V(npts,npts)
    real(8), intent(out) :: rhs(npts)
    integer              :: npts ! number of points

    ! Local variables
    real(8) :: u, gamma
    real(8) :: aa, bb, cc
    real(8) :: h    
    integer :: m, n, i, j

    ! Mesh spacing
    h = (b-a)/dble(npts+1)

    ! Problem parameters
    u     = 1.0d0
    gamma = 0.1d0    
    aa    = -u/(2.0d0*h) - gamma/(h*h)
    bb    = 2.0d0*gamma/(h*h)
    cc    = u/(2.0d0*h) - gamma/(h*h)

    ! Prepare matrix
    V = 0.0d0    
    m = npts
    n = npts
    do i = 1, m
       do j = 1, n
          if (i .eq. j-1) then
             V(i,j) = aa
          else if (i .eq. j+1) then           
             V(i,j) = cc
          else if (i .eq. j) then           
             V(i,i) = bb
          else
             ! skip zeros
          end if
       end do
    end do
   
    ! Assemble the RHS
    rhs = 0.0d0
    rhs(1) = 1.0d0
    rhs(npts) = 0.0d0

  end subroutine assemble_system

end program test
