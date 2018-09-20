program test

  ! Import dependencies
  use linear_algebra, only: solve
  use linear_transport, only: exact1, assemble_system
  use nonlinear_transport, only: exact2, exact3, exact4!, assemble_system2,  assemble_system3,  assemble_system4 

  implicit none

  call solve_linear_transport()
  !call solve_nonlinear_transport()

contains

subroutine solve_linear_transport()
  
  implicit none

  ! Problem setup
  integer, parameter :: npts = 200

  ! Matrices and vectors
  real(8), allocatable, dimension(:,:) :: A
  real(8), allocatable, dimension(:)   :: phi, b

  ! Physics parameters
  integer :: flag, i, j

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

  ! Print first order coeffs
  do i = 1, npts
     !write(*,"(100g15.5)") (A(i,j), j = 1, npts)
  enddo

  ! Solve the system
  phi = solve(A, b) 

  ! Write output
  do i = 1, npts
     write(11, *) dble(i)/dble(npts), phi(i), exact1(dble(i)/dble(npts)), &
          & exact2(dble(i)/dble(npts)), exact3(dble(i)/dble(npts)), exact4(dble(i)/dble(npts))
  end do

  close(11)

  deallocate(A,phi,b)

end subroutine solve_linear_transport

end program test
