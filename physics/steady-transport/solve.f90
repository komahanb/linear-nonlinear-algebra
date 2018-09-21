program test

  ! Import dependencies
  use linear_algebra, only: solve
  use linear_transport, only: exact1, assemble_system
  use nonlinear_transport, only: exact2, exact3, exact4!, assemble_system2,  assemble_system3,  assemble_system4 
  use direct_linear_solve, only: thomas

  implicit none

  call solve_linear_transport()
  
  ! call solve_nonlinear_transport()
  ! call test_thomas(8)

contains
  
  subroutine solve_linear_transport()

    implicit none

    ! Problem setup
    integer, parameter :: npts = 100
    logical, parameter :: sparse = .false.

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
    !filename = 'solution'//"-npts-" // trim(adjustl(strnpts)) //'.dat'
    filename = 'solution.dat'
    
    open(11, file=filename)
    write(11, *) "phi ", "u"

    if (sparse .eqv. .true.) then    
       allocate(A(npts, 3))  
    else
       allocate(A(npts, npts))  
    end if
    allocate(b(npts), phi(npts))
    call assemble_system(0.0d0, 1.0d0, npts, A, b, sparse)
    
    do i = 1, npts
       write(*,'(10f10.4)') (A(i,j), j = 1, 3), b(i)
    end do
    
    ! Solve the system
    if (sparse .eqv. .true.) then
       call thomas(A, b)
       phi = b
    else
       phi = solve(A, b)
    end if

    ! Write output
    do i = 1, npts
       write(11, *) dble(i)/dble(npts), phi(i), &
            & exact1(dble(i)/dble(npts)), &
            & exact2(dble(i)/dble(npts)), &
            & exact3(dble(i)/dble(npts)), &
            & exact4(dble(i)/dble(npts))
    end do

    close(11)

    deallocate(A,phi,b)

  end subroutine solve_linear_transport
  
  subroutine test_thomas(n)

    implicit none

    ! Arguments
    integer, intent(in) :: n 

    ! Local variables
    integer, parameter :: bandwidth = 3
    real(8)            :: A(n, bandwidth), b(n), x(n)
    integer            :: i, j

    ! Setup the matrix
    A(1,:) = [0.0d0, 4.0d0, -1.0d0]
    do concurrent(i=2:n-1)
       A(i,:) = [-1.0d0, 4.0d0, -1.0d0]
    end do
    A(n,:) = [-1.0d0, 4.0d0, 0.0d0]
    
    ! Setup right hand side
    b = 0.0d0
    b(n) = 16.0d0

    do i = 1, n
       write(*,'(10f10.6)') (A(i,j), j = 1, bandwidth), b(i)
    end do
    
    ! Solve the linear system Ax = b
    call thomas(A, b)

    print *, "solution is"
    do i = 1, n
       write(*,'(f10.6)') b(i)
    end do

  end subroutine test_thomas

end program
