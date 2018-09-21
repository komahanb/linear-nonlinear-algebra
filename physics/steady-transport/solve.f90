program test

  ! Import solvers
  use linear_algebra, only: solve
  use direct_linear_solve, only: tdma

  ! Import physics
  use linear_transport, only: exact1, assemble_system
  use nonlinear_transport, only: exact2, exact3, exact4, assemble_residual_jacobian

  implicit none

  ! Problem setup
  logical, parameter :: sparse = .false.
  integer, parameter :: sizes(10) = [10,20,30,40,50,60, 70, 80, 90, 100]

  ! Solve the linear transport
  call solve_linear_transport(sparse, sizes, exact1)

  ! Solve all three nonlinear transport cases
  call solve_nonlinear_transport(sparse, sizes, 1, exact2)
  call solve_nonlinear_transport(sparse, sizes, 2, exact3)
  call solve_nonlinear_transport(sparse, sizes, 3, exact4)

contains

  !===================================================================!
  ! Setup of nonlinear transport equation
  !===================================================================!
  
  subroutine solve_nonlinear_transport(sparse, sizes, case, exact)
    
    interface
       pure real(8) function exact(x)
         real(8), intent(in) ::  x
       end function exact
    end interface

    ! Problem setup
    logical, intent(in) :: sparse
    integer, intent(in) :: sizes(:)
    integer, intent(in) :: case

    real(8), parameter :: tau_r = 1.0d-12
    real(8), parameter :: tau_a = 1.0d-9
    integer, parameter :: max_it = 100

    ! Matrices and vectors
    real(8), allocatable, dimension(:,:) :: A
    real(8), allocatable, dimension(:)   :: phi, Q, R

    integer :: i, j
    real(8) :: xi, rmse, walltime, time_start, time_end

    integer :: npts
    character(len=100) :: filename
    character(len=100) :: strnpts
    character(len=100) :: strcase

    ! Open file  
    open(12, file='nonlinear_summary.dat')
    write(12, *) "npts ", "h ", "rmse ", "wall_time"

    do j = 1, size(sizes)

       npts = sizes(j)
       
       ! Create filename and open
       write(strnpts,*) int(npts)
       write(strcase,*) int(case)
       filename = 'nonlinear-solution' // '-npts-' // trim(adjustl(strnpts)) // "-case-"// trim(adjustl(strcase)) //'.dat'
       open(11, file=filename)
       write(11, *) "x ", "phihat ", "phi"

       call cpu_time(time_start)

       allocate(R(npts))
       allocate(A(npts, npts))

       ! Initial solution
       allocate(phi(npts))  
       do i = 1, npts
          xi = dble(i)/dble(npts+1)
          phi(i) = exact(xi)
       end do

       ! Set the source term
       allocate(Q(npts))
       if (case.eq.1) then
          Q = 0.0d0
       else if (case.eq.2) then 
          Q = 0.1d0
       else 
          do i = 1, npts
             xi = dble(i)/dble(npts+1)
             Q(i) = 0.1d0*xi
          end do
       end if

       ! Solve the nonlinear system
       call newton(sparse, tau_r, tau_a, max_it, npts, phi, Q, A, R)
       call cpu_time(time_end)
       walltime = time_end - time_start

       ! Compare solution to exact and write output
       rmse = 0.0d0
       do i = 1, npts
          xi = dble(i)/dble(npts+1)
          write(11, *) xi, phi(i), exact(xi)
          rmse = rmse + (exact(xi)-phi(i))**2.0d0
       end do
       rmse = rmse/sqrt(dble(npts))
       write(12, *) npts, 1.0d0/dble(npts+1), rmse, walltime

       close(11)

       deallocate(R,A,phi,Q)

    end do

    close(12)

  end subroutine solve_nonlinear_transport
  
  !===================================================================!
  ! Newton method for solving nonlinear system that features quadratic
  ! convergence
  !===================================================================!
  
  subroutine newton(sparse, tau_r, tau_a, max_it, npts, phi, Q, jac, res )

    logical, intent(in)    :: sparse
    real(8), intent(in)    :: tau_r, tau_a
    integer, intent(in)    :: max_it
    integer, intent(in)    :: npts
    real(8), intent(inout) :: phi(:)
    real(8), intent(in)    :: Q(:)
    real(8), intent(inout) :: res(:)
    real(8), intent(inout) :: jac(:,:)
    
    ! local variables
    real(8) :: r0 
    real(8) :: s(size(res))
    integer :: iter
  
    call assemble_residual_jacobian(sparse, npts, phi, Q, jac, res)

    ! Initial residual
    r0 = norm2(res)

    iter = 0
    do while (norm2(res) > tau_r*r0 + tau_a .and. iter .le. max_it)

       ! Increment the iteration count
       iter = iter + 1

       ! Residual and jacobian
       call assemble_residual_jacobian(sparse, npts, phi, Q, jac, res)
       
       ! Solve the linear system based on matrix storage format
       if (sparse .eqv. .true.) then
          call tdma(jac, res, s)
       else
          s = solve(jac, res)
       end if

       ! Apply the update
       phi = phi + s

       ! print newton solve details
       print *, iter, norm2(res)

    end do

  end subroutine newton

  !===================================================================!
  ! Setup of linear transport equation
  !===================================================================!

  subroutine solve_linear_transport(sparse, sizes, exact)

    interface
       pure real(8) function exact(x)
         real(8), intent(in) ::  x
       end function exact
    end interface

    ! Problem setup
    logical, intent(in) :: sparse
    integer, intent(in) :: sizes(:)
  
    ! Matrices and vectors
    real(8), allocatable, dimension(:,:) :: A
    real(8), allocatable, dimension(:)   :: phi, b

    real(8) :: xi, rmse, walltime, time_start, time_end

    ! Physics parameters
    integer :: i, j
    integer :: npts

    ! Filename
    character(len=100) :: filename
    character(len=100) :: strnpts

    ! Open file  
    open(12, file='linear_summary.dat')
    write(12, *) "npts ", "h ", "rmse ", "wall_time"

    do j = 1, size(sizes)

       npts = sizes(j)

       ! Create filename and open
       write(strnpts,*) int(npts)
       filename = 'linear-solution' // '-npts-' // trim(adjustl(strnpts)) // '.dat'    
       open(11, file=filename)
       write(11, *) "x ", "phihat ", "phi"

       call cpu_time(time_start)

       if (sparse .eqv. .true.) then    
          allocate(A(npts, 3))  
       else
          allocate(A(npts, npts))  
       end if
       allocate(b(npts), phi(npts))
       call assemble_system(0.0d0, 1.0d0, npts, A, b, sparse)

       ! Solve the system
       if (sparse .eqv. .true.) then
          call tdma(A, b, phi)
       else
          phi = solve(A, b)
       end if

       call cpu_time(time_end)
       walltime = time_end - time_start

       ! Write output
       rmse = 0.0d0
       do i = 1, npts
          xi = dble(i)/dble(npts+1)
          write(11, *) xi, phi(i), exact(xi)
          rmse = rmse + (exact(xi)-phi(i))**2.0d0
       end do
       rmse = rmse/sqrt(dble(npts))
       write(12, *) npts, 1.0d0/dble(npts+1), rmse, walltime

       ! Free resources
       close(11)
       deallocate(A, phi, b)

    end do

    close(12)

  end subroutine solve_linear_transport

end program
