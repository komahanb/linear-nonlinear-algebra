program test

  ! Import dependencies
  use lapack, only : dgetrf
  use linear_algebra 
  use direct_linear_solve, only : dlufactor, dluppfactor
  
  implicit none

  integer, parameter :: sizes(6) = [200, 400, 800, 1000, 2000, 4000]
  integer :: i, npts

  ! Run for all problem sizes
  do i = 1, size(sizes)

     npts = sizes(i)

     call gesolve(npts)
     call geppsolve(npts)
     !call exact(npts)

  end do

!  call time_report(.false.)
!  call time_report(.true.)

contains

  subroutine geppsolve(npts)

    ! Problem setup
    integer, intent(in) :: npts

    ! Matrices and vectors
    real(8), allocatable, dimension(:,:) :: A, L, U, P
    real(8), allocatable, dimension(:)   :: x, xtmp, b, y

    ! Physics parameters
    real(8), parameter :: lambdas(2) = [0.0d0, 10.0d0]
    real(8) :: lambda
    integer :: flag, i, k
     
    ! Filename
    character(len=50) :: filename
    character(len=50) :: strlambda, strnpts

    ! Run for each lambda
    do k = 1, size(lambdas)

       lambda = lambdas(k)

       ! Create filename and open
       write(strlambda,*) int(lambda)
       write(strnpts,*) int(npts)
       filename = 'gepp-solution-lam-'// trim(adjustl(strlambda)) &
            & //"-npts-" // trim(adjustl(strnpts)) //'.dat'
       open(11,file=filename)
       write(11,*) "x ", "u"

       ! Allocate matrices
       allocate(A(npts, npts))
       allocate(L, U, P, source = A)
       A = 0.0d0
       L = 0.0d0
       U = 0.0d0
       P = 0.0d0

       ! Allocate vectors
       allocate(x(npts))
       allocate(xtmp, b, y, source = x)
       x    = 0.0d0
       xtmp = 0.0d0 
       b    = 0.0d0

       ! Solve x
       call assemble_system(0.0d0, 1.0d0, npts, A, b, x, lambda)
       call dluppfactor(A, L, U, P, flag)
       b = matmul(P,b)
       call fwdsub(L, b, y, flag)
       call backsub(U, y, x, flag)

       ! Write output
       do i = 1, npts
          write(11, *) dble(i)/dble(npts), x(i)
       end do

       close(11)

       deallocate(A,L,U,P,x,xtmp,b,y)

    end do

  end subroutine geppsolve

  pure real(8) function uhat(x, lambda)

    real(8), intent(in) :: x
    real(8), intent(in) :: lambda

    real(8) :: alpha

    alpha = sqrt(lambda)
    uhat =  2.d0*exp(-alpha*x)*(exp(2.d0*alpha)-exp(2.d0*alpha*x))/(exp(2.d0*alpha)-1.d0)

  end function uhat

  subroutine gesolve(npts)

    ! Problem setup
    integer, intent(in) :: npts

    ! Matrices and vectors
    real(8), allocatable, dimension(:,:) :: A, L, U, P
    real(8), allocatable, dimension(:)   :: x, xtmp, b, y

    ! Physics parameters
    real(8), parameter :: lambdas(2) = [0.0d0, 10.0d0]
    real(8) :: lambda
    integer :: flag, i, k
     
    ! Filename
    character(len=50) :: filename
    character(len=50) :: strlambda, strnpts

    ! Run for each lambda
    do k = 1, size(lambdas)

       lambda = lambdas(k)

       ! Create filename and open
       write(strlambda,*) int(lambda)
       write(strnpts,*) int(npts)
       filename = 'ge-solution-lam-'// trim(adjustl(strlambda)) &
            & //"-npts-" // trim(adjustl(strnpts)) //'.dat'
       open(11,file=filename)
       write(11,*) "x ", "u"

       ! Allocate matrices
       allocate(A(npts, npts))
       allocate(L, U, P, source = A)
       A = 0.0d0
       L = 0.0d0
       U = 0.0d0
       P = 0.0d0

       ! Allocate vectors
       allocate(x(npts))
       allocate(xtmp, b, y, source = x)
       x    = 0.0d0
       xtmp = 0.0d0 
       b    = 0.0d0

       ! Solve x
       call assemble_system(0.0d0, 1.0d0, npts, A, b, x, lambda)
       call dlufactor(A, L, U, flag)
       call fwdsub(L, b, y, flag)
       call backsub(U, y, x, flag)

       ! Write output
       do i = 1, npts
          write(11, *) dble(i)/dble(npts), x(i)
       end do

       close(11)

       deallocate(A,L,U,P,x,xtmp,b,y)

    end do

  end subroutine gesolve

  subroutine time_report(pivot)

    ! Problem setup
    integer :: npts
    logical :: pivot

    ! Matrices and vectors
    real(8), allocatable, dimension(:,:) :: A, L, U, P
    real(8), allocatable, dimension(:)   :: x, xtmp, b, y

    ! Timing variables
    real(8) :: total_start, total_finish, total_time
    real(8) :: assembly_start, assembly_finish, assembly_time
    real(8) :: factor_start, factor_finish, factor_time
    real(8) :: elim_start, elim_finish, elim_time

    integer, parameter :: sizes(6) = [200, 400, 800, 1000, 2000, 4000]
    real(8), parameter :: lambdas(2) = [0.0d0, 10.0d0]
    real(8) :: lambda
    integer :: flag, i, k

    ! Run for each lambda
    do k = 1, size(lambdas)

       lambda = lambdas(k)

       ! Write the time for each linear solve   
       if (pivot .eqv. .false.) then
          if (k .eq. 1) then 
             open(11, file='time_study_ge_lam0.dat')
          else
             open(11, file='time_study_ge_lam10.dat')
          end if
       else 
          if (k .eq. 1) then 
             open(11, file='time_study_gepp_lam0.dat')
          else
             open(11, file='time_study_gepp_lam10.dat')
          end if
       end if
       write(11, *) "size ", "factor_time ", "elim_time ", "total_time"

       ! Run for all problem sizes
       do i = 1, size(sizes)

          npts = sizes(i)

          ! Allocate matrices
          allocate(A(npts, npts))
          allocate(L, U, P, source = A)
          A = 0.0d0
          L = 0.0d0
          U = 0.0d0
          P = 0.0d0

          ! Allocate vectors
          allocate(x(npts))
          allocate(xtmp, b, y, source = x)
          x    = 0.0d0
          xtmp = 0.0d0 
          b    = 0.0d0

          call cpu_time(total_start)

          ! Assemble linear system
          call cpu_time(assembly_start)
          call assemble_system(0.0d0, 1.0d0, npts, A, b, x, lambda)
          call cpu_time(assembly_finish)
          assembly_time = assembly_finish-assembly_start
          print *, "assembled system in ", assembly_time

          call cpu_time(factor_start)
          if (pivot .eqv. .false.) then
             ! Ax = b  is now LUx = b
             call dlufactor(A, L, U, flag)
          else
             call dluppfactor(A, L, U, P, flag)
             b = matmul(P,b)
          end if
          call cpu_time(factor_finish)
          factor_time = factor_finish-factor_start
          print *, "factorized matrix in ", factor_time

          call cpu_time(elim_start)

          ! Solve Ly = b
          call fwdsub(L, b, y, flag)

          ! Solve Ux = y
          call backsub(U, y, x, flag)

          call cpu_time(elim_finish)
          elim_time = elim_finish-elim_start
          print *, "elimination completed in ", elim_time

          call cpu_time(total_finish)
          total_time = total_finish - total_start 
          print '("Time = ",f6.3," seconds.")', total_time

          ! Write output data
          write(11, *) npts, factor_time, elim_time, total_time

          deallocate(A,L,U,P,x,xtmp,b,y)

       end do

       close(11)

    end do

  end subroutine time_report

  ! Model problem to solve
  subroutine assemble_system(a, b, npts, V, rhs, u, lambda, homogenous)

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
    logical, intent(in), optional :: homogenous

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

    if(present(homogenous)) then
       ! Assemble the RHS
       rhs = 0.0d0
       rhs(1)    = 2.0d0
       rhs(npts) = 0.0d0
    else
       ! Assemble the RHS
       rhs = 0.0d0
       rhs(1)    = 2.0d0 + h*h*(3.0d0*(dble(i)*h)**2.0d0 - 0.5d0)
       rhs(npts) = 0.0d0 + h*h*(3.0d0*(dble(npts)*h)**2.0d0 - 0.5d0)
    end if


    ! Initial solution profile (straight line)
    do i = 1, npts
       u(i) = 2.0d0 + ((0.0d0+2.0d0)/(b-a))*(dble(i)*h)
    end do

  end subroutine assemble_system

end program test
