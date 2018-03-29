module system

  implicit none
  
contains 

  !---------------------------------------------------------------------!
  ! Assemble Ax = b (Case (a))
  !---------------------------------------------------------------------!
  
  subroutine assemble_system1(npts, A, b, x0)

    implicit none

    integer, intent(in)  :: npts
    real(8), intent(out) :: A(npts,npts)
    real(8), intent(out) :: b(npts)
    real(8), intent(out) :: x0(npts)

    integer :: M, N
    integer :: i, j, k

    ! Zero all entries
    A = 0.0d0

    ! Size of the linear system
    M = npts ! nrows
    N = npts ! ncols

    ! Assemble matrix
    A(M,1) = 1.0d0
    A(1,N) = 1.0d0
    forall(i=1:M-1) A(i,i+1) = 1.0d0

    ! Assemble right hand side
    b = 0.0d0
    b(1) = 1.0d0

    ! Initial guess
    x0 = 0.0d0

  end subroutine assemble_system1

  !---------------------------------------------------------------------!
  ! Assemble Ax = b (Case (b))
  !---------------------------------------------------------------------!
  
  subroutine assemble_system2(npts, A, b, x0)

    implicit none

    integer, intent(in)  :: npts
    real(8), intent(out) :: A(npts,npts)
    real(8), intent(out) :: b(npts)
    real(8), intent(out) :: x0(npts)

    integer :: M, N
    integer :: i, j, k

    ! Zero all entries
    A = 0.0d0

    ! Size of the linear system
    M = npts ! nrows
    N = npts ! ncols

    ! Assemble matrix
    forall(i=1:M) A(i,i) = 1.0d0 ! diagonal
    j = 0 
    do i = 1, M-1, 2
       j = j + 1
       A(i,i+1) = dble(j) - 1.0d0
    end do

    ! Assemble right hand side
    b = 0.0d0
    b(1) = 1.0d0

    ! Initial guess
    x0 = 0.0d0

  end subroutine assemble_system2

  !-------------------------------------------------------------------!
  !Assemble Ax = b (GMRES)
  !-------------------------------------------------------------------!
  
  subroutine assemble_system3(npts, A, b, x0)

    implicit none

    integer, intent(in)  :: npts
    real(8), intent(out) :: A(npts,npts)
    real(8), intent(out) :: b(npts)
    real(8), intent(out) :: x0(npts)

    integer :: M, N
    integer :: i, j, k

    ! Zero all entries
    A = 0.0d0

    ! Size of the linear system
    M = npts ! nrows
    N = npts ! ncols

    ! Assemble matrix
    A(1,N) = 1.0d0
    forall(i=2:M) A(i,i-1) = 1.0d0

    ! Assemble right hand side
    b = 0.0d0
    b(1) = 1.0d0

    ! Initial guess
    x0 = 0.0d0

  end subroutine assemble_system3

end module system

!!$program test
!!$
!!$  use system 
!!$  use linear_algebra
!!$
!!$  implicit none
!!$
!!$  integer, parameter :: npts = 1000
!!$  real(8), parameter :: max_tol = 1.0d-8
!!$  integer, parameter :: max_it = 100000
!!$  real(8) :: x(npts), b(npts), A(npts,npts), P(npts, npts)
!!$  integer :: iter, flag, i, j
!!$  real(8) :: tol
!!$  complex(8) :: eigs(npts)
!!$
!!$  call assemble_system1(npts, A, b, x)
!!$  !do i = 1,npts
!!$  !   write(*,*) ( A(i,j), j=1,npts )
!!$  !enddo
!!$
!!$  eigs = sqrt(eigvals(matmul(A, transpose(A)))       )
!!$  !do i = 1, npts
!!$ !    !print *, realpart(eigs(i)), imagpart(eigs(i)), abs(eigs(i))
!!$ ! end do
!!$  print *, maxval(realpart(eigs)), minval(realpart(eigs))
!!$
!!$end program test
