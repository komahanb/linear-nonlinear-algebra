program test

  ! Import dependencies
  use lapack, only : dgetrf
  !use linear_algebra 
  use direct_linear_solve, only : dlufactor, dluppfactor
  
  implicit none

!!$
!!$  integer, parameter :: sizes(6) = [200, 400, 800, 1000, 2000, 4000]
!!$  integer :: i, npts
!!$
!!$  ! Run for all problem sizes
!!$  do i = 1, size(sizes)
!!$
!!$  end do

contains

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
