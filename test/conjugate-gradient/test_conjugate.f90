program test  
  !call check_conjugate  
end program test

subroutine check_conjugate

  use linear_algebra

  implicit none

  print *, "conjugate gradients"
  problem4: block

    integer, parameter :: npts = 1000
    real(8), parameter :: max_tol = 1.0d-3
    integer, parameter :: max_it  = 100000

    real(8) :: x(npts), xtmp(npts), b(npts), A(npts,npts), P(npts, npts)
    integer :: iter, flag, i, j
    real(8) :: tol, omega(20)
    real(8) :: scale = 1.0d0

    !call assemble_system(0.0d0, 1.0d0, npts, A, b, x, P)
    !xtmp = solve(A,b)

    call assemble_system(0.0d0, 1.0d0, npts, A, b, x, P)

    !call dcg(A, b, max_it, max_tol, x, iter, tol, flag)    
    !call dsd(A, b, max_it, max_tol, x, iter, tol, flag)    
    call dpcg(A, P,  b, max_it, max_tol, x, iter, tol, flag)

    print *, 'sd', tol, iter 

    open(11, file='checkcg.dat')
    do i = 1, npts
       write(11, *) dble(i-1)/dble(npts), x(i), xtmp(i)
    end do
    close(11)

  end block problem4

end subroutine check_conjugate

subroutine assemble_system(a, b, npts, V, rhs, x, P)

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

end subroutine assemble_system
