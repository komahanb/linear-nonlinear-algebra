program test
  
  call check_classical
  
end program test

subroutine check_classical

  use linear_algebra

  implicit none  

  ! Test forward and backward substitutions
  test_sub : block
    
    real(8), allocatable :: L(:,:), U(:,:)
    real(8), allocatable :: b(:), x1(:), x2(:)
    integer :: i, j, m, n, flag

    ! Size of the linear system
    m = 100
    n = 100

    ! Memory allocations
    allocate(L(m,n))
    allocate(U(m,n))
    allocate(b(m))
    allocate(x1(m))
    allocate(x2(m))

    ! Fill in with values
    call random_number(L)
    call random_number(U)
    call random_number(b)

    ! Zero the upper/lower triangular part
    do j = 1, m 
       do i = 1, n
          if (i .gt. j) then
             U(i,j) = 0.0d0
          else if (i .eq. j) then           
          else
             L(i,j) = 0.0d0
          end if
       end do
    end do

    ! Check forward sub
    call fwdsub(L, b, x1, flag)
    print *, flag, x1 - solve(L,b)

    ! Check backward sub
    call backsub(U, b, x2, flag)
    print *, flag, x2 - solve(U, b)

    ! Clean up memory
    deallocate(L, U, b, x1, x2)

  end block test_sub
  
  test_iter_sol: block

    real(8), allocatable :: x(:), bb(:), A(:,:)
    integer :: iter, flag    
    real(8) :: tol,omega

    allocate(A(2,2), bb(2), x(2))

    A(1,1) = 12.0d0
    A(2,1) = 5.0d0

    A(1,2) = 3.0d0
    A(2,2) = 7.0d0

    bb(1) = 11.0d0
    bb(2) = 13.0d0

    x(1) = 0.0d0
    x(2) = 0.0d0

    call dseidel(A, bb, 1.0d0, 1000, 1.0d-8, x, iter, tol, flag)
    print *, x, iter, tol, flag, solve(A, bb)

    call dsor(A, bb, 1.1d0, 1000, 1.0d-8, x, iter, tol, flag)
    print *, x, iter, tol, flag, solve(A, bb)

    call djacobi(A, bb, 1000, 1.0d-8, x, iter, tol, flag)
    print *, x, iter, tol, flag, solve(A, bb)

    deallocate(A, bb, x)

  end block test_iter_sol

  print *, "problem5"
  problem5: block

    integer, parameter :: npts = 1000
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

    call assemble_system(0.0d0, 1.0d0, npts, A, b, x)
    call dsor(A, b, 1.99d0, max_it, max_tol, x, iter, tol, flag)
    print *, 'sor', 1.99d0, tol, iter 
!!$      
!!$    print *, 'SOR'   
!!$    do j = 1, 15
!!$       omega(j) = 1.84d0 + dble(j)/100.0d0
!!$       print *, omega(j)
!!$       call assemble_system(0.0d0, 1.0d0, npts, A, b, x)
!!$       call dsor(A, b, omega(j), max_it, max_tol, x, iter, tol, flag)
!!$       print *, 'sor', omega(j), tol, iter    
!!$    end do

    open(11, file='sor.dat')
    do i = 1, npts
       write(11, *) dble(i)/dble(npts), x(i), xtmp(i)
    end do
    close(11)    

    print *, 'gauss jacobi'
    call assemble_system(0.0d0, 1.0d0, npts, A, b, x)
    call djacobi(A, b, max_it, max_tol, x, iter, tol, flag)
    print *, 'jacobi', tol, iter    
    open(11, file='jacobi.dat')
    do i = 1, npts
       write(11, *) dble(i)/dble(npts), x(i), xtmp(i)
    end do
    close(11)

  end block problem5

end subroutine check_classical

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
