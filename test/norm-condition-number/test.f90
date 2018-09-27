program test

  ! Import dependencies
  use lapack, only : dgetrf
  use linear_algebra , only : eig, eigvals, svdvals, svd

  implicit none

  integer, parameter   :: ntrials = 100  
  real(8), allocatable :: L(:,:)
  real(8)              :: l1norm(ntrials), l2norm(ntrials), linfnorm(ntrials)
  real(8)              :: ratio(ntrials,2)
  real(8)              :: condn(ntrials)
  real(8)              :: avg(2), avgcond
  integer              :: matsize
  integer              :: trial
  
  ! Open file handle for output
  open(11, file='problem4.dat')  
  write(11,'(4A20)') 'm', 'L1Norm/L2Norm', 'L2Norm/LinfNorm', 'kappa'
  
  ! For increasing matrix size
  do matsize = 1000, 10000, 1000

     ! Memory allocations
     allocate(L(matsize,matsize))

     ! Conduct 100 random trials for statistics
     do trial = 1, ntrials

        call get_random_lower(L)
                            
        ! Get the 2 norm condition number
        condn(trial) = cond(L,2)

        l1norm(trial)   = matnorm(L, 1)
        l2norm(trial)   = matnorm(L, 2)
        linfnorm(trial) = matnorm(L, 0)

     end do

     ratio(:,1) = l1norm/l2norm
     ratio(:,2) = l2norm/linfnorm

     ! Average the values
     avg(1) = sum(ratio(:,1))/dble(ntrials)
     avg(2) = sum(ratio(:,2))/dble(ntrials)
     avgcond = sum(condn)/dble(ntrials)

     ! Write output
     write(11,*) matsize, avg(1), avg(2), avgcond
     print *,  matsize, avg(1), avg(2), avgcond

     ! Free resources
     deallocate(L)

  end do

  close(11)

contains
  
  function random_normal(mean,stdev) result(c)

    real(8), intent(in) :: mean,stdev
    real(8)             :: c, r, theta
    real(8)             :: temp(2)
    real(8), parameter  :: pi = 4.d0*datan(1.d0)

    call random_number(temp)

    r = (-2.0d0*log(temp(1)))**0.5    
    theta = 2.0d0*pi*temp(2)

    c = mean+stdev*r*sin(theta)

  end function random_normal

  subroutine get_random_lower(L)

    real(8), intent(inout) :: L(:,:)
    integer :: i, j, m, n

    m = size(L,1)
    n = size(L,2)
  
    do i = 1, m 
       do j = 1, i
          L(i,j) = random_normal(0.0d0, sqrt(dble(m)))
       end do
    end do

  end subroutine get_random_lower  

  !===================================================================!
  ! Implement vector p-norm
  !===================================================================!
  
  pure real(8) function vecnorm(x, p)
    
    real(8), intent(in) :: x(:)
    integer, intent(in) :: p
    integer             :: length, k

    ! initialize value
    vecnorm = 0.0d0
    
    length = size(x)
    if ( p .eq. 0 ) then
       vecnorm = maxval(abs(x))
    else
       do k = 1, length
          vecnorm = vecnorm + abs(x(k))**dble(p)
       end do
       vecnorm = vecnorm**(1.0d0/dble(p))
    end if

  end function vecnorm

  !===================================================================!
  ! Implement Matrix p-norm
  !===================================================================!
  
  pure real(8) function matnorm(A, p)

    real(8), intent(in) :: A(:,:)
    integer, intent(in) :: p
    integer             :: length, k
    real(8) :: w(size(A,1))

    ! Create a unit vector
    w = 1.0d0
    w = w/vecnorm(w,p)

    ! Get the matrix vector product
    w = matmul(A,w)

    matnorm = vecnorm(w,p)
    
  end function matnorm

  !===================================================================!
  ! Implement Matrix p-cond
  !===================================================================!
  
  real(8) function cond(A, p)

    real(8), intent(in) :: A(:,:)
    integer, intent(in) :: p
    integer             :: length, k
    real(8) :: w(size(A,1))

    ! Get singular values
    w = svdvals(A)

    ! find the condition number
    cond = maxval(w)/minval(w)
    
  end function cond

end program test
