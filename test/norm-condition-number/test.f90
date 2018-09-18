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

  integer, parameter :: npoints = 100000
  real(8) :: array(npoints)
  integer :: i

  ! Open file handle for output
  open(11, file='problem4.dat')  
  write(11,'(4A20)') 'm', 'L1Norm/L2Norm', 'L2Norm/LinfNorm', 'kappa'
  
  ! For increasing matrix size
  do matsize = 100, 1000, 100

     ! Memory allocations
     allocate(L(matsize,matsize))

     ! Conduct 100 random trials for statistics
     do trial = 1, ntrials

        call get_random_lower(L)     ! fix martrix entries

        ! Get the 2 norm condition number
        condn(trial) = cond(L, 2)
        condn(trial) = dble(trial) ! fix condition number evaluation

        l1norm(trial)   = norm(L(:,1), 1) ! fix matrix norms
        l2norm(trial)   = norm(L(:,1), 2)
        linfnorm(trial) = norm(L(:,1), 0)

     end do

     ratio(:,1) = l1norm/l2norm
     ratio(:,2) = l2norm/linfnorm

     ! Average the values
     avg(1) = sum(ratio(:,1))/dble(ntrials)
     avg(2) = sum(ratio(:,2))/dble(ntrials)
     avgcond = sum(condn)/dble(ntrials)
     write(11,*) matsize, avg(1), avg(2), avgcond
     print *,  matsize, avg(1), avg(2), avgcond

     deallocate(L)

  end do

  close(11)

  stop

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
  
  pure real(8) function norm(x, p)

    real(8), intent(in) :: x(:)
    integer, intent(in) :: p
    integer             :: length, k

    length = size(x)
    if ( p .eq. 0 ) then
       norm = maxval(abs(x))
    else
       do k = 1, length
          norm = norm + abs(x(k))**p
       end do
       norm = norm**(1.0d0/dble(p))
    end if

  end function norm

end program test
