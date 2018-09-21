module linear_transport

implicit none

contains

  pure real(8) function exact1(x)

    real(8), intent(in) :: x
    real(8), parameter  :: L     = 1.0d0
    real(8), parameter  :: gamma = 1.d0
    real(8), parameter  :: U     = 1.0d0

    exact1 = 1.0d0 - 1.0d0*(exp(x*U/gamma)-1.0d0)/(exp(L*U/gamma)-1.0d0)
    
  end function exact1
  
  ! Model problem to solve
  subroutine assemble_system(a, b, npts, V, rhs, sparse)

    logical, intent(in)  :: sparse 
    real(8), intent(in)  :: a, b ! bound of domain
    real(8), intent(out) :: V(:,:)
    real(8), intent(out) :: rhs(:)
    integer              :: npts ! number of points

    ! Local variables
    real(8) :: aa, bb, cc
    real(8) :: h    
    integer :: m, n, i, j

    real(8), parameter  :: L     = 1.0d0
    real(8), parameter  :: gamma = 1.0d0
    real(8), parameter  :: U     = 1.0d0

    ! Mesh spacing
    h = (b-a)/dble(npts+1)

    ! Problem parameters
    aa = -u/(2.0d0*h) - gamma/(h*h)
    bb = 2.0d0*gamma/(h*h)
    cc = u/(2.0d0*h) - gamma/(h*h)
    
    ! Prepare matrix
    V = 0.0d0
    m = npts
    n = npts
    if (sparse .eqv. .true.) then           

       V(1,:) = [0.0d0, bb, cc]
       do concurrent(i=2:n-1)
          V(i,:) = [aa, bb, cc]
       end do
       V(n,:) = [aa, bb, 0.0d0]
       
    else

       do i = 1, m
          do j = 1, n
             if (i .eq. j-1) then
                V(i,j) = aa
             else if (i .eq. j) then           
                V(i,i) = bb                
             else if (i .eq. j+1) then           
                V(i,j) = cc
             else
                ! skip zeros
             end if
          end do
       end do

    end if

    ! Assemble the RHS
    rhs = 0.0d0
    rhs(1) = 1.0d0
    rhs(npts) = 0.0d0

  end subroutine assemble_system

end module linear_transport
