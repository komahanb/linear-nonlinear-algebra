module nonlinear_transport

implicit none

contains
  
  pure real(8) function exact4(x)

    real(8), intent(in) :: x

    exact4 = -1.0d0 + sqrt(4.0d0-(8.0d0/3.0d0)*x-x*x*x/3.0d0)

  end function exact4

  pure real(8) function exact3(x)

    real(8), intent(in) :: x

    exact3 = -1.0d0 + sqrt(4.0d0-2.0d0*x-x*x)
    
  end function exact3
  
  pure real(8) function exact2(x)

    real(8), intent(in) :: x

    exact2 = -1.0d0 + sqrt(4.0d0-3.0d0*x)

  end function exact2
 
  ! Model problem to solve
  subroutine assemble_residual_jacobian(sparse, npts, phi, Q, V, R)

    logical, intent(in) :: sparse 
    integer, intent(in) :: npts 
    real(8), intent(in) :: phi(:) ! current state
    real(8), intent(in) :: Q(:) ! source

    real(8), intent(inout) :: V(:,:) ! jacbian 
    real(8), intent(inout) :: R(:) ! residual

    real(8), parameter :: a = 0.0d0, b = 1.0d0 ! bound of domain
    real(8), parameter :: L = 1.0d0

    real(8), parameter :: phi0 = 1.0d0
    real(8), parameter :: phiN = 0.0d0

    ! Local variables
    real(8) :: aa, bb, cc
    real(8) :: h    
    integer :: m, n, i, j

    m = npts
    n = npts

    ! Mesh spacing
    h = (b-a)/dble(npts+1)

    ! Assemble the Residual for supplied source
    R = 0.0d0

    ! First node
    aa = 1.0d0 + phi(1)
    bb = phi0 - 2.0d0*(1.0d0+phi(1)) - phi(2)
    cc = 1.0d0 - phi(1) + 2.0d0*phi(2)
    R(1) = 0.0d0 + &
         & - ( phi(2) - 2.0d0*phi(1) + phi0 &
         & + phi(2)*phi(2) &
         & - phi(2)*phi(1) &
         & - phi(1)*phi(1) &
         & + phi(1)*phi0 &
         & + h*h*Q(1)/0.1d0 )

    ! Interior nodes
    do i = 2, m-1
       aa = 1.0d0 + phi(i)
       bb = phi(i-1) - 2.0d0*(1.0d0+phi(i)) - phi(i+1)
       cc = 1.0d0 - phi(i) + 2.0d0*phi(i+1)
       R(i) = -(phi(i+1) - 2.0d0*phi(i) + phi(i-1) &
            & + phi(i+1)*phi(i+1) &
            & - phi(i+1)*phi(i) &
            & - phi(i)*phi(i) &
            & + phi(i)*phi(i-1) &
            & + h*h*Q(i)/0.1d0)
    end do

    ! Last node
    aa = 1.0d0 + phi(npts)
    bb = phi(npts-1) - 2.0d0*(1.0d0+phi(npts)) - phiN
    cc = 1.0d0 - phi(npts) + 2.0d0*phiN
    R(npts) = 0.0 + & 
         & - ( phiN - 2.0d0*phi(npts) + phi(npts-1) &
         & + phiN*phiN &
         & - phiN*phi(npts) &
         & - phi(npts)*phi(npts) &
         & + phi(npts)*phi(npts-1) &
         & + h*h*Q(npts)/0.1d0)

    ! Assemble the jacobian
    V = 0.0d0

    if (sparse .eqv. .false.) then
       
       do i = 1, m

          ! Create matrix coefficients
          if (i.eq.1) then
             ! First node
             aa = 1.0d0 + phi(1)
             bb = phi0 - 2.0d0*(1.0d0+phi(1)) - phi(2)
             cc = 1.0d0 - phi(1) + 2.0d0*phi(2)
          else if (i.eq.npts) then 
             ! Last node
             aa = 1.0d0 + phi(npts)
             bb = phi(npts-1) - 2.0d0*(1.0d0+phi(npts)) - phiN
             cc = 1.0d0 - phi(npts) + 2.0d0*phiN
          else
             ! Interior nodes
             aa = 1.0d0 + phi(i)
             bb = phi(i-1) - 2.0d0*(1.0d0+phi(i)) - phi(i+1)
             cc = 1.0d0 - phi(i) + 2.0d0*phi(i+1)
          end if

          do j = 1, n
             if (i .eq. j-1) then
                ! upper diagonal
                V(i,j) = cc
             else if (i .eq. j) then
                ! diagonal
                V(i,i) = bb
             else if (i .eq. j+1) then
                ! lower diagonal
                V(i,j) = aa
             else
             end if
          end do
       end do
       
    else

       do concurrent(i=1:n)

          if (i.eq.1) then

             ! First node
             aa = 1.0d0 + phi(1)
             bb = phi0 - 2.0d0*(1.0d0+phi(1)) - phi(2)
             cc = 1.0d0 - phi(1) + 2.0d0*phi(2)

             V(i,:) = [0.0d0, bb, cc]

          else if (i.eq.npts) then 

             ! Last node
             aa = 1.0d0 + phi(npts)
             bb = phi(npts-1) - 2.0d0*(1.0d0+phi(npts)) - phiN
             cc = 1.0d0 - phi(npts) + 2.0d0*phiN

             V(i,:) = [aa, bb, 0.0d0]

          else

             ! Interior nodes
             aa = 1.0d0 + phi(i)
             bb = phi(i-1) - 2.0d0*(1.0d0+phi(i)) - phi(i+1)
             cc = 1.0d0 - phi(i) + 2.0d0*phi(i+1)

             V(i,:) = [aa, bb, cc]

          end if

       end do

    end if

  end subroutine assemble_residual_jacobian

end module nonlinear_transport
