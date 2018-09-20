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

end module nonlinear_transport
