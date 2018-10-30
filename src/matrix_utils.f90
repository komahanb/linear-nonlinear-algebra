module matrix_utils

  use iso_fortran_env, only : dp => REAL64

contains
 
  !===================================================================!
  ! Setup Toepliz matrix with supplied columns
  !===================================================================!

  subroutine toepliz(col, A)

    real(dp), allocatable :: A(:,:)
    real(dp) :: col(:)
    integer :: m, i, j

    m = size(A, 1)
    do j = 1, m
       do i = 1, m
          if (i.eq.j) then
             A(i,j) = col(1)
          else if (i+1.eq.j .or.i-1.eq.j) then
             A(i,j) = col(2)
             A(j,i) = col(2)
          else if (i+2.eq.j .or.i-2.eq.j) then
             A(i,j) = col(3)
             A(j,i) = col(3)
          else if (i+3.eq.j .or.i-3.eq.j) then
             A(i,j) = col(4)
             A(j,i) = col(4)
          else if (i+4.eq.j .or.i-4.eq.j) then
             A(i,j) = col(5)
             A(j,i) = col(5)
          end if
       end do
    end do

  end subroutine toepliz
  
  subroutine print(matrix)

    real(dp), intent(in) :: matrix(:,:)
    integer :: i, j, m, n

    m = size(matrix, 1)
    n = size(matrix, 2)

    do i = 1, m
       write(*,'(100g15.5)') ( matrix(i,j), j=1,n )
    end do
  end subroutine print

end module matrix_utils
