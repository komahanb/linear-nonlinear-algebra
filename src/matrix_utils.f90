module matrix_utils

contains
 
  !===================================================================!
  ! Setup Toepliz matrix with supplied columns
  !===================================================================!

  subroutine toepliz(col, A)

    real(8), allocatable :: A(:,:)
    real(8) :: col(:)
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

end module matrix_utils
