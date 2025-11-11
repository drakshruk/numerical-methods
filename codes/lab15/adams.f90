subroutine adams(y0, x0, h, n, eps, func, y1)
    implicit none
    integer :: i, n
    real(8) :: y0(0:2), x0, h, eps
    real(8), external :: func
    real(8) :: y1(0:n), y
    
    y1(0:2) = y0(0:2)
    do i = 3, n
        y = 0d0
        y1(i) = y1(i-1) + h*(23d0/12*func(x0+(i-1)*h, y1(i-1)) - 4d0/3*func(x0+(i-2)*h, y1(i-2)) + 5d0/12*func(x0+(i-3)*h,y1(i-3)))
        ! do while(abs(y1(i)-y) > eps)
        !     y = y1(i)
            y1(i) = y1(i-1) + h*(5d0/12*func(x0+i*h, y1(i)) + 2d0/3*func(x0+(i-1)*h, y1(i-1)) - 1d0/12*func(x0+(i-2)*h,y1(i-2)))
        ! end do
    end do
end subroutine