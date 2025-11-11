subroutine milln(y0, x0, h, n, eps, func, y1)
    implicit none
    integer :: i, n
    real(8) :: y0(0:2), x0, h, eps
    real(8), external :: func
    real(8) :: y1(0:n), y
    
    y1(0:2) = y0(0:2)
    do i = 3, n
        y = 0d0
        y1(i) = y1(i-2) + h/3d0*(7d0*func(x0+(i-1)*h, y1(i-1)) - 2d0*func(x0+(i-2)*h, y1(i-2)) + func(x0+(i-3)*h,y1(i-3)))
        ! do while(abs(y1(i)-y) > eps)
        !     y = y1(i)
            y1(i) = y1(i-2) + h/3d0*(func(x0+(i-2)*h, y1(i-2)) + 4d0*func(x0+(i-1)*h, y1(i-1)) + func(x0+i*h,y1(i)))
        ! end do
    end do
end subroutine