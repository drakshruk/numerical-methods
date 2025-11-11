subroutine runge(y0,x0,h,func,n,y1)
    implicit none
    integer :: i, n
    real(8) :: y0, x0, h
    real(8), external :: func
    real(8) :: y1(0:n), x1, k(1:4)

    y1(0) = y0
    x1 = x0

    do i = 1,n+1
        k(1) = h*func(x1, y1(i-1))
        k(2) = h*func(x1+h*5d-1, y1(i-1)+k(1)*5d-1)
        k(3) = h*func(x1+h*5d-1, y1(i-1)+k(2)*5d-1)
        k(4) = h*func(x1+h, y1(i-1)+k(3))
        x1 = x0 + h*i
        y1(i) = y1(i-1) + 1d0/6*(k(1) + 2*k(2) + 2*k(3) + k(4))
    end do
end subroutine