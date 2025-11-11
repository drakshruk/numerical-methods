subroutine runge(u, v, n, h, x0, f1, f2)
    integer :: n
    real(8) :: u(0:n), v(0:n), h, x0, k1(1:4), k2(1:4)
    real(8), external :: f1, f2

    do i = 1,n
        k1(1) = h*f1(x0 + h*i, u(i-1), v(i-1))
        k2(1) = h*f2(x0 + h*i, u(i-1), v(i-1))

        k1(2) = h*f1(x0 + h*i+h*5d-1, u(i-1)+k1(1)*5d-1, v(i-1)+k2(1)*5d-1)
        k2(2) = h*f2(x0 + h*i+h*5d-1, u(i-1)+k1(1)*5d-1, v(i-1)+k2(1)*5d-1)

        k1(3) = h*f1(x0 + h*i+h*5d-1, u(i-1)+k1(2)*5d-1, v(i-1)+k2(2)*5d-1)
        k2(3) = h*f2(x0 + h*i+h*5d-1, u(i-1)+k1(2)*5d-1, v(i-1)+k2(2)*5d-1)

        k1(4) = h*f1(x0 + h*i+h, u(i-1)+k1(3), v(i-1)+k2(3))
        k2(4) = h*f2(x0 + h*i+h, u(i-1)+k1(3), v(i-1)+k2(3))
        
        u(i) = u(i-1) + 1d0/6*(k1(1) + 2*k1(2) + 2*k1(3) + k1(4))
        v(i) = v(i-1) + 1d0/6*(k2(1) + 2*k2(2) + 2*k2(3) + k2(4))
    end do
end subroutine