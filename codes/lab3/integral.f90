real(8) function integral(a, b, N, f1, f2, weight)
    implicit none
    real(8) a, b, dx, f1(0:), f2(0:), weight(0:)
    integer N, i
    dx = (b - a) / N
    integral = 0.
    do i = 1, N
      integral = integral + (f1(i-1)*f2(i-1)*weight(i-1) + f1(i)*f2(i)*weight(i)) * dx / 2
    end do
end function integral

real(8) function polyValue(x, n, A)
    real(8) :: x, A(0:n)
    integer n
    polyValue = 0.
    do i = 0, n
        polyValue = polyValue + A(i)*x**i
    end do
end function polyValue

real(8) function fPoly(x,a)
    real(8) :: x, a
    fPoly = abs((x+1.5)**2 - 2*a*(x+1.5) + 1.5)
end function fPoly