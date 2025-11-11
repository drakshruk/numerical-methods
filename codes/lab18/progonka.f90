subroutine progonka(studNum, Ag, Bg, alpha, beta, h, n, p, q, f, y)
    integer :: n, studNUm
    real(8) :: Ag, Bg, alpha(0:1), beta(0:1), h, y(0:n), x(0:n)
    real(8) :: mi, ni, Fi, c(0:n-1), d(0:n-1)
    real(8), external :: p, q, f

    x = [(1d0 + h*i, i = 0,n)]

  
    mi = (-4+2*h*h*q(x(1)))/(2+h*p(x(1)))
    ni = (2-h*p(x(1)))/(2+h*p(x(1)))
    Fi = 2*f(x(1), studNum)/(2+h*p(x(1)))

    c(0) = (mi+4)*alpha(1)/(2d0*h*alpha(0)+alpha(1)*(ni-3))
    d(0) = h*(2d0*Ag+alpha(1)*f(x(1), studNum)*h)/(mi+4)/alpha(1)
    do i = 1, n-1
        mi = (-4+2d0*h*h*q(x(i)))/(2+h*p(x(i)))
        ni = (2-h*p(x(i)))/(2+h*p(x(i)))
        Fi = 2d0*f(x(i), studNum)/(2+h*p(x(i)))
        c(i) = 1/(mi-ni*c(i-1))
        d(i) = h*h*Fi-ni*c(i-1)*d(i-1)
    end do

    y(n) = (2d0*Bg*h-(4+c(n-2))*c(n-1)*d(n-1)+c(n-2)*d(n-2))/(4d0*h-3-4d0*c(n-1)-c(n-2)*c(n-1))
    do i = n-1, 0, -1
        y(i) = c(i)*(d(i)-y(i+1))
    end do
end subroutine