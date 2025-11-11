real(8) function trapez(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    trapez = 0d0
    do i = 1, n
        trapez = trapez + (y(i-1) + y(i)) * h / 2
    end do
end function

real(8) function lRect(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    lRect = 0d0

    do i = 0, n-1
        lRect = lRect + y(i) * h
    end do
end function

real(8) function rRect(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    rRect = 0d0

    do i = 1, n
      rRect = rRect + y(i) * h
    end do
end function

real(8) function cRect(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    cRect = 0d0

    do i = 1, n
      cRect = cRect + y(i) * h
    end do
end function

real(8) function simpson(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    if(mod(n,2) /= 0) n = n - 1
    
    h = (b - a) / n
    simpson = 0d0
    do i = 0, n/2-1
        simpson = simpson + y(2*i) + 4*y(2*i+1) + y(2*i+2)
    end do
    simpson = simpson * h / 3
end function

real(8) function newton(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    newton = 0d0
    do i = 1, n/3-1
        newton = newton + 2*y(3*i) + 3*(y(3*i-1) + y(3*i-2))
    end do
    newton = newton + y(0) + y(n) + 3*(y(n-1) + y(n-2))
    newton = newton * 3 * h / 8
end function

real(8) function f1(x, p)
    implicit none
    real(8) x, p
    f1 = (x**p/(1+x*x) + 2*x**(-p+2)/(-p+1)/(1+x*x)**2)
end function

real(8) function f2(x, p)
    implicit none
    real(8) x, p
    f2 = cos(x)**2*cos(p*x)
end function

real(8) function f3(x, p)
    implicit none
    real(8) x, p
    f3 = x**p*(log(1/x))**2
end function

real(8) function f4(x, p)
    implicit none
    real(8) x, p
    f4 = 1./((p-1)*x**2-2*x+1)**1.5
end function