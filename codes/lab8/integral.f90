real(8) function simpson(Y,a,b,n)
    implicit none
    integer n, i
    real(8) a, b, h, Y(0:n)

    h = (b - a) / n
    simpson = 0d0
    do i = 0, n/2-1
        simpson = simpson + y(2*i) + 4*y(2*i+1) + y(2*i+2)
    end do
    simpson = simpson * h / 3
end function

real(8) function f1(x,n)
    implicit none
    real(8) x
    integer n
    f1 = x**(-0.5 + 0.02*n)*sin((0.5 + 0.02*n)*log(x))/log(x)
end function

real(8) function f2(x,n)
    implicit none
    real(8) x
    integer n
    f2 = x**(1/3.)/(1+x*sqrt(3.)+0.02*n+x*x)
end function

real(8) function f3(x,n)
    implicit none
    real(8) x
    integer n
    f3 = log(1/x)*exp(-(3+0.02*n)*x)*sin((4+0.02*n)*x)
end function