real(8) function y(x)
    implicit none
    real(8) x
    y = exp(sin(100*x + 1))
end function
! 100*cos(100*x+1)*exp(sin(100*x + 1))
real(8) function yRealDer(x)
    implicit none
    real(8) x, y
    yRealDer = 100*cos(100*x+1)*y(x)
end function

real(8) function yder(x,h)
    implicit none
    real(8) x, h, y
    yder = 1./2./h*(y(x) - y(x-2*h))
end function

real(8) function runge(x,h,r,p)
    implicit none
    real(8) :: x, h, r, yder
    integer p
    runge = yder(x,h) + (yder(x,h) - yder(x,r*h))/(r**p-1)
end function