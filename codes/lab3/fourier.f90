subroutine fourier(m)
    implicit none
    real(8) :: Ai(0:m), Bi(0:m), pi, x0, h, a, y0, x00, l
    real(8), allocatable :: func(:), weight(:), func1(:), func2(:)
    integer n, m, ni, j, i, studNum
    interface integral_interface
        real(8) function integral(a, b, N, f1, f2, weight)
            implicit none
            real(8) a, b, f1(0:), f2(0:), weight(0:)
            integer N
        end function integral
    end interface integral_interface

    n = 40
    ni = 40
    studNum = 1
    allocate(weight(0:ni), func(0:ni), func1(0:ni), func2(0:ni))
    weight(0:ni) = 1d0
    a = 2d-1 + 0.002*studNum
    pi = 4 * atan(1d0)
    l = 2d0
    x0 = -l
    h = 2*l/ni
    do j = 0, ni
        func(j) = fFourier(x0+j*h, a)
    end do
    do i = 0, m
        do j = 0, ni
            func1(j) = cos(i*(x0 + j*h)*pi/l)
            func2(j) = sin(i*(x0 + j*h)*pi/l)
        end do
        Ai(i) = 1 / l * integral(-l, l, ni, func, func1, weight)
        Bi(i) = 1 / l * integral(-l, l, ni, func, func2, weight)
    end do
    write(*,*) "fourier", m
    write(*,*) " x         f       fourier        |fourier-f|        |fourier-f|/|f|"
    x00 = -21d-1
    do j = 0, 40
        x00 = x00 + 0.1
        y0 = 0d0
        do i = 0, m
            y0 = y0 + Ai(i) * cos(x00*i) + Bi(i) * sin(x00*i)
        end do
        write(*,'(5(f12.8))') x00, y0, fFourier(x00, a), abs(fFourier(x00, a)-y0), abs(fFourier(x00, a) - y0)/fFourier(x00, a)
    end do
    contains
    
    real(8) function fFourier(x,a)
        implicit none
        real(8) :: x, a
        x = x + 2
        x = x - int(x/4) - 2
        if(x >= -2 .and. x < -a) then
            fFourier = 0
        else if(x >= -a .and. x < 0) then
            fFourier = x + a
        else if(x >= 0 .and. x < a) then
            fFourier = -x + a
        else if(x >= a .and. x < 2) then
            fFourier = 0
        end if
    end function fFourier
end subroutine