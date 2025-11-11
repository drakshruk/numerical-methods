program main
    integer :: studNum, n
    real(8) :: u0, v0, x0, a, b, h
    real(8), allocatable :: u1(:), v1(:)

    studNum = 14
    a = 1d0
    b = 1.1d0
    h = 1d-4
    n = 1000

    x0 = 1d0
    u0 = 1d0
    v0 = 6d0

    allocate(u1(0:n), v1(0:n))
    u1(0) = u0
    v1(0) = v0
    
    call runge(u1, v1, n, h, x0, fu, fv)

    do i = 0, n
        write(*,*) x0+h*i, realU(x0+h*i), u1(i), abs(u1(i)-realU(x0+h*i)), abs(realU(x0+h*i)-u1(i))/realU(x0+h*i)
    end do
    write(*,*)
    do i = 0, n
        write(*,*) x0+h*i, realV(x0+h*i), v1(i), abs(realV(x0+h*i)-v1(i)), abs(realV(x0+h*i)-v1(i))/realV(x0+h*i)
    end do
    write(*,*)
    
    call scheme(u1, v1, h, n, studNum)
    do i = 0, n
        write(*,*) x0+h*i, realU(x0+h*i), u1(i), abs(u1(i)-realU(x0+h*i)), abs(realU(x0+h*i)-u1(i))/realU(x0+h*i)
    end do
    write(*,*)
    do i = 0, n
        write(*,*) x0+h*i, realV(x0+h*i), v1(i), abs(realV(x0+h*i)-v1(i)), abs(realV(x0+h*i)-v1(i))/realV(x0+h*i)
    end do
    
    pause
    contains
    real(8) function fu(x,u,v)
        real(8) :: x, u, v
        fu = v
    end function

    real(8) function fv(x,u,v)
        real(8) :: x, u, v
        fv = 4d0*u
    end function

    real(8) function realU(x)
        real(8) :: x
        realU = 2d0*exp(2d0*x-2) - exp(2-2d0*x)
    end function

    real(8) function realV(x)
        real(8) :: x
        realV = 4d0*exp(2d0*x-2) + 2d0*exp(2-2d0*x)
    end function
end program