program main
    integer :: studNum, n
    real(8) :: y0, x0, a, b
    real(8) :: h, eps, c
    real(8), allocatable :: y1(:)

    studNum = 1
    n = 100
    a = 1d0
    b = 2d0
    c = 1 + 1d-3*studNum
    h = (b-a)/n

    x0 = 1d0
    y0 = 1d0 + c
    eps = 1d-4

    allocate(y1(0:n))
    call runge(y0,x0,h,f,n,y1)
    do i = 0, n
        write(*,*) x0+h*i, y1(i), abs(y1(i)-corF(x0+h*i)),abs(y1(i)-corF(x0+h*i))/corF(x0+h*i)
    end do
    write(*,*)
    call adams(y1(0:2),x0,h,n,eps,f,y1)
    do i = 0, n
        write(*,*) x0+h*i, y1(i), abs(y1(i)-corF(x0+h*i)),abs(y1(i)-corF(x0+h*i))/corF(x0+h*i)
    end do
    write(*,*)
    call runge(y0,x0,h,f,n,y1)
    call milln(y1(0:2),x0,h,n,eps,f,y1)
    do i = 0, n
        write(*,*) x0+h*i, y1(i), abs(y1(i)-corF(x0+h*i)), abs(y1(i)-corF(x0+h*i))/corF(x0+h*i)
    end do
    
    do i = 0, n
        write(*,*) corF(x0+h*i)
    end do
    pause
contains
    real(8) function f(x,y)
        real(8) :: x, y
        f = -x/4d0/c - 5*c/4d0/x/x + y*y/4d0/c
    end function
    
    real(8) function corF(x)
        real(8) :: x
        corF = (c + x*sqrt(x))/x
    end function
end program