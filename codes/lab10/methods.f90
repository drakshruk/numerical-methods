subroutine dichotomy(a,b,N,f,eps, res)
    implicit none
    real(8) a, b, a1, b1, eps, res
    integer N
    real(8), external :: f
    
        ! plot exp(2*x) + 3*x - 1 w l lw 2 lc 'blue' t '1', 0 w l lc 'black' t '', 'dichotomy_1.txt' w l lc 'red' t '2'
        ! plot x**4 - (1+a)*x**3 - (4-a)*x*x + (4+4*a)*x - 4*a - sqrt(x-1) w l lc 'blue' t '1', 0 w l lc 'black' t '', 'dichotomy_2.txt' w p lc 'red' t '2'
    if(abs(f((b+a)/2,N)) < eps) then
        res = (a+b)/2
        return
    else if(abs(f(a,N)) < eps) then
        res = a
        return
    else if(abs(f(b,N)) < eps) then
        res = b
        return
    else if(f((b+a)/2,N)*f(a,N) < 0d0) then
        a1 = a
        b1 = (b+a)/2
    else
        b1 = b
        a1 = (b+a)/2
    end if
    write(*,*) (b+a)/2, 0
    write(*,*) (b1+a1)/2, 0

    do while(abs(b1-a1) > eps)
        if(abs(f((b1+a1)/2,N)) <= eps) then
            res = (a1+b1)/2
            return
        else if(f((b1+a1)/2,N)*f(a1,N) < 0d0) then
            a1 = a1
            b1 = (b1+a1)/2
        else if(f((b1+a1)/2,N)*f(b1,N) < 0d0) then
            b1 = b1
            a1 = (b1+a1)/2
        end if
        write(*,*) (b1+a1)/2, 0
    end do
    res = (a1+b1)/2
    write(*,*)
end subroutine dichotomy

subroutine newton(a,b,N,f,fder,eps, res)
    implicit none
    real(8) a, b, eps, x, x1, res
    real(8) x0, y0
    integer N
    real(8), external :: f, fder
    if(abs(f(a,N)) < eps) then
        res = a
        return
    else if(abs(f(b,N)) < eps) then
        res = b
        return
    end if
    x = b
    x1 = x - f(x,N)/fder(x,N)
    do while(abs(x1-x) > eps)
        x1 = x
        x = x1 - f(x1,N)/fder(x1,N)
        ! plot exp(2*x) + 3*x - 1 w l lw 2 lc 'blue' t '1', 0 w l lc 'black' t '', 'newton_1.txt' w l lc 'red' t '2'
        ! plot x**4 - (1+a)*x**3 - (4-a)*x*x + (4+4*a)*x - 4*a - sqrt(x-1) w l lw 2 lc 'blue' t '1', 0 w l lc 'black' t '', 'newton_2.txt' w l lc 'red' t '2'
        ! вывод для графика
            x0 = (x1-x)/sqrt((x1-x)**2+(f(x1,n)-f(x,n))**2)
            y0 = (f(x1,n)-f(x,n))/sqrt((x1-x)**2+(f(x1,n)-f(x,n))**2)

            write(*,*) x - x0*10, f(x,N) - y0*10
            write(*,*) x - x0*5, f(x,N) - y0*5
            write(*,*) x - x0*3, f(x,N) - y0*3
            write(*,*) x - x0, f(x,N) - y0
            write(*,*) x, f(x,N)
            write(*,*) x1 + x0, f(x1,N) + y0
            write(*,*) x1 + x0*3, f(x1,N) + y0*3
            write(*,*) x1 + x0*5, f(x1,N) + y0*5
            write(*,*) x1 + x0*10, f(x1,N) + y0*10
            write(*,*)
            write(*,*)
        !
    end do
    res = x
end subroutine newton

subroutine secant(a,b,N,f,eps, res)
    implicit none
    real(8) a, b, eps, x, x1, x2, res
    integer N
    real(8), external :: f
    if(abs(f(a,N)) < eps) then
        res = a
        return
    else if(abs((f(b,N))) < eps) then
        res = b
        return
    end if
    x = a
    x1 = b
    x2 = x1 - f(x1,N)/(f(x1,N)-f(x,N))*(x1-x)
    
    write(*,*) x, f(x,N)
    write(*,*) x1, f(x1,N)
    write(*,*) x2, f(x2,N)
    do while(abs(x2-x1) > eps)
        x = x2
        x2 = x2 - f(x2,N)/(f(x2,N)-f(x1,N))*(x2-x1)
        x1 = x
        ! вывод для графика
            write(*,*) x, f(x,N)
        !
    end do
    res = x2
end subroutine secant

subroutine iterations(a,b,N,f,eps, res)
    real(8) a, b, eps, x1, x2, res
    integer N, ni
    real(8), external :: f
    ni = (b-a)/0.01
    if(abs(f(a,N) - a) < eps) then
        res = a
        return
    else if(abs(f(b,N) - b) < eps) then
        res = b
        return
    end if
    
    ! plot 1/3*(1-exp(2*x)) w l lc 'blue' title '1', x w l lc 'red' title '2', 'simpiter_1.txt' w p lc 'black' title '3'
    ! plot ((1+a)*x**3 + (4-a)*x*x - (4+4*a)*x + 4*a + sqrt(x-1))**(1d0/4) w l lc 'blue' title '1', x w l lc 'red' title '2', 'simpiter_1.txt' w l lc 'black' title '3'
    x1 = f(a,N)
    x2 = f(x1,N)
    write(*,*) a, f(a,N)
    write(*,*) x1, x1
    write(*,*) x1, f(x1,N)
    write(*,*) x2, x2
    ! write(*,*) x2, f(x1,N)
    ! write(*,*) x2, x2
    do while(abs(x2-x1) > eps)
        x1 = x2
        x2 = f(x2,N)
        write(*,*) x1, f(x1,N)
        write(*,*) x2, x2
    end do
    res = x2
    write(*,*)
    write(*,*)
end subroutine iterations