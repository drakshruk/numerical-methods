real(8) function f(x, i)
    real(8) x
    integer i
    if(i == 1) then 
        f = 1
    else if(i == 2) then 
        f = x
    else if(i == 3) then 
        f = x*x
    else if(i == 4) then 
        f = 1/x
    else if(i == 5) then 
        f = log(x)
    else if(i == 6) then 
        f = exp(x)
    end if
end function

subroutine calc(a,b,f1,f2, y,n)
    implicit none
    real(8) a, b, f1(1:n), f2(1:n), y(1:n)
    real(8), allocatable :: coef(:)
    integer j,n
    allocate(coef(1:n))
    coef(:) = 0
    do j = 1, n
        coef(1) = coef(1) + f1(j)*f1(j)
        coef(2) = coef(2) + f1(j)*f2(j)
        coef(3) = coef(3) + f1(j)*f2(j)
        coef(4) = coef(4) + f2(j)*f2(j)
        coef(5) = coef(5) + f1(j)*y(j)
        coef(6) = coef(6) + f2(j)*y(j)
    end do
    a = (coef(5) * coef(4) - coef(6) * coef(3)) / (coef(1)*coef(4) - coef(2)*coef(3))
    b = (coef(1) * coef(6) - coef(2) * coef(5)) / (coef(1)*coef(4) - coef(2)*coef(3))
end subroutine

program mnk
    implicit none
    integer i, j, n
    real(8),allocatable :: X(:), Y(:), Func(:,:)
    real(8),allocatable :: a(:), b(:), s(:)
    real(8), external :: f
    real(8) smallS
    interface interface
        subroutine calc(a,b,f1,f2,y,n)
            implicit none
            real(8) a, b, f1(1:n), f2(1:n), y(1:n)
            integer n
        end subroutine
    end interface interface
    n = 6
    allocate(X(1:n),Y(1:n),Func(1:n,1:n),a(1:n),b(1:n),s(1:n))
    X = [0d0, 5d0, 10d0,15d0,20d0,25d0]
    Y = [21d0,39d0,51d0,63d0,70d0,90d0]
    do i = 1, n
        do j = 1, n
            Func(i,j) = f(X(j),i)
            if(abs(abs(Func(i,j)) - 10*9) > 10**9) then 
                Func(i,j) = 0d0
            end if
        end do
    end do

    call calc(a(1), b(1),Func(2:2,1:n), Func(1:1,1:n),Y,n) !ax + b
    call calc(a(2), b(2),Func(1:1,1:n), Func(5:5,1:n),Y,n) !a + blnx
    call calc(a(3), b(3),Func(1:1,1:n), Func(4:4,1:n),Y,n) !a + b/x
    call calc(a(4), b(4),Func(3:3,1:n), Func(1:1,1:n),Y,n) !axx + b
    call calc(a(5), b(5),Func(3:3,1:n), Func(2:2,1:n),Y,n) !axx + bx
    call calc(a(6), b(6),Func(1:1,1:n), Func(6:6,1:n),Y,n) !a + bex
    s(:) = 0d0
    do j = 1, n
        s(1) = s(1) + (a(1) * Func(2,j) + b(1) * Func(1,j) - y(j))**2
        s(2) = s(2) + (a(2) * Func(1,j) + b(2) * Func(5,j) - y(j))**2
        s(3) = s(3) + (a(3) * Func(1,j) + b(3) * Func(4,j) - y(j))**2
        s(4) = s(4) + (a(4) * Func(3,j) + b(4) * Func(1,j) - y(j))**2
        s(5) = s(5) + (a(5) * Func(3,j) + b(5) * Func(2,j) - y(j))**2
        s(6) = s(6) + (a(6) * Func(1,j) + b(6) * Func(6,j) - y(j))**2
    end do
    smallS = s(1)
    do i = 1,n
        if(s(i) < smallS) then
            smallS = s(i)
        end if
    end do
    write(*,*) 's =', s(1)
    write(*,*) 'f1 = ',a(1),'  x + ',b(1)
    write(*,*) 's =', s(2)
    write(*,*) 'f2 = ',a(2),'    + ',b(2),'lnx'
    write(*,*) 's =', s(3)
    write(*,*) 'f3 = ',a(3),'    + ',b(3),'/x'
    write(*,*) 's =', s(4)
    write(*,*) 'f4 = ',a(4),'x^2 + ',b(4)
    write(*,*) 's =', s(5)
    write(*,*) 'f5 = ',a(5),'x^2 + ',b(5),'x'
    write(*,*) 's =', s(6)
    write(*,*) 'f6 = ',a(6),'    + ',b(6),'e^x'
    write (*,*) 'smallest S =', smallS

    do i = 1, n
        write(*,*)  x(i), a(1)*x(i) + b(1)
    end do
    write(*,*)
    do i = 1, n
        write(*,*)  x(i), y(i)
    end do
    pause
end program