program main
    implicit none
    real(8), allocatable :: Y(:)
    real(8) :: intValue(1:3), a(1:3), b(1:3)
    real(8) :: simpson, h1, h2, h3, f1, f2, f3, eps, a1, h22, integral(1:3)
    integer n, i, j, studNum
    studNum = 1
    intValue(1:3) = [0.785398163351, 1.25495772545, 0.23702174104]
    eps = 1d-4
    a(1) = (eps/3*(0.5+0.02*studNum))**(1/(0.5+0.02*studNum))
    b(1) = (1-eps/3)**(1/(0.5+0.02*studNum))
    a(2) = 0d0
    b(2) = (3/eps)**1.5
    a(3) = sqrt(2*eps/3/(4-0.02*studNum))
    b(3) = log(3/eps/(3+0.02*studNum)**2)/(3+0.02*studNum)
    
    n = 13000
    allocate(Y(0:n))
    h1 = (b(1)-a(1))/n
    Y(0:n) = [(f1(a(1)+h1*i,studNum), i = 1,n)]
    integral(1) = simpson(Y,a(1),b(1),n)
    write(*,*) integral(1), intValue(1), abs(intValue(1)-integral(1)), abs(intValue(1)-integral(1))/intValue(1)
    deallocate(Y)

    n = 150000000
    allocate(Y(0:n))
    h2 = ((a(2)+b(2))/2-a(2))/n
    do i = 0, n
        Y(i) = f2(a(2)+h2*i,studNum)
    end do
    integral(2) = simpson(Y,a(2), (a(2)+b(2))/2,n)
    h2 = (b(2)-(a(2)+b(2))/2)/n
    do i = 0, n
        Y(i) = f2((a(2)+b(2))/2+h2*i,studNum)
    end do
    integral(2) = integral(2) + simpson(Y,(a(2)+b(2))/2, b(2),n)
    write(*,*) integral(2), intValue(2), abs(intValue(2)-integral(2)), abs(intValue(2)-integral(2))/intValue(2)
    deallocate(Y)

    n = 530
    allocate(Y(0:n))
    h3 = (b(3)-a(3))/n
    Y(0:n) = [(f3(a(3)+h3*i,studNum), i = 1,n)]
    integral(3) = simpson(Y,a(3),b(3),n)
    write(*,*) integral(3), intValue(3), abs(intValue(3)-integral(3)), abs(intValue(3)-integral(3))/intValue(3)
    pause
end program