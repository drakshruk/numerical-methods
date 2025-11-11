program main
    integer :: n
    real(8)  :: intVal(1:3), integral(1:3), a(1:3), b(1:3), c(1:3), d(1:3)
    real(8)  :: f1, f2, f3, trapez, cRect, simpson
    real(8), allocatable ::  X(:), Y(:), func(:,:)
    interface interface
        subroutine mk(n)
            integer n
        end subroutine
    end interface
    intVal = [ 4.66666666667,0.0408219945203,0.107250184799]
    a = [0,3,0]
    b = [1,4,1]
    c = [0,1,1]
    d = [2,2,2]

    n = 1000
    allocate(X(0:n),Y(0:n),func(0:n,0:n))
    X = [(a(1)+(b(1)-a(1))/n*i, i = 0,n)]
    Y = [(c(1)+(d(1)-c(1))/n*i, i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f1(x(j),y(i))
        enddo
    enddo
    integral(1) = trapez((b(1)-a(1))/n,(d(1)-c(1))/n,n,n,func)
    integral(2) = simpson((b(1)-a(1))/n,(d(1)-c(1))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(1), integral(1), abs(intVal(1)-integral(1)), abs(intVal(1)-integral(1))/intVal(1) 
    write(*,*)
    write(*,*) intVal(1), integral(2), abs(intVal(1)-integral(2)), abs(intVal(1)-integral(2))/intVal(1)

    X = [(a(1)+(b(1)-a(1))/n*(i+0.5), i = 0,n)]
    Y = [(c(1)+(d(1)-c(1))/n*(i+0.5), i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f1(x(i),y(j))
        enddo
    enddo
    integral(3) = cRect((b(1)-a(1))/n,(d(1)-c(1))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(1), integral(3), abs(intVal(1)-integral(3)), abs(intVal(1)-integral(3))/intVal(1) 
    deallocate(x,y,func)

    n = 1000
    allocate(X(0:n),Y(0:n),func(0:n,0:n))
    X = [(a(2)+(b(2)-a(2))/n*i, i = 0,n)]
    Y = [(c(2)+(d(2)-c(2))/n*i, i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f2(x(j),y(i))
        enddo
    enddo
    integral(1) = trapez((b(2)-a(2))/n,(d(2)-c(2))/n,n,n,func)
    integral(2) = simpson((b(2)-a(2))/n,(d(2)-c(2))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(2), integral(1), abs(intVal(2)-integral(1)), abs(intVal(2)-integral(1))/intVal(1) 
    write(*,*)
    write(*,*) intVal(2), integral(2), abs(intVal(2)-integral(2)), abs(intVal(2)-integral(2))/intVal(1) 

    X = [(a(2)+(b(2)-a(2))/n*0.5+(b(2)-a(2))/n*i, i = 0,n)]
    Y = [(c(2)+(d(2)-c(2))/n*0.5+(d(2)-c(2))/n*i, i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f2(x(i),y(j))
        enddo
    enddo
    integral(3) = cRect((b(2)-a(2))/n,(d(2)-c(2))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(2), integral(3), abs(intVal(2)-integral(3)), abs(intVal(2)-integral(3))/intVal(1) 
    deallocate(x,y,func)

    n = 1000
    allocate(X(0:n),Y(0:n),func(0:n,0:n))
    X = [(a(3)+(b(3)-a(3))/n*i, i = 0,n)]
    Y = [(c(3)+(d(3)-c(3))/n*i, i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f3(x(j),y(i))
        enddo
    enddo
    integral(1) = trapez((b(3)-a(3))/n,(d(3)-c(3))/n,n,n,func)
    integral(2) = simpson((b(3)-a(3))/n,(d(3)-c(3))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(3), integral(1), abs(intVal(3)-integral(1)), abs(intVal(3)-integral(1))/intVal(1) 
    write(*,*)
    write(*,*) intVal(3), integral(2), abs(intVal(3)-integral(2)), abs(intVal(3)-integral(2))/intVal(1) 

    X = [(a(3)+(b(3)-a(3))/n*0.5+(b(3)-a(3))/n*i, i = 0,n)]
    Y = [(c(3)+(d(3)-c(3))/n*0.5+(d(3)-c(3))/n*i, i = 0,n)]
    do i = 0, n
        do j = 0,n
            func(i,j) = f3(x(i),y(j))
        enddo
    enddo
    integral(3) = cRect((b(3)-a(3))/n,(d(3)-c(3))/n,n,n,func)
    write(*,*)
    write(*,*) intVal(3), integral(3), abs(intVal(3)-integral(3)), abs(intVal(3)-integral(3))/intVal(1) 
    deallocate(x,y,func)
    write(*,*)
    call mk(10000)
    write(*,*)
    call mk(100000)
    write(*,*)
    call mk(1000000)
    write(*,*)
    pause
end program