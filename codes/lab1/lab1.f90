program main
    implicit none    
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=30) :: output_string
        end function to_lx
    end interface func
    integer n, i, j
    real(8) h, x0, tmp
    real(8), allocatable :: X(:), Y(:), P(:)
    n = 70
    allocate(X(0:n), Y(0:n), P(0:n))
    h = 2./n
    x0 = 1
    X = [(-1 + h*i, i = 0,n)]
    Y = [(exp(X(i)), i = 0,n)]
    call lagrange(P, Y, X, n)
    call NewtonForwards(P, Y, X, n)
    call NewtonBackwards(P, Y, X, n)
    call Newton(P, Y, X, n)
    call gorner(x0, P, n)

    deallocate(X, Y)
    n = 15
    allocate(X(0:n), Y(0:n))
    X = [(-1d0 + 2d0/15*i, i = 0, n)]
    Y = [(exp(X(i)), i = 0, n)]
    write(*,*) " x_i  y_i  p_i  abs_pogr   otn_pogr"

    do j = 0, n
        tmp = sum([(P(i)*x(j)**(7-i), i = 0, 7)])
        write(*,*) x(j), Y(j), tmp, (abs(tmp - Y(j))), (abs(tmp - Y(j))/abs(y(j)))
    end do
    pause
end program