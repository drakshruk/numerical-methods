subroutine L(A, X, n, k)
    implicit none
    integer i, k, n, j
    real(8) A(0:n), X(0:n), tmp
    real(8), allocatable :: c(:)
    allocate(c(0:n))
    A(0:n) = 0d0
    A(n) = 1d0
    tmp = 1d0
    do i = 0, n
        if(i /= k) then
            c(0:n) = A(0:n)
            A(0:n-1) = A(1:n)
            A(n) = 0d0
            do j = 0, n
                A(j) = A(j) + c(j) * (-1.0) * x(i) 
            end do
            tmp = tmp * (x(k) - x(i))
        end if
    end do
    A(0:n) = A(0:n) / tmp
    deallocate(c)
end subroutine

subroutine lagrange(P, Y, X, n)
    implicit none
    integer n, i,j
    real(8) :: X(0:n), Y(0:n), P(0:n), A(0:n)
    P(0:n) = 0d0
    do i = 0, n
        call L(A, X, n, i)
        P(0:n) = P(0:n) + A(0:n) * Y(i)
    end do
    write(*,*) '=====Lagrange====='
    do i = n, 0, -1
        write(*,'(f12.6)') P(i)
    end do
    write(*,*) " x_i        P_i"
    do j = 0, n
        write(*,*) x(j), sum([(P(i)*x(j)**(n-i), i = 0, n)])
    end do
end subroutine