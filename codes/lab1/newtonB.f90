subroutine NPB(A, X, n, k, h)
    implicit none
    integer i, n, j, k
    real(8) A(0:n), X(0:n), h
    real(8), allocatable :: c(:)
    allocate(c(0:n))
    A(0:n) = 0d0
    A(n) = 1
    do i = 0, k - 1
        c(0:n) = A(0:n)
        A(0:n-1) = A(1:n)
        A(n) = 0d0
        do j = 0, n
            A(j) = (A(j) + c(j) * (-x(n) + h * i))/(h*(i+1))
        end do
    end do
    deallocate(c)
end subroutine

subroutine NewtonBackwards(P, Y, X, n)
    implicit none
    integer n, i, j
    real(8) :: X(0:n), Y(0:n), P(0:n), A(0:n), C(0:n,0:n), h
    call difMat(Y, C, n)
    write(*,*) '=====NewtonBackwards====='
    P(0:n) = 0d0
    P(n) = Y(n)
    h = X(1) - X(0)
    do i = 1, n
        call NPB(A, X, n, i, h)
        P(0:n) = P(0:n) + A(0:n) * C(i,n-i)
    end do
    do i = n, 0, -1
        write(*,'(f12.6)') P(i)
    end do
    write(*,*) " x_i        P_i"
    do j = 0, n
        write(*,*) x(j), sum([(P(i)*x(j)**(n-i), i = 0, n)])
    end do
end subroutine