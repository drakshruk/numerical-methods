subroutine difMat(Y, A, n)
    implicit none
    integer n, i
    real(8) :: Y(0:n), A(0:n,0:n)
    A(0:n,0:n) = 0d0
    do i = 0, n
        if(i == 0) then
            A(i,0:n) = y(0:n)
        else if(i == 1) then
            A(i,0:n-1) = y(1:n) - y(0:n-1)
            A(i,n) = 0d0
        else
            A(i,0:n-1) = A(i-1, 1:n) - A(i-1,0:n-1)
            A(i,n) = 0d0
        end if
    end do
end subroutine

subroutine NPF(A, X, n, k, h)
    implicit none
    integer i, n, j, k
    real(8) A(0:n), X(0:n), h
    real(8), allocatable :: c(:)
    allocate(c(0:n))
    A(0:n) = 0d0
    A(n) = 1d0
    do i = 0, k - 1
        c(0:n) = A(0:n)
        A(0:n-1) = A(1:n)
        A(n) = 0d0
        do j = 0, n
            A(j) = (A(j) + c(j) * (-x(0) - h * i))/(h*(i+1))
        end do
    end do
    deallocate(c)
end subroutine

subroutine NewtonForwards(P, Y, X, n)
    implicit none
    integer n, i,j
    real(8) :: X(0:n), Y(0:n), P(0:n), A(0:n), C(0:n,0:n), h
    call difMat(Y, C, n)
    P(0:n) = 0d0
    P(n) = Y(0)
    h = X(1) - X(0)
    do i = 1, n
        call NPF(A, X, n, i, h)
        P(0:n) = P(0:n) + A(0:n) * C(i,0)
    end do
    write(*,*) '=====NewtonForwards====='
    do i = n, 0, -1
        write(*,'(f12.6)') P(i)
    end do
    write(*,*) " x_i        P_i"
    do j = 0, n
        write(*,*) x(j), sum([(P(i)*x(j)**(n-i), i = 0, n)])
    end do
end subroutine