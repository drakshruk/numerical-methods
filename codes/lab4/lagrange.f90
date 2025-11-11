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
    integer n, i
    real(8) :: X(0:n), Y(0:n), P(0:n), A(0:n)
    P(0:n) = 0d0
    do i = 0, n
        call L(A, X, n, i)
        P(0:n) = P(0:n) + A(0:n) * Y(i)
    end do
end subroutine

real(8) function polyValue(x, n, A)
    implicit none
    integer n, i
    real(8) :: x, A(0:n)
    polyValue = 0d0
    do i = 0, n
        polyValue = polyValue + A(i)*x**(n-i)
    end do
end function polyValue

subroutine polyDer(P, n, pDer)
    implicit none
    integer n, i
    real(8) P(0:n), pDer(0:n)
    pDer(1:n) = [(P(i)*(n-i), i = 0, n-1)]
    pDer(0) = 0d0
end subroutine
