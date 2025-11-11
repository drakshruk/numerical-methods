subroutine spDifMat(X, Y, n, A)
    implicit none
    integer n, i, j
    real(8) X(0:n), Y(0:n), A(0:n,0:n)
    do i = 0, n
        do j = 0, n
            if( j < n - i) then 
                if( i == 0) then
                    A(j,i) = (Y(j+1) - Y(j)) / (X(j+1) - X(j))
                else
                    A(j,i) = (A(j+1,i-1) - A(j, i-1)) / (X(j+i+1) - X(j))
                end if
            else
                A(j,i) = 0d0
            end if
        end do
    end do
end subroutine

subroutine NP(A, X, n, k)
    implicit none
    integer i, n, j, k
    real(8) A(0:n), X(0:n)
    real(8), allocatable :: c(:)
    allocate(c(0:n))
    A(:) = 0d0
    A(n) = 1d0
    do i = 0, k - 1
        c(0:n) = A(0:n)
        A(0:n-1) = A(1:n)
        A(n) = 0
        do j = 0, n
            A(j) = A(j) + c(j) * (-1) * x(i) 
        end do
    end do
    deallocate(c)
end subroutine

subroutine Newton(P, Y, X, n)
    implicit none
    integer n, i, j
    real(8) :: X(0:n), Y(0:n), P(0:n), A(0:n), C(0:n,0:n)
    call spDifMat(X,Y,n,C)
    write(*,*) '=====Newton====='
    P(:) = 0
    P(n) = Y(0)
    do i = 1, n
        call NP(A, X, n, i)
        P(0:n) = P(0:n) + A(0:n) * C(0,i-1)
    end do
    do i = n, 0, -1
        write(*,'(f12.6)') P(i)
    end do
end subroutine