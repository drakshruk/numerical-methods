subroutine seidel(A1, B1, X1, n, eps)
    implicit none
    integer :: n, i, j, counter
    real(8), intent(in) :: A1(1:n,1:n), B1(1:n)
    real(8) :: X1(1:n), X(1:n), A(1:n,1:n), B(1:n), eps
    X1(1:n) = 0d0
    counter = 0
    do i = 1, n
        A(i,1:n) = [(-A1(i,j) / A1(i,i), j = 1,n)]
        A(i,i) = 0d0
    end do
    B(1:n) = [(B1(i)/A1(i,i), i = 1,n)]
    X1(1:n) = B(1:n)
    do while(norm2(X1(1:n) - X(1:n)) > eps)
        X(1:n) = X1(1:n)
        X1(1) = B(1) + sum([(A(1,j)*X(j), j = 2, n)])
        do i = 2, n
            X1(i) = B(i) + sum([(A(i,j)*X(j), j = i+1, n)]) + sum([(A(i,j)*X1(j), j = 1, i-1)])
        end do
        counter = counter + 1
    end do
    write(*,*) "Seidel counter:"
    write(*,*) counter
end subroutine