subroutine gauss(A1, B1, X, n)
    implicit none
    integer :: n,i,j,k,s
    real(8), intent(in) :: A1(1:n,1:n), B1(1:n)
    real(8) :: X(1:n), tmp, C(1:n)
    real(8) :: A(1:n,1:n+1), M(1:n,1:n)
    A(1:n,1:n) = A1(1:n,1:n)
    A(1:n,n+1) = B1(1:n)
    do i = 1, n - 1
        M(1:n,1:n) = 0d0
        if(abs(A(i,i)) < 1d-6) then
            do j = i, n
                tmp = 0d0
                if(abs(A(i,j)) > tmp) then 
                    tmp = A(i,j)
                    s = j
                end if
            end do
            if(abs(A(s,i)) > 1d-6) then
                C(1:n) = A(s,1:n)
                A(s,1:n) = A(i,1:n)
                A(i,1:n) = C(1:n)
            else
                write(*,*) 'matrix is zero'
                pause
                stop
            end if
        end if
        do j = 1,n
            M(j,j) = 1d0
        end do
        M(i:n,i) = [(-A(j,i)/A(i,i), j = i,n)]
        m(i,i) = 1d0
        A(1:n,1:n+1) = matmul(M(1:n,1:n),A(1:n,1:n+1))
    end do
    do j = 1,n
        A(j,1:n+1) = A(j,1:n+1) / A(j,j)
    end do
    do j = n, 2, -1
        do k = j-1, 1, -1
            A(k,j:n+1) = A(k,j:n+1) - A(k,j) * A(j,j:n+1)
        end do
        X(j) = A(j,n+1)
    end do
    X(1) = A(1,n+1)
end subroutine