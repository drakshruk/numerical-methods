subroutine gauss(AB, n, X)
    implicit none
    integer n, i, j
    real(8) :: AB(0:n,0:n+1), X(0:n)
    do i = 0, n - 1
        ab(i,i:n+1) = ab(i,i:n+1) / ab(i,i)
        do j = i + 1, n
            ab(j,i:n+1) = ab(j,i:n+1) - ab(j,i) * ab(i,i:n+1)
        end do
    end do
    ab(n,n:n+1) = ab(n,n:n+1) / ab(n, n)

    do i = n, 1, -1
        do j = i-1, 0, -1
            ab(j,i:n+1) = ab(j,i:n+1) - ab(j,i) * ab(i,i:n+1)
        end do
        X(i) = ab(i,n+1)
    end do
    X(0) = ab(0,n+1)

end subroutine gauss

real(8) function f(x)
    implicit none
    real(8) x
    f = exp(2*sqrt(x) + 1./2./sqrt(x))
end function

real(8) function fder(x)
    implicit none
    real(8) x
    fder = (1./sqrt(x) - 1./4./sqrt(x**3))*exp(2*sqrt(x) + 1./2./sqrt(x))
end function

subroutine MNK(X, n, der, x0)
    implicit none
    integer n, i, j
    real(8) :: x(0:n), c(0:n), fi(0:n,0:n+1), der, x0, f
    do j = 0, n
        fi(j,0:n) = [((x(i)-x(0))**j, i = 0,n)]
        fi(j,n+1) = j*(x0-x(0))**(j-1)
    end do
    call gauss(fi,n,c)
    write(*,*)
    do i = 0, n
        write(*,'(f12.8)') c(i)
    end do
    der = sum([(c(i)*f(x(i)), i = 0, n)])
end subroutine