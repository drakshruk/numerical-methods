recursive subroutine Cheb(A, n, k)
    integer n, k
    real(8) :: A(0:n), B(0:n)
    if (k == 0) then
        A(:) = 0
        A(0) = 1
    else if(k==1) then
        A(:) = 0
        A(1) = 1
    else
        call Cheb(A,n,k-2)
        call Cheb(B,n,k-1)
        B(1:k+1) = B(0:k) * 2
        B(0) = 0
        A(0:n) = B(0:n) - A(0:n)
    end if
end subroutine

subroutine calcChebCoefficients(m, pogr)
    implicit none
    real(8) :: pogr(1:5)
    integer m, ni, i, j
    real(8), allocatable :: poly(:), func(:), weight(:)
    real(8) :: ai(0:m), x0, h, a, al(0:m), ap(0:m), y, x
    interface chebPoly_interface
        subroutine Cheb(A, n, k)
            real(8) A(0:n)
            integer k, n
        end subroutine
    end interface chebPoly_interface
    interface lezhValue_interface
        real function polyValue(x, n, A)
            real(8) :: x, A(0:n)
            integer n
        end function polyValue
    end interface lezhValue_interface
    interface fPoly_interface
        real function fPoly(x,a)
            real(8) :: x, a
        end function fPoly
    end interface fPoly_interface
    interface integral_interface
        real function integral(a, b, N, f1, f2, weight)
            implicit none
            real(8) a, b, f1(0:), f2(0:), weight(0:)
            integer N
        end function integral
    end interface integral_interface
    ni = 100000
    allocate(weight(0:ni), func(0:ni), poly(0:ni))
    a = 1.5 - 0.02*9
    x0 = -1.
    h = 2./ni
    func = [(fPoly(x0 + i*h, a),i = 0,ni)]
    weight(0:ni) = [(1/sqrt(1-(x0+i*h)**2), i = 0, ni)]
    weight(0) = 1.
    weight(ni) = 1.
    ap(:) = 0.
    do i = 0, m-1
        call Cheb(al,m,i)
        poly = [(polyValue(x0+j*h, m, al), j = 0, ni)]
        Ai(i) = integral(-1d0, 1d0, ni, func, poly, weight) / integral(-1d0, 1d0, ni, poly, poly, weight)
        ap(0:m) = ap(0:m) + al(0:m) * ai(i)
    end do
    x = -1.5
    do j = 1, 5
        y = 0.
        x = x + 0.5
        do i = 0, m - 1
            y = y + ap(i)*x**i
        end do
        pogr(j) = y
        write(*,'(5(f12.8))') x+1.5, y, fPoly(x, a), abs(fPoly(x, a)-y), abs(fPoly(x, a)-y)/fPoly(x, a)
    end do
end subroutine