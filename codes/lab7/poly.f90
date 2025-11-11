recursive subroutine Lezh(A, n, k)
    integer n, k
    real(8) :: A(0:n), B(0:n)
    A(0:n) = 0d0
    if (k == 0) then
        A(0:n) = 0d0
        A(0) = 1d0
    else if(k==1) then
        A(0:n) = 0d0
        A(1) = 1d0
    else
        call Lezh(A,n,k-2)
        call Lezh(B,n,k-1)
        B(1:k+1) = B(0:k) * (2*k - 1)/(k)
        B(0) = 0d0
        A(0:n) = A(0:n) * (k - 1) / k
        A(0:n) = B(0:n) - A(0:n)
    end if
end subroutine

subroutine lezhCoef(X,n,C)
    integer n
    real(8) :: X(1:n), C(1:n), A(0:n), derivPolyValue
    call Lezh(A,n,n)
    C(1:n) = [(2/((1-x(i)*x(i))*derivPolyValue(x(i),n,A)**2), i = 1,n)]
end subroutine

real(8) function lezhFunc(x, b)
    real(8) :: x, b
    lezhFunc = 2/sqrt((3+x)**3)/sqrt(12+4*x+2*b)
end function lezhFunc

recursive subroutine Cheb(A, n, k)
    integer n, k
    real(8) :: A(0:n), B(0:n)
    if (k == 0) then
        A(:) = 0d0
        A(0) = 1d0
    else if(k==1) then
        A(:) = 0d0
        A(1) = 1d0
    else
        call Cheb(A,n,k-2)
        call Cheb(B,n,k-1)
        B(1:k+1) = B(0:k) * 2d0
        B(0) = 0
        A(0:n) = B(0:n) - A(0:n)
    end if
end subroutine

subroutine chebCoef(x,n,C)
    integer n
    real(8) :: X(1:n), C(1:n)
    C(1:n) = [(4d0*atan(1d0)/n, i = 1,n)]
end subroutine

real(8) function chebFunc(x, b)
    real(8) :: x, b
    chebFunc = sqrt((1d0-x*x)/(3d0+x-2d-2*b)/(5d0+x))
end function chebFunc

recursive subroutine Lag(A,n,k)
    integer n, k
    real(8) :: A(0:n), B(0:n)
    if (k == 0) then
        A(0:n) = 0d0
        A(0) = 1d0
    else if(k==1) then
        A(0:n) = 0d0
        A(0) = 1d0
        A(1) = -1d0
    else
        call Lag(A,n,k-2)
        call Lag(B,n,k-1)
        A(0:n) = A(0:n)*(-k+1)/(k)
        A(0:n) = A(0:n) + B(0:n)*(2*k-1)/(k)
        B(1:n) = B(0:n-1)/k
        B(0) = 0d0
        A(0:n) = A(0:n) - B(0:n)
    end if
end subroutine

subroutine lagCoef(x,n,C)
    integer n
    real(8) :: X(1:n), C(1:n), A(0:n+1), polyValue
    call Lag(A,n+1,n+1)
    C(1:n) = [(x(i)/(n+1)**2/polyValue(x(i),n+1,A)**2, i = 1,n)]
end subroutine

real(8) function lagFunc(x, b)
    real(8) :: x, b
    lagFunc = exp(b*x)/sqrt(x)
end function lagFunc

recursive subroutine Herm(A,n,k)
    integer n, k
    real(8) :: A(0:n), B(0:n)
    if (k == 0) then
        A(0:n) = 0d0
        A(0) = 1d0
    else if(k==1) then
        A(0:n) = 0d0
        A(1) = 2d0 !если физич то  = 2
    else
        call Herm(A,n,k-2)
        call Herm(B,n,k-1)
        B(1:n) = B(0:n-1)
        B(0) = 0d0
        A(0:n) = (B(0:n) - A(0:n)*(k-1))*2 !если физич то всё *2
    end if
end subroutine

subroutine hermCoef(x,n,C)
    integer n
    real(8) :: X(1:n), C(1:n), A(0:n), polyValue
    call Herm(A,n,n-1)
    C(1:n) = [(2**(n-1)*gamma(n+1d0)*sqrt(4*atan(1d0))/(n*polyValue(x(i),n,A))**2, i = 1,n)]
end subroutine

real(8) function hermFunc(x, b)
    real(8) :: x, b
    hermFunc = 1d0/sqrt(4d0*atan(1d0))*cos(2d0*b*x/(b+1d-2))
end function hermFunc

real(8) function polyValue(x, n, A)
    integer n
    real(8) :: x, A(0:n)
    polyValue = 0d0
    do i = 0, n
        polyValue = polyValue + A(i)*x**i
    end do
end function polyValue

real(8) function derivPolyValue(x,n,A)
    integer n
    real(8) :: x, A(0:n)
    derivPolyValue = 0d0
    do i = 1, n
        derivPolyValue = derivPolyValue + i*A(i)*x**(i-1)
    end do
end function derivPolyValue