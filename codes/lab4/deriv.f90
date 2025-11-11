real(8) function derL(Y, x0, x, h, n)
    implicit none
    integer i, n, k
    real(8) :: Y(0:n), x0, x, h, tmp, qN, q
    q = (x - x0) / h
    derL = 0d0
    do i = 0, n
        tmp = sum([(qN(n, q) / (q - i) / (q - k), k = 0, n)])
        tmp = tmp - qN(n, q) / (q - i) / (q - i)
        derL = derL + tmp/h* Y(i) / gamma(i+1d0) / gamma(n-i+1d0) * (-1)**(n-i)
    end do
end function derL

real(8) function derL2(Y, X, x0, n)
    implicit none
    integer n, i
    real(8) :: Y(0:n), X(0:n), P(0:n), p1(0:n), x0
    derL2 = 0d0
    call lagrange(P,Y,X,n)
    call polyDer(P ,n,p1)
    call polyDer(p1,n, p)
    do i = 0, n
        derL2 = derL2 + p(i)*x0**(n-i)
    end do
end function

real(8) function derN(Y, x0, x, h, n)
    implicit none
    integer n, i, j
    real(8) :: difM(0:n,0:n), Y(0:n), x0, x, h, q,  qN, tmp
    q = (x - x0) / h
    call difMat(Y, difM, n)
    derN = 0d0
    do i = 0, n
        tmp = sum([(qN(i-1,q) / (q-j), j = 0, i-1)])
        tmp = tmp * difM(i,0) / gamma(i+1d0)
        derN = derN + tmp/h
    end do
end function

real(8) function derN2(Y, x0, x, h, n)
    implicit none
    integer n, k, i, j
    real(8) :: Y(0:n), x0, x, h, difM(0:n,0:n), q,  qN, tmp
    q = (x - x0) / h
    call difMat(Y, difM, n)
    derN2 = 0d0
    do i = 2, n
        tmp = 0d0
        do j = 0, i-1
            do k = 0, i-1
                if(k /= j) then
                    tmp = tmp + qN(i-1,q)/(q-j)/(q-k)
                end if
            end do
        end do
        tmp = tmp * difM(i, 0) / gamma(i+1d0)
        derN2 = derN2 + tmp/h/h
    end do
end function

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

real(8) function qN(n, q)
    implicit none
    integer n, i
    real(8), intent(in) :: q
    qN = 1d0
    do i = 0, n
        qN = qN * (q - i)
    end do
end function qN

real(8) function derS(Y, X)
    implicit none
    integer i, n, k
    real(8) :: Y(-2:2), X(-2:2), difM(4,4), h
    difM(:,1) = y(-1:2) - y(-2:1)
    do i = 2, 4
        difM(1:5-i,i) = difM(2:6-i,i-1) - difM(1:5-i,i-1)
    end do
    h = x(1)-x(0)
    derS = 1/h*((difM(2,1)+difM(3,1))/2 - (difM(1,3)+difM(2,3))/12)
end function

real(8) function derS2(Y, X)
    implicit none
    integer i, n, k
    real(8) :: Y(-2:2), X(-2:2), difM(4,4), h
    difM(:,1) = y(-1:2) - y(-2:1)
    do i = 2, 4
        difM(1:5-i,i) = difM(2:6-i,i-1) - difM(1:5-i,i-1)
    end do
    h = x(1)-x(0)
    derS2 = 1/h/h*(difM(2,2) - difM(1,4)/12)
end function

subroutine stirlingTable(diff_table,Y,n)
    implicit none
    integer n, i, j
    real(8) :: diff_table(0:n,0:n), Y(0:n)

    diff_table(0:n,0) = y

    do i = 1, n-1
        do j = 0, n-i-1
            diff_table(j,i) = diff_table(j + 1,i - 1) - diff_table(j,i - 1)
        end do
    end do
end subroutine
