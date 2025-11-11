program main
    real(8), allocatable :: X(:), Y(:), Aa(:,:)
    real(8) h, x0, a, d1, d2
    real(8) dl1,dl2,dn1,dn2,ds1,ds2
    integer n
    n = 19
    allocate(X(0:n), Y(0:n),Aa(0:n,0:n))
    h = 1e-2
    x0 = 7e-2
    a = 1e-2
    X = [(0.01*i, i = 0,n)]
    Y = [(sqrt(x(i)+a), i = 0,n)]
    d1 = 1./2./sqrt(x0+a)
    d2 = -1./4./sqrt(x0+a)**3
    write(*,*) '0.07l:'
    dl1 = derL(Y, X(0), x0, h, n)
    dl2 = derL2(Y(4:10), X(4:10), x0, 6)
    write(*,*) '0.07l:'
    write(*,('(6(f14.10))')) dl1, dl2, abs(d1-dl1), abs(d2-dl2), abs(d1-dl1)/d1, abs(d2-dl2)/abs(d2)
    dn1 = derN(Y, X(0), x0, h, n)
    dn2 = derN2(Y, X(0), x0, h, n)
    write(*,*) '0.07n:'
    write(*,('(6(f14.10))')) dn1, dn2, abs(d1-dn1), abs(d2-dn2), abs(d1-dn1)/d1, abs(d2-dn2)/abs(d2)
    ds1 = derS(Y(5:9), X(5:9))
    ds2 = derS2(Y(5:9), X(5:9))
    write(*,*) '0.07s:'
    write(*,('(6(f14.10))')) ds1, ds2, abs(d1-ds1), abs(d2-ds2), abs(d1-ds1)/d1, abs(d2-ds2)/abs(d2)

    x0 = x0*2.
    d1 = 1./2./sqrt(x0+a)
    d2 = -1./4./sqrt(x0+a)**3
    write(*,*) '0.14:'
    dl1 = derL(Y, X(0), x0, h, n)
    dl2 = derL2(Y(11:17), X(11:17), x0, 6)
    write(*,*) '0.14l:'
    write(*,('(6(f14.10))')) dl1, dl2, abs(d1-dl1), abs(d2-dl2), abs(d1-dl1)/d1, abs(d2-dl2)/abs(d2)
    dn1 = derN(Y, X(0), x0, h, n)
    dn2 = derN2(Y, X(0), x0, h, n)
    write(*,*) '0.14n:'
    write(*,('(6(f14.10))')) dn1, dn2, abs(d1-dn1), abs(d2-dn2), abs(d1-dn1)/d1, abs(d2-dn2)/abs(d2)
    ds1 = derS(Y(12:16), X(12:16))
    ds2 = derS2(Y(12:16), X(12:16))
    write(*,*) '0.14s:'
    write(*,('(6(f14.10))')) ds1, ds2, abs(d1-ds1), abs(d2-ds2), abs(d1-ds1)/d1, abs(d2-ds2)/abs(d2)

    pause
end program