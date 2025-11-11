subroutine mk(n)
    implicit none
    integer n, nIN, i
    real(8) :: X, Y, Z, intValue(1:3), f1mk, f2mk, f3mk, tmp, trueInt(1:3)
    trueInt = [0.261799, 0.500000, 0.034074]
    intValue(1:3) = 0d0
    nIN = 0
    do i = 1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(Y)
        call RANDOM_NUMBER(Z)
        if(X >= Y) intValue(1) = intValue(1) + f1mk(x,y)
    end do
    do i = 1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(Y)
        call RANDOM_NUMBER(Z)
        if(X <= Y**2) intValue(2) = intValue(2) + f2mk(x,y)
    end do
    do i = 1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(Y)
        call RANDOM_NUMBER(Z)
        if(X+Y+Z <= 1) intValue(3) = intValue(3) + f3mk(x,y,z)
    end do
    intValue(1:3) = intValue(1:3) / n
    write(*,'(*(f12.6))') trueInt
    write(*,'(*(f12.6))') intValue
    write(*,*) abs(intValue-trueInt)
    write(*,*) abs(intValue-trueInt)/trueInt
end subroutine

real(8) function f1mk(x,y)
    real(8) x,y
    f1mk = sqrt(x*x-y*y)
end function

real(8) function f2mk(x,y)
    real(8) x,y
    f2mk = exp(x/y)
end function

real(8) function f3mk(x,y,z)
    real(8) x,y,z
    f3mk = 1/(x+y+z+1.)**3
end function