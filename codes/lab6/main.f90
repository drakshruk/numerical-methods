program main
    implicit none
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    real(8), allocatable :: Y(:), YC(:), y1(:), y2(:)
    real(8) :: intValue(1:4), integrals(1:7)
    real(8) :: trapez, rRect, cRect, lRect, simpson, newton, h, a, b, f1,f2,f3,f4, p, eps
    integer n, i, j, studNum, n1, n2
    studNum = 1
    eps = 1d-4
    intValue(1:4) = [33.3456735358, 0.00063296505, 0.2426475369, 1.9851362]
    n = 20
    allocate(Y(0:n), YC(0:n))
    !1!====================================================================
    b = 1
    a = 0d0
    h = (b-a)/n
    p = 0.969 + 0.001*studNum
    y = [(f1(a+h*i, p), i = 0,n)]
    yc =[(f1(((a+h*i)+(a+h*(i-1)))/2, p), i = 0,n)]
    integrals(1:4) = [intValue(1),0.5/(1-p)+simpson(Y,a,b,n),0.5/(1-p)+trapez(Y,a,b,n),0.5/(1-p)+cRect(yc,a,b,n)]
    integrals(5:7) = [0.5/(1-p)+lRect(Y,a,b,n),0.5/(1-p)+rRect(Y,a,b,n),0.5/(1-p)+newton(Y,a,b,n+1)]
    write(*,'(*(f12.8))') integrals
    do i = 1, 7
        do j = 1, 7
            write(*,'(*((f12.8),a5))') abs(integrals(j) - integrals(i)), '&'
        end do
        write(*,*) '\\ \hline'
    end do
    !2!====================================================================
    b = atan(1.)*2
    a = 0d0
    h = (b-a)/n
    p = 10.2 + 0.01*studNum
    y = [(f2(a+h*i, p), i = 0,n)]
    yc =[(f2(((a+h*i)+(a+h*(i-1)))/2, p), i = 0,n)]
    integrals(1:7) = [intValue(2),simpson(Y,a,b,n),trapez(Y,a,b,n),cRect(yc,a,b,n),lRect(Y,a,b,n),rRect(Y,a,b,n),newton(Y,a,b,n+1)]
    write(*,'(*(f12.8))') integrals
    do i = 1, 7
        do j = 1, 7
            write(*,'(*((f12.8),a5))') abs(integrals(j) - integrals(i)), '&'
        end do
        write(*,*) '\\ \hline'
    end do
    !3!====================================================================
    allocate(y1(0:10000), y2(0:10000))
    p = 1 + 0.02*studNum
    n1 = 4
    b = 1
    a = 1d-10
    h = (b-a)/n1
    y1(:) = 0.0
    y1(0:n1) = [(f3(a+h*i, p), i = 0,n1)]
    do while(abs(intValue(3)-simpson(y1(0:n1),a,b,n1))>eps)
        n1 = n1 + 2
        h = (b-a)/n1
        y1(0:n1) = [(f3(a+h*i, p), i = 0,n1)]
    end do
    write(*,'(f12.8)') h
    
    n2 = 2
    h = (b-a)/n2
    y2(:) = 0.0
    y2(0:n2) = [(f3(a+h*i, p), i = 0,n2)]
    do while(abs(intValue(3)-trapez(y2(0:n2),a,b,n2))>eps)
        n2 = n2 + 1
        h = (b-a)/n2
        y2(0:n2) = [(f3(a+h*i, p), i = 0,n2)]
    end do
    write(*,'(f12.8)') h
    write(*,'(*(f12.8))') simpson(y1,a,b,n1), trapez(y2,a,b,n2), intValue(3)
    write(*,'(*(f12.8))') abs(simpson(y1,a,b,n1)-intvalue(3)), abs(trapez(y2,a,b,n2)-intvalue(3))
    write(*,'(*(f12.8))') abs(simpson(y1,a,b,n1)-intvalue(3))/intValue(3), abs(trapez(y2,a,b,n2)-intvalue(3))/intValue(3)
    !4!====================================================================
    write(*,*)

    p = 3 + 0.01*studNum
    n1 = 4
    b = 1
    a = 0.0000000001
    h = (b-a)/n1
    y1(:) = 0.0
    y1(0:n1) = [(f4(a+h*i, p), i = 0,n1)]
    do while(abs(intValue(4)-simpson(y1(0:n1),a,b,n1))>eps)
        n1 = n1 + 2
        h = (b-a)/n1
        y1(0:n1) = [(f4(a+h*i, p), i = 0,n1)]
    end do
    write(*,'(f12.8)') h

    n2 = 2
    h = (b-a)/n2
    y2(:) = 0.0
    y2(0:n2) = [(f4(a+h*i, p), i = 0,n2)]
    do while(abs(intValue(4)-trapez(y2(0:n2),a,b,n2))>eps)
        n2 = n2 + 1
        h = (b-a)/n2
        y2(0:n2) = [(f4(a+h*i, p), i = 0,n2)]
    end do
    write(*,'(f12.8)') h
    write(*,'(*(f12.8))') simpson(y1,a,b,n1), trapez(y2,a,b,n2), intValue(4)
    write(*,'(*(f12.8))') abs(simpson(y1,a,b,n1)-intvalue(4)), abs(trapez(y2,a,b,n2)-intValue(4))
    write(*,'(*(f12.8))') abs(simpson(y1,a,b,n1)-intvalue(4))/intValue(4), abs(trapez(y2,a,b,n2)-intValue(4))/intValue(4)
    pause
end program


function to_lx(number) result(output_string)
    implicit none
    real(8), intent(in) :: number
    character(50) :: output_string
    character(50) :: buffer
    integer :: e_position, dot_position
    character(20) :: mantissa, exponent
    logical :: negative_exponent
    character(20) :: formatted_mantissa
    integer :: i
    
    write(buffer, '(ES23.15E3)') number
    
    e_position = index(buffer, 'E')
    if (e_position == 0) then
        e_position = index(buffer, 'e')
    end if
    
    if (e_position == 0) then
        write(output_string, '(F0.6)') number
        output_string = trim(adjustl(output_string))
        return
    end if
    
    mantissa = buffer(1:e_position-1)
    exponent = buffer(e_position+1:)
    
    dot_position = index(mantissa, '.')
    
    if (dot_position > 0) then
        if (len_trim(mantissa) > dot_position + 6) then
            mantissa = mantissa(1:dot_position+6)
        end if
        mantissa = adjustl(mantissa)
        if (mantissa(len_trim(mantissa):len_trim(mantissa)) == '.') then
            mantissa = trim(mantissa) // '0'
        end if
    end if
    
    exponent = adjustl(exponent)
    
    negative_exponent = exponent(1:1) == '-'
    if (negative_exponent) then
        exponent = exponent(2:)
    end if
    
    if (exponent(1:1) == '+') then
        exponent = exponent(2:)
    end if
    
    if (negative_exponent) then
        output_string = '$' // trim(adjustl(mantissa)) // '\cdot 10^{-' // trim(adjustl(exponent)) // '}$'
    else
        output_string = '$' // trim(adjustl(mantissa)) // '\cdot 10^{' // trim(adjustl(exponent)) // '}$'
    end if
end function to_lx
