program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: studNum, n
    real(8), allocatable :: y1(:), c1(:), y2(:), c2(:), y3(:), x1(:), x2(:), x3(:)
    real(8) :: h, a, tmp, tmp1, tmp2, tmp3
    character(50) :: format

    format = '(f12.6, a3, f12.6, a3, f12.6, a3, a30, a3, a30)'
    studNum = 1
    n = 10
    h = 1d-2
    allocate(y1(0:10), c1(0:5), x1(0:10))
    x1 = [(2d0+h*i, i = 0,10)]
    a = 8d-2 + 2d-2*studNum

    c1(0) = 2d0
    c1(1) = 4d0
    c1(2) = -(5d0*a+12d0)/20d0
    c1(3) = (96d0-35d0*a)/600d0
    c1(4) = (125d0*a*a+1505d0*a-1248d0)/24000d0
    c1(5) = (11232d0-22620d0*a-16250d0*a*a)/600000d0
    write(*,*) 'task1'
    do j = 0, 10
        y1(j) = sum([(c1(i)*(x1(j)-2d0)**i/gamma(i+1d0), i = 0,5)])
        tmp = abs(rY1(x1(j), studNum)-y1(j))
        write(*,format) x1(j),"&",rY1(x1(j), studNum),"&",y1(j),"&",to_lx(tmp),"&",to_lx(tmp/rY1(x1(j), studNum))
    end do

    allocate(c2(0:5), y2(0:10), x2(0:10))
    a = -3d0 + 1d-2*studNum
    tmp = t(2d0,studNum)
    tmp1 = dert(2d0,studNum)
    tmp2 = der2t(2d0,studNum)
    tmp3 = 18d0*a+12d0
    x2 = [(2d0+h*i, i = 0,10)]

    c2(0) = 1/3d0
    c2(1) = -1/3d0
    c2(2) = -(2d0*tmp*c2(1)+tmp2*c2(0))/tmp
    c2(3) = -((3d0*tmp2*tmp-6d0*tmp1*tmp1)*c2(1)+(tmp*tmp3-3d0*tmp1*tmp2)*c2(0))/tmp/tmp
    c2(4) = (4d0*tmp*tmp*tmp3-24d0*tmp*tmp1*tmp2+24d0*tmp1**3)*c2(1)
    c2(4) = -(c2(4)+ (-4d0*tmp*tmp1*tmp3-6d0*tmp*tmp2*tmp2+12d0*tmp1*tmp1*tmp2)*c2(0))/tmp**3
    c2(5) = ((20d0*tmp*tmp*tmp2-20*tmp*tmp1*tmp1)*tmp3-60d0*tmp*tmp1*tmp2*tmp2+60d0*tmp1**3*tmp2)*c2(0)
    c2(5) = (c2(5) + (40d0*tmp*tmp*tmp1*tmp3+30d0*tmp*tmp*tmp2*tmp2-180d0*tmp*tmp1*tmp1*tmp2+120d0*tmp1**4)*c2(1))/c2(0)**4
    
    write(*,*) 'task2'
    do j = 0, 10
        y2(j) = sum([(c2(i)*(x2(j)-2d0)**i/gamma(i+1d0), i = 0,5)])
        tmp = abs(rY2(x1(j), studNum)-y2(j))
        write(*,format) x2(j),"&",rY2(x1(j), studNum),"&",y2(j),"&",to_lx(tmp),"&",to_lx(tmp/rY2(x1(j), studNum))
    end do

    write(*,*) 'task3'
    call num3
    pause
    contains
    
    real(8) function rY1(x, studNum)
        implicit none
        real(8) :: x, c1, c2, tmp
        integer :: studNum
        tmp = sqrt(2d0*studNum+4d0)*log(2d0)/10d0
        c1 = 2d0*sqrt(2d0*studNum+4d0)*sin(tmp)+76d0*cos(tmp)
        c2 = 2d0*sqrt(2d0*studNum+4d0)*cos(tmp)-76d0*sin(tmp)
        
        tmp = 2d0**(2d-1)*sqrt(2d0*studNum+4d0)
        c1 = c1/tmp
        c2 = c2/tmp

        tmp = sqrt(4d-2+2d-2*studNum)*log(x)
        rY1 = x**(2d-1)*(c1*sin(tmp) + c2*cos(tmp))
    end function

    real(8) function t(x, studNum)
        implicit none
        real(8) :: x, a
        integer :: studNum
        a = -3d0 + 1d-2*studNum
        t = (3d0*a+2d0)*x*x*x-(a+3d0)*x*x+(a+3d0)*x+(2d0*a-1d0)
    end function
    
    real(8) function dert(x, studNum)
        implicit none
        real(8) :: x, a
        integer :: studNum
        a = -3d0 + 1d-2*studNum
        dert = (3d0*a+2d0)*x*x*3d0-(a+3d0)*x*2d0+(a+3d0)
    end function
    
    real(8) function der2t(x, studNum)
        implicit none
        real(8) :: x, a
        integer :: studNum
        a = -3d0 + 1d-2*studNum
        der2t = (3d0*a+2d0)*x*6d0-(a+3d0)*2d0
    end function

    real(8) function rY2(x, studNum)
        implicit none
        real(8) :: x, c1, c2
        integer :: studNum
        c1 = ((3d0*a+2d0)*12d0-(a+3d0)*4d0+(a+3d0))/3d0 - t(2d0, studNum)/3d0
        c2 = t(2d0, studNum)/3d0 - 2d0*c1
        
        rY2 = (c1*x+c2)/t(x,studNum)
    end function
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