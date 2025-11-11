program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: studNum
    real(8) :: x(0:10), y(0:10), h, b, tmp
    character(60) :: format
    studNum = 1
    h = 1d-1
    x = [(0d0 + h*i, i = 0,10)]
    b = 1d0 + 0.001d0*studNum
    y(0) = 1d0/(b-2d0)
    y(1) = 1d0/(b-2d0)
    y(1) = y(1) + (b-1d0)/(b-2d0)*h + (2d0*b-3d0)/(b-2d0)*h*h/2d0 + (3d0*b-5d0)/(b-2d0)*h**3/6d0 + (4d0*b-7d0)/(b-2d0)/24d0*h**4

    do i = 2, 10
        y(i) = y(i-1)*(24d0+10d0*h*h*(2d0*b-3d0))/(12d0-h*h*(2d0*b-3d0)) - y(i-2) 
        y(i) = y(i) - 2d0*h*h*(b-2d0)*(x(i)*exp(x(i)) + 10d0*x(i-1)*exp(x(i-1)) + x(i-2)*exp(x(i-2)))/(12d0-h*h*(2d0*b-3d0))
    end do

    format = '(f12.6, a3, f12.6, a3, f12.6, a3, a30, a3, a30)'
    tmp = abs(rY(x(1))-y(1))
    write(*,format) x(1)," & ", rY(x(1))," & ", y(1)," & ", to_lx(tmp),"&", to_lx(tmp/abs(rY(x(1))))
    tmp = abs(rY(x(5))-y(5))
    write(*,format) x(5)," & ", rY(x(5))," & ", y(5)," & ", to_lx(tmp),"&", to_lx(tmp/abs(rY(x(5))))
    tmp = abs(rY(x(10))-y(10))
    write(*,format) x(10)," & ", rY(x(10))," & ", y(10)," & ", to_lx(tmp),"&", to_lx(tmp/abs(rY(x(10))))
    pause
    contains
    real(8) function rY(x)
        real(8) :: x, b
        b = 1d0 + 1d-3*studNum
        rY = (x+1d0/(b-2d0))*exp(x)
    end function
end program

function to_lx(number) result(output_string)
    implicit none
    real(kind=8), intent(in) :: number
    character(len=50) :: output_string
    character(len=50) :: buffer
    integer :: e_position, dot_position
    character(len=20) :: mantissa, exponent
    logical :: negative_exponent
    character(len=20) :: formatted_mantissa
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
