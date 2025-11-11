program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: studNum, n, i
    real(8) :: y0, x0, a, b
    real(8) :: h, eps, c, tmp
    real(8), allocatable :: y(:), x(:)
    character(50) :: format

    format = '(f12.3, a3, f12.6, a3, f12.6, a3, a30, a3, a30)'

    studNum = 1
    n = 100
    a = 1d0
    b = 2d0
    c = 1 + 1d-3*studNum
    h = (b-a)/n

    x0 = 1d0
    y0 = 1d0 + c
    eps = 1d-4

    allocate(y(0:n), x(0:n))
    call runge(y0,x0,h,f,n,y)
    x = [(x0+h*i, i=0,n)]
    do i = 3, n
        y(i) = y(i-1) + h/12d0*(5d0*f(x(i-3),y(i-3) - 16d0*f(x(i-2),y(i-2)) + 23d0*f(x(i-1),y(i-1))))
        tmp = (5*f(x(i-3),y(i-3))+15d0*f(x(i-2),y(i-2))+39d0*f(x(i-1),y(i-1))+13d0*f(x(i),y(i)))
        y(i) = (y(i-1)+y(i-2)+y(i-3))/3d0 + h/36d0*tmp
    end do

    do i = 0, n
        tmp = abs(y(i)-corF(x(i)))
        write(*,*) x(i),'&',corF(x(i)),'&', y(i),'&', tmp,'&',tmp/corF(x(i))
    end do
    
    do i = 0, n
        write(*,*) x0+h*i, corF(x0+h*i)
    end do
    pause
contains
    real(8) function f(x,y)
        real(8) :: x, y
        f = -x/4d0/c - 5*c/4d0/x/x + y*y/4d0/c
    end function
    
    real(8) function corF(x)
        real(8) :: x
        corF = (c + x*sqrt(x))/x
    end function
end program

subroutine runge(y0,x0,h,fun,n,y1)
    implicit none
    integer :: i, n
    real(8) :: y0, x0, h
    real(8), external :: fun
    real(8) :: y1(0:n), x1, k(1:4)

    y1(0) = y0
    x1 = x0

    do i = 1,n+1
        k(1) = h*fun(x1, y1(i-1))
        k(2) = h*fun(x1+h*5d-1, y1(i-1)+k(1)*5d-1)
        k(3) = h*fun(x1+h*5d-1, y1(i-1)+k(2)*5d-1)
        k(4) = h*fun(x1+h, y1(i-1)+k(3))
        x1 = x0 + h*i
        y1(i) = y1(i-1) + 1d0/6*(k(1) + 2*k(2) + 2*k(3) + k(4))
    end do
end subroutine

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
