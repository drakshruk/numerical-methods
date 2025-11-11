program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: n, i, studNum
    real(8), allocatable :: y(:), x(:), rY(:)
    real(8) :: mi, ni, Fi, Ag, Bg, alpha(0:1), beta(0:1), h, a
    real(8) :: y1(0:10), y2(0:10), tmp
    character(50) :: format

    format = '(f12.6, a3, f12.6, a3, f12.6, a3, a30, a3, a30)'
    studNum = 1
    a = 1d0 + 0.001d0*studNum
    h = 1d-2
    n = 100

    alpha = [ 1d0, 1d0 ]
    beta = [ 2d0, -1d0 ]
    Ag = a*(sin(1d0)+cos(1d0)) + 6d0 + cos(1d0)*5d-1
    Bg = a*(sin(2d0)+2*cos(2d0)) + 5d0/8*sin(2d0) - cos(2d0)/4d0-3d0/4d0

    allocate(y(0:n), x(0:n), rY(0:n))
    x = [(1d0 + h*i, i = 0,n)]
    ry = [(realY(x(i), studNum), i = 0,n)]
    call progonka(studNum, Ag, Bg, alpha, beta, h, n, p, q, f, y)

    tmp = abs(ry(i)-y(i))

    write(*,*) 'h = 1d-2'
    do i = 0, n, n/10
        write(*,format) x(i),"&", ry(i),"&", y(i),"&", to_lx(tmp),"&", to_lx(tmp/abs(ry(i)))
    end do

    y1 = [(y(n/10*i), i = 0, 10)]

    deallocate(x, y, rY)
    h = 1d-3
    n = 1000
    allocate(y(0:n), x(0:n), rY(0:n))
    x = [(1d0 + h*i, i = 0,n)]
    ry = [(realY(x(i), studNum), i = 0,n)]
    call progonka(studNum, Ag, Bg, alpha, beta, h, n, p, q, f, y)

    write(*,*) 'h = 1d-3'
    do i = 0, n, n/10
        write(*,format) x(i),"&", ry(i),"&", y(i),"&", to_lx(tmp),"&", to_lx(tmp/abs(ry(i)))
    end do

    y2 = [(y(n/10*i), i = 0,10)]
    do i = 0,10
        write(*,'(*(a30))') to_lx(abs(y1(i)-y2(i)))
    end do
    pause
contains

real(8) function p(x)
    real(8) :: x
    p = -1d0/x
end function

real(8) function q(x)
    real(8) :: x
    q = -8d0/x/x
end function

real(8) function f(x, studNum)
    real(8) :: x, a
    integer :: studNum
    a = 1d0 + 0.001d0*studNum
    f = -a/x/x*(x*x*cos(x)-x*sin(x)+8*cos(x)) - (x*x*sin(x) + 3*x*cos(x)+5*sin(x))/2/x/x/x
end function

real(8) function realY(x, studNum)
    real(8) :: x, a
    integer :: studNum
    a = 1d0 + 0.001d0*studNum
    realY = sin(x)/2d0/x + a*cos(x) + (2*a*sin(1d0) + 5)/5*x**4 - 1d0/x/x
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
