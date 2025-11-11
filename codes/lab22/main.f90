program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    real(8) :: a, hx, ht, x0, x1, t0, t1, pi, k, tmp
    integer :: studNum, N, M
    real(8), allocatable :: u(:,:)
    character(50) :: format
    
    pi = 3.1415926535d0
    studNum = 1
    a = 2d0 + 1d-2*studNum

    x0 = 0d0; x1 = 1d0
    t0 = 0d0; t1 = a/2d0+1d0
    ht = 1d-2; hx = 1d-1
    N = int((t1-t0)/ht + 5d-1); M = int((x1-x0)/hx + 5d-1)
    k = ht/a/hx

    allocate(u(0:N,0:M))
    u(0,0:M) = [(nUsl1(x0+hx*i), i = 0,M)]
    u(0:N,0) = [(gUsl1(t0+ht*i), i = 0,N)]
    u(0:N,M) = [(gUsl2(t0+ht*i), i = 0,N)]

    u(1,0:M) = [(u(0,i) - ht*ht*pi*pi/a/a*sin(pi*(x0+hx*i)), i = 0,M)]

    do i = 1, N-1
        do j = 1, M-1
            u(i+1,j) = k*k*(u(i,j+1)-2d0*u(i,j)+u(i,j-1)) + 2d0*u(i,j) - u(i-1,j)
        end do
    end do
    
    format = '((f12.6, a3, f12.6, a3, f12.6, a3 a30, a3, a30))'
    do i = 0, M
        tmp = uReal(x0+i*hx, t1)
        write(*,format) x0+i*hx," & ",tmp," & ",u(N,i)," & ",to_lx(abs(tmp-u(N,i)))," & ",to_lx(abs(tmp-u(N,i))/abs(tmp))
    end do

    pause
    contains
    real(8) function nUsl1(x)
        real(8) :: x
        nUsl1 = 2d0*sin(pi*x)
    end function
    
    real(8) function nUsl2(x)
        real(8) :: x
        nUsl2 = 0d0
    end function
    
    real(8) function gUsl1(t)
        real(8) :: t
        gUsl1 = 0d0
    end function
    
    real(8) function gUsl2(t)
        real(8) :: t
        gUsl2 = 0d0
    end function

    real(8) function uReal(x,t)
        real(8) :: x, t
        uReal = 2d0*cos(pi/a*t)*sin(pi*x)
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