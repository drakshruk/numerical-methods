program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: studNum, N, M, i, j
    real(8) :: hx, ht, x0, x1, t0, t1
    real(8) :: sgm, tmp, k, a, temp
    real(8), allocatable :: u(:,:), c(:), d(:)
    character(50) :: format

    x0 = 0d0;    x1 = 1d0
    t0 = 0d0;    t1 = 20d0
    studNum = 1
    ht = 2d-3; hx = 1d-1
    N = int((t1-t0)/ht + 5d-1); M = int((x1-x0)/hx + 5d-1)
    sgm = 5d-1
    k = ht/hx/hx
    a = -28d-1 + 5d-1*studNum

    allocate(u(0:M,0:N), c(0:M), d(0:M))
    c(0:M) = 0d0;    d(0:M) = 0d0
    u(0:M,0) = 2d0*a
    u(0,0:N) = 3d0*a
    u(M,0:N) = a

    do j = 0, N
        tmp = 1d0+2d0*sgm*k
        c(1) = sgm * k / tmp;
        d(1) = (u(1,j) + (1d0 - sgm) * k * (u(2,j) - 2d0 * u(1,j) + u(0,j)) + sgm * k * 3d0 * a) / tmp;
        do i = 2, M-1
            tmp = 1d0 + 2d0 * sgm * k - sgm * k * c(i - 1);
            c(i) = sgm * k / tmp;
            d(i) = (u(i,j) + (1d0 - sgm) * k * (u(i + 1,j) - 2d0 * u(i,j) + u(i-1,j)) + sgm * k * d(i - 1)) / tmp;
        end do

        do i = M-1, 1, -1
            u(i,j+1) = c(i)*u(i+1,j+1) + d(i)
        end do
        tmp = 0d0
    end do

    format = '((f12.6, a3, f12.6, a3, a30, a3, a30))'
    do i = 0, M
        temp = abs(uReal(x0+hx*i)-u(i,N))
        write(*,format) u(i,N), ' & ', uReal(x0+hx*i), ' & ', to_lx(temp), ' & ', to_lx(temp/abs(uReal(x0+hx*i)))
    end do
    write(*,*)
    pause
    contains
    real(8) function uReal(x)
        real(8) :: x
        uReal = a*(3d0-2d0*x)
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