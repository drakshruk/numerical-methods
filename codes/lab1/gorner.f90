subroutine gorner(x, P, n)
    implicit none
    integer n, m, i, j
    real(8) :: x, P(0:n), bn, y
    real(8) t, t1
    m = 1000000
    write(*,*) '============gorner=========='
    call cpu_time(t)
    do j = 0, m
        bn = P(0)
        do i = 1, n
            bn = bn * x + P(i)
        end do
    end do
    call cpu_time(t1)
    write(*,*) 'Gorner time =', t1 - t
    call cpu_time(t1)
    do j = 0, m
        y = 0
        do i = 0, n
            y = y + P(i) * x**(n - i)
        end do
    end do
    call cpu_time(t)
    write(*,*) 'Regular time = ', t - t1
end subroutine

function to_lx(number) result(output_string)
    implicit none
    real(kind=8), intent(in) :: number
    character(len=30) :: output_string
    character(len=30) :: buffer
    integer :: e_position, dot_position
    character(len=20) :: mantissa, exponent
    logical :: negative_exponent
    character(len=20) :: formatted_mantissa
    integer :: i
    
    write(buffer, '(ES23.15E3)') number
    
    ! Find the position of 'E' in the string
    e_position = index(buffer, 'E')
    if (e_position == 0) then
        e_position = index(buffer, 'e')  ! Try lowercase e if uppercase not found
    end if
    
    ! If no exponent found, return the number as is (with 6 decimal digits)
    if (e_position == 0) then
        write(output_string, '(F0.6)') number
        output_string = trim(adjustl(output_string))
        return
    end if
    
    ! Extract mantissa and exponent parts
    mantissa = buffer(1:e_position-1)
    exponent = buffer(e_position+1:)
    
    ! Find the decimal point in the mantissa
    dot_position = index(mantissa, '.')
    
    ! Format mantissa to have exactly 6 digits after decimal point
    if (dot_position > 0) then
        ! Keep only 6 digits after decimal
        if (len_trim(mantissa) > dot_position + 6) then
            mantissa = mantissa(1:dot_position+6)
        end if
        mantissa = adjustl(mantissa)
        ! If we removed all decimals, keep the decimal point
        if (mantissa(len_trim(mantissa):len_trim(mantissa)) == '.') then
            mantissa = trim(mantissa) // '0'
        end if
    end if
    
    ! Remove leading zeros from exponent
    exponent = adjustl(exponent)
    
    ! Check if exponent is negative
    negative_exponent = exponent(1:1) == '-'
    if (negative_exponent) then
        exponent = exponent(2:)
    end if
    
    ! Remove any leading '+' in exponent
    if (exponent(1:1) == '+') then
        exponent = exponent(2:)
    end if
    
    ! Construct the output string
    if (negative_exponent) then
        output_string = '$' // trim(adjustl(mantissa)) // '\cdot 10^{-' // trim(adjustl(exponent)) // '}$'
    else
        output_string = '$' // trim(adjustl(mantissa)) // '\cdot 10^{' // trim(adjustl(exponent)) // '}$'
    end if
end function to_lx