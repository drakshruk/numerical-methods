program main
    implicit none
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    real(8), allocatable :: x(:)
    real(8) :: der, x0, fder, h, yRealDer, runge, r, yder, eps, tmp
    integer n, i
    n = 4
    allocate(x(0:n))
    x(0:n) = [(0.25 + (0.05+0.001)*i, i = 0, n)]
    x0 = 3d-1 + 0.001
    write(*,*) 'mnk'
    write(*,'(f12.8)') x0
    call MNK(X,n,der,x0)
    write(*,*)
    write(*,'(4(f12.8))') x0, der, fder(x0), abs(der - fder(x0)), abs(der - fder(x0))/abs(fder(x0))
    x0 = 0.4 + 0.003
    write(*,'(f12.8)') x0
    call MNK(X,n,der,x0)
    write(*,*)
    write(*,'(4(f12.8))') x0, der, fder(x0), abs(der - fder(x0)), abs(der - fder(x0))/abs(fder(x0))
    write(*,*) 'runge'
    x0 = 4*atan(1.)/200
    h = 4*atan(1.)/800
    eps = 1d-4
    r = 0.1
    do while(abs(yRealDer(x0) - runge(x0,h,r,2)) >= eps)
        write(*,'(*(f12.6))') r
        tmp = abs(yRealDer(x0)-yder(x0,h*r))
        write(*,*) 'P'
        write(*,'(*(f12.6))') x0, yder(x0,h*r)!, to_lx(tmp), to_lx(tmp/abs(yRealDer(x0)))
        der = runge(x0,h,r,2)
        tmp = abs(yRealDer(x0)-der)
        write(*,*) 'Y'
        write(*,'(*(f12.6))') x0, der!, to_lx(tmp), to_lx(tmp/abs(yRealDer(x0)))
        r = r/10
    end do
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
