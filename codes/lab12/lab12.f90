program lab12
    implicit none
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    real(8) :: A1(1:5,1:5), B1(1:5), A2(1:5,1:5), B2(1:5), tmp(1:3,1:3)
    real(8) ::  X1(1:5), X2(1:5), X2Sei(1:5), X2Iter(1:5)
    real(8) :: a, eps
    integer :: studNum, i, j
    studNum = 1
    eps = 1d-3
    a = 0.01 + 0.002*studNum
    A1(1,1:5) = [a*a, 11d-1, 12d-1, 13d-1, 14d-1]
    A1(2,1:5) = [1.4*a, 2.*a*a, 1.1+a, 1.2+a, 1.3+a]
    A1(3,1:5) = [1.3*a, 1.4*a, 3.*a*a, 1.1+a, 1.2+2*a]
    A1(4,1:5) = [1.2*a, 1.3*a, 1.4*a, 4.*a*a, 1.1+3*a]
    A1(5,1:5) = [1.1*a, 1.2*a, 1.3*a, 1.4*a, 5.*a*a]
    B1(1:3) = [3.21*a*a + 1.934, 19.58*a*a - 1.326*a - 8.253, 16.74*a*a - 0.631*a - 13.251]
    B1(4:5) = [-17.16*a*a + 3.061*a - 7.821, -35.55*a*a + 16.527*a]
    call gauss(A1, B1, X1, 5)
    ! write(*,*) 'a1'
    ! do i = 1, 5
    !     write(*,'(*(f12.6))') A1(1:5,i)
    ! end do
    ! write(*,*) 'b1'
    ! write(*,'(f12.6)') B1
    ! write(*,*) 'x1'
    ! write(*,'(f12.6)') X1

    a = 0.02 + 0.001*studNum
    A2(1,1:5) = [10d0,      a + 0.3,            a + 0.2,            a + 0.1,            a + 1d0]
    A2(2,1:5) = [10.*a-3.,  a*a+9.91,           a*a+0.9*a+0.24,     a*a+0.8*a+0.17,     a*a+1.7*a+1.7]
    A2(3,1:5) = [10.*a-2.,  a*a+10.1*a-3.06,    2.*a*a+9.87,        2.*a*a+0.8*a+0.22,  2.*a*a+3.5*a+0.22]
    A2(4,1:5) = [10.*a-1.,  a*a+10.2*a-2.03,    2.*a*a+10.2*a-3.08, 3.*a*a+9.86,        3.*a*a+6.4*a+2.6]
    A2(5,1:5) = [10.*a-10., a*a+9.3*a-20.3,     2.*a*a+7.5*a-30.8,  3.*a*a+4.6*a-41.4,  4.*a*a-20.]
    B2(1:3) = [35.07*a + 25.07, 35.07*a*a + 49.379*a + 40.122, 69.9*a*a + 88.55*a + 63.1691]
    B2(4:5) = [103.2*a*a + 155.906*a + 169.1416, 124.38*a*a + 71.759*a - 979.664]

    ! write(*,'(*(f12.6))') transpose(A2)
    ! write(*,'(*(f12.6))') B2

    A2(5,1:5) = A2(5,1:5)+ 2.9*A2(1,1:5) + 3.8*A2(2,1:5) + 4.*A2(3,1:5)+ A2(4,1:5)*4.
    B2(5) = B2(5) + 2.9*B2(1) + 3.8*B2(2) + 4.*B2(3) + 4.*B2(4)
    ! write(*,*) 'a2'
    ! do i = 1, 5
    !     write(*,'(*(f12.6))') A2(1:5,i)
    ! end do
    ! write(*,*) 'b2'
    ! write(*,'(f12.6)') b2

    call gauss(A2, B2, X2, 5)
    call iterations(A2,B2,X2Iter,5,eps)
    call seidel(A2,B2,X2Sei,5,eps)

    ! write(*,*) 'x2 iter'
    ! write(*,'(*(f12.6))') X2Iter
    ! write(*,'(*(f12.6))') X2
    ! do i = 1, 5
    ! write(*,*) to_lx(abs(X2(i)-X2Iter(i)))
    ! end do

    ! write(*,*) 'x2 sei'
    ! write(*,'(*(f12.6))') X2Sei
    ! write(*,'(*(f12.6))') X2
    ! do i = 1, 5
    ! write(*,*) to_lx(abs(X2(i)-X2Sei(i)))
    ! end do

    tmp = 0d0
    tmp(1,1) = 1d0
    tmp(2,3) = 1d0
    tmp(3,2) = 1d0
    call gauss(tmp,[(0d0,i=1,3)],[(0d0,i=1,3)],3)
    write(*,*)
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