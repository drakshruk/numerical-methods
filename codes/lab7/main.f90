program main
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer n
    real(8) :: iVal(1:4)
    real(8) :: xLezh(1:7), xCheb(1:7), xLag(1:7), xHerm(1:7)
    real(8) :: cLezh(1:7), cCheb(1:7), cLag(1:7), cHerm(1:7)
    real(8) :: lezhFunc, chebFunc, lagFunc, hermFunc
    real(8) :: lezhB, chebB, lagB, hermB
    real(8) :: lezhF, chebF, lagF, hermF
    character(50) :: format
    format = '(f12.6, a3, f12.6, a3, a30, a3, a30)'
    n = 7
    iVal = [0.229315645309d0, 0.531589823524d0, 12.5331413712d0, 0.352989137746d0]
    lezhF = 0d0
    chebF = 0d0
    lagF  = 0d0
    hermF = 0d0
    lezhB = 1d0 + 0.01d0
    chebB = 1d0
    lagB  = 1d0 - 0.02d0
    hermB = -0.5d0 + 0.01d0
    xLezh = [-0.94911d0,-0.74153d0,-0.40585d0,0d0,0.40585d0,0.74153d0,0.94911d0]
    xCheb = [(cos((2d0*i-1d0)/n*2d0*atan(1d0)), i = 1,n)]
    xLag  = [0.19304d0, 1.02666d0, 2.56788d0, 4.90031d0, 8.18271d0, 12.73223d0, 19.39782d0]
    xHerm = [-2.65196d0, -1.67355d0, -0.81629d0, 0d0, 0.81629d0, 1.67355d0, 2.65196d0]
    call lezhCoef(xLezh,n,cLezh)
    call chebCoef(xCheb,n,cCheb)
    call lagCoef(xLag,n,cLag)
    call HermCoef(xHerm,n,cHerm)
    do i = 1, n
        lezhF = lezhF + cLezh(i)*lezhFunc(xLezh(i),lezhB)
        chebF = chebF + cCheb(i)*chebFunc(xcheb(i),chebB)
        lagF = lagF + cLag(i)*lagFunc(xlag(i),lagB)
        hermF = hermF + cHerm(i)*hermFunc(xherm(i),hermB)
    end do
    write(*,format) lezhF, ' & ', iVal(1), ' & ', to_lx(abs(iVal(1)-lezhF)), ' & ', to_lx(abs(iVal(1)-lezhF)/iVal(1)) 
    write(*,format) chebF, ' & ', iVal(2), ' & ', to_lx(abs(iVal(2)-chebF)), ' & ', to_lx(abs(iVal(2)-chebF)/iVal(2)) 
    write(*,format) lagF,  ' & ', iVal(3), ' & ', to_lx(abs(iVal(3)-lagF)),  ' & ', to_lx(abs(iVal(3)-lagF)/iVal(3)) 
    write(*,format) hermF, ' & ', iVal(4), ' & ', to_lx(abs(iVal(4)-hermF)), ' & ', to_lx(abs(iVal(4)-hermF)/iVal(4)) 
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