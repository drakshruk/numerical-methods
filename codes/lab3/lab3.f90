program lab3
    integer n
    real(8) :: pogr1(1:39), pogr2(1:39)
    interface fourier_interface
        subroutine fourier(m)
            integer m
        end subroutine
    end interface fourier_interface
    interface calcLezh_interface
        subroutine calcLezhCoefficients(m, pogr)
            integer m
            real(8) :: pogr(1:5)
        end subroutine
    end interface calcLezh_interface
    interface calcCheb_interface
        subroutine calcChebCoefficients(m, pogr)
            integer m
            real(8) :: pogr(1:5)
        end subroutine
    end interface calcCheb_interface
    n = 8
    call fourier(n)
    ! write(*,*)'lezh'
    ! call calcLezhCoefficients(n, pogr1)
    ! write(*,*)'cheb'
    ! call calcChebCoefficients(n, pogr1)
    write(*,*)
    write(*,*)
    n = 9
    call fourier(n)
    ! write(*,*)'lezh'
    ! call calcLezhCoefficients(n, pogr2)
    ! write(*,*)'cheb'
    ! call calcChebCoefficients(n, pogr2)
    ! write(*,'(f12.8)') abs(pogr2-pogr1)
    pause
    end program