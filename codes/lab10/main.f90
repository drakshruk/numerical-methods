program main
    implicit none
    integer :: studNum, i
    real(8) :: eps, root, x1, x2

    studNum = 1
    eps = 1d-3
    x1 = 0d0
    x2 = 2.18629075d0

    ! call dichotomy(-3d-1,1d0,studNum,f1,eps, root)
    ! write(*,'(*(f12.8))') root, x1, abs(x1 - root), abs(x1 - root)/abs(x1)
    ! call dichotomy(2d0,3d0,studNum,f2,eps, root)
    ! write(*,'(*(f12.8))') root, x2, abs(x2 - root), abs(x2 - root)/abs(x2)

    ! call newton(-3d-1,1d0,studNum,f1,f1_der,eps, root)
    ! write(*,'(*(f12.8))') root, x1, abs(x1 - root), abs(x1 - root)/abs(x1)
    ! call newton(2d0,3d0,studNum,f2,f2_der,eps, root)
    ! write(*,'(*(f12.8))') root, x2, abs(x2 - root), abs(x2 - root)/abs(x2)

    ! call secant(-3d-1,1d0,studNum,f1,eps, root)
    ! write(*,'(*(f12.8))') root, x1, abs(x1 - root), abs(x1 - root)/abs(x1)
    ! call secant(2d0,3d0,studNum,f2,eps, root)
    ! write(*,'(*(f12.8))') root, x2, abs(x2 - root), abs(x2 - root)/abs(x2)
    
    ! call iterations(-3d-1,1d0,studNum,f1_iter,eps, root)
    ! write(*,'(*(f12.8))') root, x1, abs(x1 - root), abs(x1 - root)/abs(x1)
    ! call iterations(2d0,3d0,studNum,f2_iter,eps, root)
    ! write(*,'(*(f12.8))') root, x2, abs(x2 - root), abs(x2 - root)/abs(x2)
    
    do i = 0, 30
        write(*,*) 2d0+0.1d0*i, f2_iter(1d0+0.1d0*i, 1)
    end do
    pause
    contains
        real(8) function f1(x,N)
            real(8) :: x
            integer :: N
            f1 = exp(2*x) + 3*x - N
        end function f1

        real(8) function f1_der(x,N)
            real(8) :: x
            integer :: N
            f1_der = 2*exp(2*x) + 3
        end function f1_der

        real(8) function f1_iter(x,N)
            real(8) :: x
            integer :: N
            f1_iter = 1d0/3*(N-exp(2*x))
        end function f1_iter

        real(8) function f2(x,N)
            real(8) :: x, a
            integer :: N
            a = 1.01 - 0.001*N
            f2 = x**4 - (1+a)*x**3 - (4-a)*x*x + (4+4*a)*x - 4*a - sqrt(x-1)
        end function f2

        real(8) function f2_der(x,N)
            real(8) :: x, a
            integer :: N
            a = 1.01 - 0.001*N
            f2_der = 4*x**3 - 3*(1+a)*x*x - 2*(4-a)*x +(4+4*a) - 0.5/sqrt(x-1)
        end function f2_der
        
        real(8) function f2_iter(x,N)
            real(8) :: x, a
            integer :: N
            a = 1.01 - 0.001*N
            f2_iter = ((1+a)*x**3 + (4-a)*x*x - (4+4*a)*x + 4*a + sqrt(x-1))**(1d0/4)
        end function f2_iter
end program
