subroutine scheme(u, v, h, n, studNum)
    integer :: n, studNum
    real(8) :: u(0:n), v(0:n), h, x0, c

    c = 1 - 8d0*h*h/studNum/studNum

    do i = 2, n
        u(i) = u(i-1)/c*16d0*h*h/studNum/studNum + u(i-2)*(1+8d0*h*h/studNum/studNum)/c
        u(i) = u(i) + 2d0*h/studNum*(studNum-2)/c*v(i-1) + 4d0*h/studNum/c*v(i-2)
        v(i) = 8d0*h*h/studNum/studNum*(studNum-2)/c*v(i-1) + (1+8d0*h*h/studNum/studNum)/c*v(i-2)
        v(i) = v(i) + 8d0*h/studNum*(studNum-1)/c*u(i-1) + 8d0*h/studNum/c*u(i-2)
    end do
end subroutine