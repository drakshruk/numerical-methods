program main
    real(8) :: a, hx, ht, x0, x1, t0, t1
    integer :: studNum, N, M
    real(8), allocatable :: u(:,:)
    
    studNum = 1
    x0 = 0d0; x1 = 1d0
    t0 = 0d0; t1 = 1d0
    ht = 2d-1; hx = 1d-1
    N = int((t1-t0)/ht + 5d-1); M = int((x1-x0)/hx + 5d-1)
    a = -5d-1 + 0.001d0*studNum

    allocate(u(0:N,0:M))
    u(0,0:M) = [(nUsl(x0+hx*i), i = 0,M)]
    u(0:N,0) = [(gUsl(t0+ht*i), i = 0,N)]

    do i = 0, N-1
        do j = 0, M-1
            u(i+1,j+1) = (a*(i+1)*ht - (u(i+1,j) - u(i,j))/ht)*(2d0*a+1d0)/a*hx + u(i+1,j)
        end do
    end do
    
    do i = 0, M
        write(*,'(*(f12.6))') u(0:N,i)
    end do

    write(*,*)
    
    do i = 0, M
        write(*,'(*(f12.6))') [(uReal(x0+hx*i,t0+ht*j), j = 0,N)]
    end do

    write(*,*)
    
    do i = 0, M
        write(*,*) [(abs(u(j,i) - uReal(x0+hx*i,t0+ht*j)), j = 0,N)]
    end do

    write(*,*)
    
    do i = 0, M
        write(*,*) [(abs(u(j,i) - uReal(x0+hx*i,t0+ht*j))/uReal(x0+hx*i,t0+ht*j), j = 0,N)]
    end do

    pause
    contains
    real(8) function nUsl(x)
        real(8) :: x
        nUsl = 3d0+(2d0*a+1d0)*x
    end function
    
    real(8) function gUsl(t)
        real(8) :: t
        gUsl = 3d0-a*t+a*t*t/2d0
    end function

    real(8) function uReal(x,t)
        real(8) :: x, t
        uReal = 3d0 + (2d0*a + 1d0)*x + a*t*t/2d0 - a*t
    end function
end program