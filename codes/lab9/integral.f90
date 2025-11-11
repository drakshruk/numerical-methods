
real(8) function cRect(hx,hy,nx,ny,F)
    implicit none
    integer :: nx, ny, i, j
    real(8) :: F(0:nx,0:ny), hx, hy
    cRect = 0d0
    do i = 0, nx
        do j = 0, ny
            cRect = cRect + F(i,j)
        end do
    end do
    cRect = cRect * hx * hy
end function cRect

real(8) function trapez(hx,hy,nx,ny,F)
    implicit none
    integer :: nx, ny, i, j
    real(8) :: F(0:nx,0:ny), hx, hy
    trapez = 0d0
    do i = 0, nx-1
        do j = 0, ny-1
            trapez = trapez + F(i,j)/2 + F(i+1,j+1)/2
        end do
    end do
    trapez = trapez * hx * hy
end function trapez

real(8) function simpson(hx,hy,nx,ny,F)
    implicit none
    integer :: nx, ny, i, j
    real(8) :: F(0:nx,0:ny), hx, hy
    simpson = 0d0
    do i = 0, nx/2-1
        do j = 0, ny/2-1
            simpson = simpson + F(2*i,2*j)  + 4*F(2*i+1,2*j)    + F(2*i+2,2*j)
            simpson = simpson + F(2*i,2*j+1) + 4*F(2*i+1,2*j+1) + F(2*i+2,2*j+1)
            simpson = simpson + F(2*i,2*j+2) + 4*F(2*i+1,2*j+2) + F(2*i+2,2*j+2)
        end do
    end do
    simpson = simpson * hx * hy / 9*2
end function simpson

real(8) function f1(x,y)
    implicit none
    real(8) x,y
    f1 = x*x+2*y
end function f1

real(8) function f2(x,y)
    implicit none
    real(8) x,y
    f2 = 1/(x+y)**2
end function f2
    
real(8) function f3(x,y)
    implicit none
    real(8) x,y
    f3 = x*x/(1+y*y)
end function f3