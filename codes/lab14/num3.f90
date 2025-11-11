subroutine num3
    interface func
        function to_lx(number) result(output_string)
            implicit none
            real(kind=8), intent(in) :: number
            character(len=50) :: output_string
        end function to_lx
    end interface func
    integer :: studNum, i, n
    real(8) :: w1, w2, w3, a, tmp, xtmp
    real(8) :: x, l, xs, lx, ls, x0, x1, s
    real(8), allocatable :: y1_0(:), y1_1(:), y1_2(:)
    real(8), allocatable :: y2_0(:), y2_1(:), y2_2(:)
    real(8), allocatable :: y3_0(:), y3_1(:), y3_2(:)
    character(50) :: format

    
    studNum = 1
    s = 0.64d0
    a = -1d0 + 0.001d0*studNum
    x0 = 0.64d0;    x1 = 0.81d0
    x0 = 0.64d0;    x1 = 0.81d0
    n = (x1-x0)/0.01d0
    allocate(y1_0(0:n), y1_1(0:n), y1_2(0:n))
    allocate(y2_0(0:n), y3_1(0:n), y2_2(0:n))
    allocate(y3_0(0:n), y2_1(0:n), y3_2(0:n))
    w1 = 0.64d0*a+0.8d0
    w2 = 0.64d0*a+1.6d0
    w3= 0.64d0*a-0.8d0

	do i = 0, n
        x = 0.64d0 + 0.01d0*i
        l = log(x / s)
        xs = sqrt(x) - sqrt(s)
        lx = log(x); ls = log(s)
		y1_0(i) = w1 + l * w2 - 3d0 * xs
		y2_0(i) = w2 + w3 * l + 4d0 * xs
		y3_0(i) = w3 + w1 * l - 3d0 * xs
		y1_1(i) = (w3 * l * l + 2d0 * (w2 * l + w1) - 8d0 * sqrt(s) * l + 10d0 * xs) / 2d0
		y2_1(i) = (w1 * l * l + (2d0 * w3 + 6d0 * sqrt(s)) * l + 2d0 * w2 - 4d0 * xs) / 2d0
		y3_1(i) = (w2 * l * l + (2d0 * w1 + 6d0 * sqrt(s)) * l + 2d0 * w3 - 18d0 * xs) / 2d0

		y1_2(i) = w1 * lx * lx * lx - (3d0 * ls * w1 - 3d0 * w3) * lx * lx
		y1_2(i) = y1_2(i) - sqrt(s) * (-9d0 * lx * lx + (18d0 * ls - 12d0) * lx - 9d0 * ls * ls + 12d0 * ls - 42d0)
		y1_2(i) = y1_2(i) - (6d0 * ls * w3 - 6d0 * w2 - 3d0 * w1 * ls * ls) * lx - 42d0 * sqrt(x)
		y1_2(i) = (y1_2(i) + 3d0 * w3 * ls * ls - 6d0 * ls * w2 - (ls * ls * ls - 6d0) * w1)/6d0

		y2_2(i) = -6d0 * w2 * lx * lx * lx + 18d0 * w2 * ls * lx * lx
		y2_2(i) = y2_2(i) + 504d0 * sqrt(x) - 18d0 * w1 * l * l - 54d0 * sqrt(s) * l * l - 36d0 * w3 * l
		y2_2(i) = y2_2(i) - 18d0 * w2 * ls * ls * l - 324d0 * sqrt(s) * l - 12d0 * ls * ls * ls * w2 - 504d0 * sqrt(s)
        y2_2(i) = w2 - (y2_2(i))/36d0

        y3_2(i) = 32.0 / 3.0 * w3 * lx * lx * lx - 32.0 * w3 * ls * lx * lx + 448.0 * sqrt(x)
		y3_2(i) = y3_2(i) + 32.0 * w2 * l * l - 128.0 * sqrt(s) * l * l + 32.0 * w3 * ls * ls * l + 64.0 * w1 * l
		y3_2(i) = y3_2(i) - 320.0 * sqrt(s) * l + 64.0 / 3.0 * w3 * ls * ls * ls - 448.0 * sqrt(s)
        y3_2(i) = w3 + (y3_2(i))/64d0
        end do

    format = '(f12.2, a3, f12.6, a3, f12.6, a3, a30, a3, a30)'
    write(*,*)    'y1_1'
    do i = 0, n
        tmp = abs(y1_r(0.64d0+i*1d-2)-y1_0(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp,'&',y1_r(xtmp),'&', y1_0(i),'&', to_lx(tmp),'&', to_lx(tmp/y1_r(xtmp))
    end do
    write(*,*)    'y1_2'
    do i = 0, n
        tmp = abs(y1_r(0.64d0+i*1d-2)-y1_1(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp, '&',y1_r(xtmp),'&', y1_1(i),'&', to_lx(tmp),'&', to_lx(tmp/y1_r(xtmp))
    end do
    write(*,*)    'y1_3'
    do i = 0, n
        tmp = abs(y1_r(0.64d0+i*1d-2)-y1_2(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp, '&',y1_r(xtmp),'&', y1_2(i),'&', to_lx(tmp),'&', to_lx(tmp/y1_r(xtmp))
    end do
    write(*,*)


    write(*,*)    'y2_1'
    do i = 0, n
        tmp = abs(y2_r(0.64d0+i*1d-2)-y2_0(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp,'&',y2_r(xtmp),'&', y2_0(i),'&', to_lx(tmp),'&', to_lx(tmp/y2_r(xtmp))
    end do
    write(*,*)    'y2_2'
    do i = 0, n
        tmp = abs(y2_r(0.64d0+i*1d-2)-y2_1(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp, '&',y2_r(xtmp),'&', y2_1(i),'&', to_lx(tmp),'&', to_lx(tmp/y2_r(xtmp))
    end do
    write(*,*)    'y2_3'
    do i = 0, n
        tmp = abs(y2_r(0.64d0+i*1d-2)-y2_2(i))
        xtmp = 0.64d0+i*0.01d0
        ! write(*,format) xtmp, '&',y2_r(xtmp),'&', y2_2(i),'&', to_lx(tmp),'&', to_lx(tmp/y2_r(xtmp))
    end do
    write(*,*)


    write(*,*)    'y3_1'
    do i = 0, n
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_0(i))
        xtmp = 0.64d0+i*0.01d0
        write(*,format) xtmp,' & ',y3_r(xtmp),' & ', y3_0(i),' & ', to_lx(tmp),' & ', to_lx(tmp/abs(y3_r(xtmp)))
    end do
    write(*,*)    'y3_2'
    do i = 0, n
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_1(i))
        xtmp = 0.64d0+i*0.01d0
        write(*,format) xtmp, ' & ',y3_r(xtmp),' & ', y3_1(i),' & ', to_lx(tmp),' & ', to_lx(tmp/abs(y3_r(xtmp)))
    end do
    write(*,*)    'y3_3'
    do i = 0, n
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_2(i))
        xtmp = 0.64d0+i*0.01d0
        write(*,format) xtmp, ' & ',y3_r(xtmp),' & ', y3_2(i),' & ', to_lx(tmp),' & ', to_lx(tmp/abs(y3_r(xtmp)))
    end do

    do i = 0, n
        xtmp = 0.64d0+i*0.01d0
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_0(i))
        write(*,*) xtmp, tmp/abs(y3_r(xtmp))
    end do
    write(*,*)
    do i = 0, n
        xtmp = 0.64d0+i*0.01d0
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_1(i))
        write(*,*) xtmp, tmp/abs(y3_r(xtmp))
    end do
    write(*,*)
    do i = 0, n
        xtmp = 0.64d0+i*0.01d0
        tmp = abs(y3_r(0.64d0+i*1d-2)-y3_2(i))
        write(*,*) xtmp, tmp/abs(y3_r(xtmp))
    end do
    pause
    contains
    real(8) function y1_r(rx)
        real(8) :: rx
        y1_r = sqrt(rx) + a*rx
    end function

    real(8) function y2_r(rx)
        real(8) :: rx
        y2_r = 2d0*sqrt(rx) + a*rx
    end function

    real(8) function y3_r(rx)
        real(8) :: rx
        y3_r = -sqrt(rx) + a*rx
    end function
end subroutine