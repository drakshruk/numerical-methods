program main
    integer :: studNum
    real(8) :: a, mat1(1:5,1:5), mat2(1:5,1:5), det, kNorm, lNorm, mNorm
    studNum = 1
    a = 0.1 - 0.001*studNum
    mat1(1,1:5) = [a*a-1d0, 0.1*a+0.1, 0.2*a+0.2, 0.3*a+0.3, 0.4*a+0.4]
    mat1(2,1:5) = [0.1*a-0.1, a*a-3.99, 0.5*a+1.02, 0.6*a+1.23, 0.7*a+1.44]
    mat1(3,1:5) = [0.2*a-0.2, 0.5*a-0.98, a*a+0.29, 0.8*a+0.36, 0.9*a+0.43]
    mat1(4,1:5) = [0.3*a-0.3, 0.6*a-1.17, 0.8*a+0.36, a*a - 14.91, a-2.74]
    mat1(5,1:5) = [0.4*a-0.4, 0.7*a-1.36, 0.9*a+0.43, a+5.26, a*a-6.54]
    write(*,'(*(f12.6))') transpose(mat1(:,:))
    call gaussDet(mat1,5, det)
    write(*,*) det

    a = 1d0 - 0.01*studNum
    mat2(1,1:5) = [5d0, -4d0/5/a, 1d0/5/a/a, -1d0/15/a/a/a, 1d0/30/a/a/a/a]
    mat2(2,1:5) = [-20d0*a, 4d0, -1d0/a, 1d0/3/a/a, -1d0/6/a/a/a]
    mat2(3,1:5) = [45d0*a*a, -9d0*a, 3d0, -1d0/a, 1d0/2/a/a]
    mat2(4,1:5) = [-60d0*a*a*a, 12d0*a*a, -4d0*a, 2d0, -1d0/a]
    mat2(5,1:5) = [30d0*a*a*a*a, -6d0*a*a*a, 2d0*a*a, -a, 1d0]
    call gaussRev(mat2,mat1, 5)

    write(*,'(*(f12.6))') transpose(mat2(:,:))
    write(*,'(*(f12.6))') transpose(mat1(:,:))

    write(*,'(*(f12.6))') mNorm(mat2,5), lNorm(mat2,5), kNorm(mat2,5)
    write(*,'(*(f12.6))') mNorm(mat1,5), lNorm(mat1,5), kNorm(mat1,5)
    pause
end program