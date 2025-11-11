subroutine gaussDet(A, n, det)
    implicit none
    integer :: n,i,j
    real(8) :: A(1:n,1:n), M(1:n,1:n), det
    do i = 1, n - 1
        M(1:n,1:n) = 0d0
        do j = 1,n
            M(j,j) = 1d0
        end do
        M(i:n,i) = [(-A(j,i)/A(i,i), j = i,n)]
        m(i,i) = 1d0
        A(1:n,1:n) = matmul(M(1:n,1:n),A(1:n,1:n))
    end do
    det = 1d0
    do i = 1,n
        det = det*a(i,i)
    end do
end subroutine

subroutine gaussRev(A1, AREV,  n)
    implicit none
    integer :: n,i,j,k,s
    real(8), intent(in) :: A1(1:n,1:n)
    real(8) :: A(1:n,1:n), AREV(1:n,1:n), M(1:n,1:n), B(1:n,1:n), MB(1:n,1:n+1)
    B(1:n,1:n) = 0d0
    do i = 1,n
        B(i,i) = 0d0
    end do
    A(1:n,1:n) = A1(1:n,1:n)
    do i = 1, n - 1
        M(1:n,1:n) = 0d0
        do j = 1,n
            M(j,j) = 1d0
        end do
        M(i:n,i) = [(-A(j,i)/A(i,i), j = i,n)]
        m(i,i) = 1d0
        A(1:n,1:n) = matmul(M(1:n,1:n),A(1:n,1:n))
        B(1:n,1:n) = matmul(M(1:n,1:n),B(1:n,1:n))
    end do
    AREV(1:n,1:n) = 0d0
    do i = 1,n
        MB(1:n,1:n) = A(1:n,1:n)
        MB(1:n,n+1) = B(1:n,i)
        do j = 1,n
            MB(j,1:n+1) = MB(j,1:n+1)/MB(j,j)
        end do
        do j = n, 2, -1
            do k = j-1, 1, -1
                MB(k,j:n+1) = MB(k,j:n+1) - MB(k,j) * MB(j,j:n+1)
            end do
            AREV(i,j) = mb(j,n+1)
        end do
        AREV(i,1) = mb(1,n+1)
    end do
    AREV = transpose(AREV)
end subroutine

real(8) function lNorm(A,n)
    implicit none
    integer :: n,i,j
    real(8) :: A(1:n,1:n), tmp
    lNorm = sum([(abs(A(i,1)), i = 1,n)])
    do i = 1,n
        tmp = sum([(abs(A(j,i)), j = 1,n)])
        if(lNorm < tmp) then 
            lNorm = tmp
        end if
    end do
end function lNorm

real(8) function mNorm(A,n)
    implicit none
    integer :: n,i,j
    real(8) :: A(1:n,1:n), tmp
    mNorm = sum([(abs(A(1,i)), i = 1,n)])
    do i = 1,n
        tmp = sum([(abs(A(i,j)), j = 1,n)])
        if(mNorm < tmp) then 
            mNorm = tmp
        end if
    end do
end function mNorm

real(8) function kNorm(A,n)
    implicit none
    integer :: n,i,j
    real(8) :: A(1:n,1:n), tmp
    kNorm = sqrt(sum([(sum([(a(i,j)**2, i = 1,n)]), j = 1,n)]))
end function kNorm