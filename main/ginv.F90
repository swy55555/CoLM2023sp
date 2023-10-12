#include <define.h>
module pseudoinverse_module
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: pseudoinverse

contains

    function diag(v)
        real(dp), intent(in) :: v(:)
        real(dp), allocatable :: diag(:,:)
        integer :: n, i

        n = size(v)
        allocate(diag(n,n))
        diag = 0.0_dp
        do i = 1, n
            diag(i,i) = v(i)
        end do
    end function diag

    subroutine pseudoinverse(Vlam, Vlam_inv)
        ! 使用LAPACK的奇异值分解计算矩阵的伪逆

        real(dp), intent(in) :: Vlam(:,:)  ! 输入矩阵
        real(dp), intent(out) :: Vlam_inv(:,:)  ! 输出矩阵的伪逆

        integer :: m, n, lda, ldu, ldvt, info, lwork, i, j
        real(dp), allocatable :: s(:), u(:,:), vt(:,:), work(:)
        real(dp), dimension(:,:), allocatable :: s_pinv(:,:)
        real(dp), parameter :: rcond = -1.0_dp
        integer, allocatable :: iwork(:)

        m = size(Vlam, 1)
        n = size(Vlam, 2)
        lda = m
        ldu = m
        ldvt = n

        allocate(s(min(m, n)))
        allocate(u(ldu, m))
        allocate(vt(ldvt, n))
        allocate(iwork(8*min(m, n)))
        allocate(s_pinv(n,m))
        lwork = 3*min(m,n)*min(m,n) + max(max(m,n), 4*min(m,n)*min(m,n)+4*min(m,n))
        allocate(work(lwork))

        ! 使用LAPACK的dgesdd函数进行奇异值分解
        call dgesdd('A', m, n, Vlam, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
        ! 计算s的伪逆
        where (s /= 0.0_dp)
            s = 1.0_dp / s
        elsewhere
            s = 0.0_dp
        end where
        ! 将s转换为一个方阵
        
        s_pinv = 0.0_dp
        do i = 1, min(m,n)
            s_pinv(i,i) = s(i)
        end do
        !print *, size(vt, 1), size(vt, 2) , size(s_pinv, 1), size(s_pinv, 2), size(u,1), size(u,2)
        !u = u(:, :n)
        ! 计算伪逆
        Vlam_inv = matmul(matmul(transpose(vt), s_pinv), transpose(u))
        !print *, '=============================================='
        !do i = 1, n
        !    write(*, "(5F8.2)") (Vlam_inv(i, j), j = 1, m)
        !end do
    end subroutine pseudoinverse

end module pseudoinverse_module