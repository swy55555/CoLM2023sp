#include <define.h>
module calalb
    use, intrinsic :: iso_fortran_env, only: real64, real32
    use pseudoinverse_module
    !use :: pinverse, only: pinv
    IMPLICIT NONE
    contains
    subroutine readVlam(filename, Vlam, Vlam_inv)
        character(len=*), intent(in) :: filename
        real(kind=real64), dimension(:,:), intent(out) :: Vlam
        real(kind=real64), dimension(:,:), intent(out) :: Vlam_inv
        real(kind=real64),  dimension(:), allocatable :: vec
        integer :: i, j, k, m, n, stat
        k = 1
        m = 4
        n = 190

        open(unit=10, file=filename, form='unformatted', access='stream', status='old', iostat=stat)
        if (stat /= 0) stop 'Cannot open the file'

        allocate(vec(n*m), stat=stat)
        if (stat /= 0) stop 'Cannot allocate memory'

        read(10) vec
        close(10)

        do i = 1, m
            do j =1, n
                Vlam(j,i) = vec(k)
                k = k + 1
            end do
        end do

        !do i = 1, n
        !    write(*, "(5F8.2)") (Vlam(i, j), j = 1, m)
        !end do
        call pseudoinverse(Vlam,Vlam_inv)
    end subroutine readVlam
    SUBROUTINE landalb(alb,filename,Vlam,Vlam_inv,albsp)
        CHARACTER(len=*), intent(in):: filename
        real(kind=real64), dimension(:,:), intent(inout) :: Vlam
        real(kind=real64), dimension(:,:), intent(inout) :: Vlam_inv
        real(kind=real64), dimension(:,:), intent(in) :: alb
        real(kind=real64), dimension(190,2), intent(inout) :: albsp
        integer :: stat, i
        call readVlam(filename, Vlam, Vlam_inv)
        
        do i = 1, 30
            albsp(i,1) = alb(1,1)
        end do
        do i = 31, 190
            albsp(i,1) = alb(2,1)
        end do
        do i = 1, 30
            albsp(i,2) = alb(1,2)
        end do
        do i = 31, 190
            albsp(i,2) = alb(2,2)
        end do
        albsp(:,1) = matmul(matmul(albsp(:,1), transpose(Vlam_inv)), transpose(Vlam))
        albsp(:,2) = matmul(matmul(albsp(:,2), transpose(Vlam_inv)), transpose(Vlam))
        !print *,albsp(1,:)
    end SUBROUTINE
end module
        
