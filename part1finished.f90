module vars
    implicit none
    integer, parameter :: wp = selected_real_kind(15)
    !subroutine SOR variables:
    real(wp) :: w, eps1,eps2, p
    !main program variables:
    real(wp), parameter :: kx = 10, ky = 10, dx= 100, dy = 100
    real(wp) :: transmissivity, satthickness, px, py
    integer :: ixstart , iystart 
    integer, parameter :: nx = 20, ny = 21, ixend= nx, iyend = ny-1
    !other vars:
    character(len=20) :: output="case_analysis.csv", afmt = '(20(f10.4,","))'
end module
!____________________________________________________________________________________
module Succesive_over_relaxation
    use vars
    implicit none
    contains
    !SOR algorithms for 2d elliptical PDE
    subroutine SOR(u)
        implicit none
!        real(wp), intent(in) :: w, eps1,eps2
        !integer, intent(in) :: ixstart, ixend, iystart, iyend, nx, ny
        !real(wp), dimension(:,:), allocatable, intent(inout) :: u
        !local vars:
        real(wp), intent (inout), dimension(:,:), allocatable :: u
        integer :: j, k, i, iterations, maxiterations
        real(wp) :: maxupdate
        real(wp):: update
        !------------------------------------------------------------------
        w = 1
        eps1 = 0.01
        iterations = 0
        maxiterations = 500
        iterations = 0
        py = px
        do 
            maxupdate = 0 
            do i = ixstart,nx
                do j = iystart,ny-1
                    if (i == nx) then
                        update = w*(2*px*u(i-1,j) + py*u(i,j-1) + py*u(i,j+1) &
                        - 2*(px+py)*u(i,j))/(2*(px+py))
                    else
                        update = w*(px*u(i-1,j) + px*u(i+1,j) + py*u(i,j-1) &
                        + py*u(i,j+1) - 2*(px+py)*u(i,j))/(2*(px+py))
                    end if
                    if (maxupdate < abs(update)) then
                        maxupdate=abs(update)
                    end if
                    u(i,j) = u(i,j) + update
                end do
            end do
            iterations = iterations+1
            if (maxupdate<eps1) then
                write(*,*) 
                write(*,*) "SOR Subroutine Completed:"
                write(*,*)
                exit 
            else if (maxiterations<iterations) then
                write(*,*) "Too many iterations!"
                STOP
            end if
        end do
        !write it out here
!        do j = ny,1,-1
!            write(*,afmt) (u(i,j), i = 1, nx)
!            write(20,afmt) (u(i,j), i = 1, nx)
!        end do
    end subroutine SOR
end module Succesive_over_relaxation
!____________________________________________________________________________________
program well_analysis
    use vars
    use Succesive_over_relaxation
    implicit none
    integer :: i, j
    real(wp), dimension(:,:), allocatable :: u
    !------------------------------------------------------------
    open(unit=20,file=output)
    allocate(u(nx,ny))
    !initialize start and end
    ixstart = 2
    iystart = 2
    u=55
    w=1
    px = kx/(dx**2)
    py = ky/(dy**2)

    if ((abs(px-py)) < 0.0001) p = px
    do j = 1, ny
        u(1,j) = 50 !initialize the first column of the matrix
    end do
    do i = ixstart,nx
        u(i,1) = 0.0025*(i-1)*dx + 50 !initialize top
        u(i,ny) = 0.0005*(i-1)*dx + 50 !initialize bottom
    end do
    call SOR(u)
    do j = ny,1,-1
            write(*,afmt) (u(i,j), i = 1, nx)
            write(20,afmt) (u(i,j), i = 1, nx)
    end do
STOP
end program well_analysis