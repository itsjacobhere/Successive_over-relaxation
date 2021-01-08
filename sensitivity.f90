module vars
    implicit none
    integer, parameter :: wp = selected_real_kind(15)
    !subroutine SOR variables:
    real(wp) :: w, eps1,eps2, p
    !main program variables:
    real(wp) :: kx, ky, dx= 100, dy = 100, k
    real(wp) :: transmissivity, satthickness, px, py
    integer :: ixstart , iystart , iterations
    integer, parameter :: nx = 20, ny = 21, ixend= nx, iyend = ny-1
    !other vars:
    character(len=25) ::  afmt = '(20(f10.4,","))',output="case_analysis.csv"
    !kx and ky: the hydraulic conductivity
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
        integer :: j, i, maxiterations
        real(wp) :: maxupdate
        real(wp):: update
        !------------------------------------------------------------------
        !w = 1
        eps1 = 0.0001
        iterations = 0
        maxiterations = 2000
        py = px
        !write(*,*) "k=",k
        !write(*,*) "w=",w
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
                !write(*,*) "SOR Subroutine Completed:"
                write(*,*)
                exit 
            else if (maxiterations<iterations) then
                write(*,*) "Too many iterations!"
                STOP
            end if
        end do
    end subroutine SOR
end module Succesive_over_relaxation
!________________________________   ____________________________________________________
program well_analysis
    use vars
    use Succesive_over_relaxation
    implicit none
    integer :: i, j, c
    real(wp), dimension(:,:), allocatable :: u
    character(len=10) :: answer, wfile = "w.csv", test = "test.txt"
    !-------------------------------------------------------------------------------
    call system('CLS')
    allocate(u(nx,ny))
    write(*,*)"------------------------------------------------------------------------------------"
    write(*,*)
    write(*,*) "Do you want to perform the sensitivity analysis on k or w?"
    write(*,*)
    read(*,*) answer
    open(unit=200,status="replace",file=wfile)
    write(200,*) "relaxation,iterations"
    !sensitivity on iterations vs relaxation factor
    !sensitivity on hydraulic conductiviy vs head at each node
    c=100
    do c = 50,150, 1
        if (answer=="no") then
            output = "well.csv" 
            go to 100
        end if
        write(output,'(A,I0,A)') "sensitivity_",c,".csv"
        open(unit=c,file=output,status="replace")
        100 continue
        !intitializations
        ixstart = 2
        iystart = 2
        u=55
        w=1.0 !relaxation facor
        k = 10 !hydraulic 
        if (answer=="k") then
            k = k*((real(c))/100)
        else if (answer=="w") then
            w = w*((real(c))/100)
        end if
        kx = k
        ky = k
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
        if (answer =="w") write(*,'(f5.2, ",",I0)') w, iterations
        if (answer =="w") write(200,'(f5.2, ",",I0)') w, iterations
        do j = ny,1,-1
            !write(*,afmt) (u(i,j), i = 1, nx)
            write(c,afmt) (u(i,j), i = 1, nx)
        end do
        !write(*,*) "writing to:", output
        if (answer=="no") go to 200
    end do
    200 continue
STOP
end program well_analysis