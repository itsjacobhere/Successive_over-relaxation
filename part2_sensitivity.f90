module vars
    implicit none
    integer, parameter :: wp = selected_real_kind(15)
    !subroutine SOR variables:
    real(wp) :: w, eps1,eps2, p
    !main program variables:
    real(wp) :: kx = 10, ky = 10, dx= 100, dy = 100, Q
    real(wp) :: transmissivity, satthickness, px, py,tx,ty
    integer :: ixstart , iystart 
    integer, parameter :: nx = 20, ny = 21, ixend= nx, iyend = ny-1
    !other vars:
    character(len=20) :: output="case_analysis.csv", afmt = '(20(f10.4,","))', &
    output1 = "well_analysis.csv"
end module
!____________________________________________________________________________________
module Succesive_over_relaxation
    use vars
    implicit none
    contains
    !SOR algorithms for 2d elliptical PDE
    subroutine SOR(u,f)
        implicit none
        real(wp), intent (inout), dimension(:,:), allocatable :: u,f
        integer :: j,i, iterations, maxiterations
        real(wp) :: maxupdate
        real(wp):: update
        !------------------------------------------------------------------
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
                        - 2*(px+py)*u(i,j) - f(i,j))/(2*(px+py))
                    else
                        update = w*(px*u(i-1,j) + px*u(i+1,j) + py*u(i,j-1) &
                        + py*u(i,j+1) - 2*(px+py)*u(i,j) - f(i,j))/(2*(px+py))
                    end if
                    if (maxupdate < abs(update)) then
                        maxupdate=abs(update)
                    end if
                    u(i,j) = u(i,j) + update
                end do
            end do
            iterations = iterations+1
            if (maxupdate<eps1) then
                !write(*,*) "SOR Subroutine Completed:"
                exit 
            else if (maxiterations<iterations) then
                write(*,*) "Too many iterations!"
                STOP
            end if
        end do
    end subroutine SOR
end module Succesive_over_relaxation
!____________________________________________________________________________________
program well_analysis
    use vars
    use Succesive_over_relaxation
    implicit none
    integer :: i, j, c
    real(wp), dimension(:,:), allocatable :: u,f, uinitial, diff, pit
    !------------------------------------------------------------
    open(unit=20,file=output)
    open(unit=30,file=output1)
    open(unit=40,file="Q20.csv")
    open(unit=50,file="diff_wells.csv")
    allocate(u(nx,ny))
    allocate(f(nx,ny))
    allocate(diff(nx,ny))
    write(*,*) "----------------------------------------------------------------------------------------------------------"
    !initialize start and end
    f = 0
    ixstart = 2
    iystart = 2
    u=55
    w=1 !relaxation factor
    satthickness = 60
    transmissivity = kx*satthickness
    tx = transmissivity
    ty = transmissivity
    
    px = kx/(dx**2)
    py = ky/(dy**2)
    
    !-----------------------------------------------------------
    !Part One:
    if ((abs(px-py)) < 0.0001) p = px
    do j = 1, ny
        u(1,j) = 50 !initialize the first column of the matrix
    end do
    
    do i = ixstart,nx
        u(i,1) = 0.0025*(i-1)*dx + 50 !initialize top
        u(i,ny) = 0.0005*(i-1)*dx + 50 !initialize bottom
    end do
    
    call SOR(u,f) 
    
    do j = ny,1,-1
        write(*,afmt) (u(i,j), i = 1, nx)
        write(20,afmt) (u(i,j), i = 1, nx)
    end do
    !-----------------------------------------------------------
    !Part Two:
    !Solve for the new PDE in this section for Part B
    uinitial = u !store the initial matrix from part one in a new variable, will use to comapare
    write(*,*) 
    write(40,*) "transmissivity,percent, Q"
    do c = 80, 120, 1
        Q = 1!initial flow guess
        !headloss = 4.0 ! minimum decrease in head within the gravel pit
        transmissivity = kx*satthickness
        transmissivity = transmissivity*((real(c))/100)
        tx = transmissivity
        ty = transmissivity
        px = tx/(dx**2)
        py = ty/(dy**2)
        do 
            f = 0 !f is the change in head due to the pump
            !assign well locations around the gravel pit
            !gravel pit is from i = 4 to 7, and j = 14 to 18
            f(4,13) = (Q)/(dx*dy) 
            f(5,13) = (Q)/(dx*dy) !these are below gravel pit
            f(6,13) = (Q)/(dx*dy) 
            f(7,13) = (Q)/(dx*dy)
            f(5,19) = (Q)/(dx*dy) ! above gravel pit
            f(6,19) = (Q)/(dx*dy) !above gravel pit
            f(8,16) = (Q)/(dx*dy) !right side of gravel pit
            f(3,16) = (Q)/(dx*dy) !left side of gravel pit
            call sor(u,f) !gets the new u(head) matrix with the wells at each location
            !a new u matrix was just calculated with the wells
            !stores diff of new head matrix and initial head matrix, this diff should equal 4
            !diff(4:7,14:18) = uinitial(4:7,14:18) - u(4:7,14:18) ! the new value should be 4 less than the original
            diff = uinitial - u !difference between old
            !this loop checks each node within the gravel pit to make make it at least 4 less
            if ( (minval(diff(4:7,14:18)) >= 4.0)) then
                write(*,*)
                write(*,*) "Found the new U matrix with head lowered 4m in the gravel pit:"
                write(*,*)
                go to 400
            else if ( (minval(diff(4:7,14:18)) <= 4.0)) then
                Q = Q +0.005
            end if
        end do
        400 continue
        !write out the final u matrix that is 4 lower in the gravel pit
        pit = u(4:7,14:18)
        do j = ny,1,-1
            !write(*,afmt) (u(i,j), i = 1, nx)
            write(30,afmt) (u(i,j), i = 1, nx)
        end do
        write(*,*)
        write(*,*) "Difference in the old matrix and the new lowered matrix:"
        write(*,*)
        do j = ny,1,-1
            write(50,afmt) (diff(i,j), i = 1, nx)
        end do
        write(50,*) transmissivity, Q
        write(40,'(f10.4,",",I3,",",f10.4)') transmissivity,c-100,Q
        write(*,'(f10.4,",",I3,",",f10.4)') transmissivity,c-100,Q
    end do
stop
end program well_analysis