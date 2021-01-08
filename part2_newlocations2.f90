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
                !write(*,*) 
                !write(*,*) "SOR Subroutine Completed:"
                !write(*,*)
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
    integer :: i, j
    real(wp), dimension(:,:), allocatable :: u,f, uinitial, diff, gravelpit
    !------------------------------------------------------------
    open(unit=20,file=output)
    open(unit=30,file=output1)
    allocate(u(nx,ny))
    allocate(f(nx,ny))
    allocate(diff(nx,ny))
    allocate(gravelpit(nx,ny))
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
    Q = 1!initial flow guess
    !headloss = 4.0 ! minimum decrease in head within the gravel pit
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
        f(5,19) = (Q)/(dx*dy) !these two are above gravel pit
        f(6,19) = (Q)/(dx*dy) 
        f(8,16) = (Q)/(dx*dy) !right side of gravel pit
        f(3,16) = (Q)/(dx*dy) !left side of gravel pit
        call sor(u,f) !gets the new u(head) matrix with the wells at each location
        !a new u matrix was just calculated with the wells
        !stores diff of new head matrix and initial head matrix, this diff should equal 4
        diff = uinitial - u ! the new value should be 4 less than the original
        !this loop checks each node within the gravel pit to make make it at least 4 less
        do i = 4, 7 !range of the gravel pit square
            do j = 14,18
                if (0 < (diff(i,j) - 4.0) .and. (diff(i,j) - 4.0) < 0.001) then !if the diff array is 4 or more at each node, success
                    write(*,*) "done iterating"
                    write(*,*) "final Q=",Q
                    go to 400
                !else if (diff(i,j) > 4) then ! dont need this because diff must be greater than 4
                    !Q = Q-1
                else if (diff(i,j) <= 4.0) then 
                    Q = Q+ 0.005 !increase the flow to make sure each node is at least 4 less
                    !write(*,*) "Q= ",Q
                    ! a higher flow rate will mean a lower head because it is pumping out
                end if
            end do
        end do
    end do
    400 continue
    !write out the final u matrix that is 4 lower in the gravel pit
    gravelpit = u(4:7,14:18)
    do j = ny,1,-1
        write(*,afmt) (u(i,j), i = 1, nx)
        write(30,afmt) (u(i,j), i = 1, nx)
    end do
    write(*,*) "Difference:"
    do j = ny,1,-1
        write(*,afmt) (diff(i,j), i = 1, nx)
        !write(30,afmt) (diff(i,j), i = 1, nx)
    end do

STOP
end program well_analysis