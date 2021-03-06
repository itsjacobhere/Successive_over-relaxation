! By Jacob Turner on 11/14/18
! Part 1: uses successive over relaxation algorithm to calculate a  steady-state head matrix in an aquifer
! by evaluating an elliptical Partial differential equation
! Part 2: uses the old matrix from part 1 to determine a flow that will lower the head in a the gravel pit
! by at least 4 m. the gravel pit is in nodes (4:7.14:18). the Q is searched for from an initial guess.
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
    
    !kx/ky = hydraulic conductivity
    !tx/ty = transmissivity
    !dx/dy = x and y location increment
    !w = relaxation factor
    !Q = flow needed to lower head in Head matrix
    !u = Head matrix calculated by SOR sub
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
        !py = px
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
                !completed
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
    integer :: i, j, s = 0
    real(wp), dimension(:,:), allocatable :: u,f, uinitial, diff
    !------------------------------------------------------------
    open(unit=20,file=output)
    open(unit=30,file=output1)
    allocate(u(nx,ny))
    allocate(f(nx,ny))
    allocate(diff(nx,ny))

    write(*,*) "----------------------------------------------------------------------------------------------------------"
    !initialize values
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
    !write out original u head matrix here
    write(*,*) "The original head matrix is shown below:"
    write(*,*) "=============================================================================================="
    do j = ny,1,-1
        write(*,afmt) (u(i,j), i = 1, nx)
        write(20,afmt) (u(i,j), i = 1, nx)
    end do
    !-----------------------------------------------------------
    !Part Two:
    !Solve for the new PDE in this section for Part B
    uinitial = u !store the initial matrix from part one in a new variable, will use to comapare
    Q = 1! m3/day!initial flow guess, start low and increase until conditions are satisfied.
    !headloss = 4.0 ! minimum decrease in head within the gravel pit
    transmissivity = kx*satthickness
    !transmissivity = transmissivity*((real(c))/100)
    tx = transmissivity
    ty = transmissivity
    px = tx/(dx**2)
    py = ty/(dy**2)
    do 
        s = s + 1
        f = 0 !f matrix should be zero except where the wells are located
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
        call sor(u,f) !gets the new u(head) matrix with the wells factored in at each location
        !a new u matrix was just calculated with the wells
        !stores diff of new head matrix and initial head matrix, this diff should equal 4
        !diff(4:7,14:18) = uinitial(4:7,14:18) - u(4:7,14:18) ! the new value should be 4 less than the original
        diff = uinitial - u !difference between old and the matrix just found (u)
        !this loop checks each node within the gravel pit to make make it at least 4 less
        if ( (minval(diff(4:7,14:18)) >= 4.0)) then
            write(*,*)
            write(*,*) "Found the new U matrix with head lowered 4m in the gravel pit:"
            go to 400
        else if ( (minval(diff(4:7,14:18)) < 4.0)) then
            Q = Q + 0.005 ! change Q by a small amount until the min val in pit is decreased by at least 4 m
        end if
    end do
    400 continue
    !finished Just write out the results to output and to screen
    write(*,*) "=============================================================================================="
    !write out the final u matrix that is 4 lower in the gravel pit
    do j = ny,1,-1
        write(*,afmt) (u(i,j), i = 1, nx)
        write(30,afmt) (u(i,j), i = 1, nx)
    end do
    write(*,*)
    write(*,*) "Difference in the old matrix and the new lowered matrix:"
    write(*,*) "==============================================================================================="
    do j = ny,1,-1
        write(*,afmt) (diff(i,j), i = 1, nx)
    end do
    write(*,*) "Final flow needed to lower the head in the gravel pit by 4 (m) =",Q,"m3/day"
    STOP
end program well_analysis