!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! The subroutines in this module take input from ths user for the various 
!! parameters needed to calculate the flux of the box and the flux of the box 
!! with the sphere inside. It writes the various results to results_basic.dat 
!! and results_advanced.dat
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! read_input
!! read_advanced_input
!!----------------------------------------------------------------------
!! Included functions:
!!
!! read_real
!! read_integer
!-----------------------------------------------------------------------
module read_write
use types
use neutron_flux

implicit none

private
public :: read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine reads input from the user and sends it to other 
!! subroutines to make sure it is valid for the problem and sends it to 
!! other subroutines to write the resultant flux into a file. 
!!----------------------------------------------------------------------
!! Output:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(out) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(out) ::  n_grid, m_samples

    
    print *, 'This program calculates neutron flux of a reactor' 
    print *, 'The reactor is assumed to be a rectangular shape with' 
    print *, 'Depth D, width W, height H with a detector placed y0 from' 
    print *, 'One of the corners and x0 away from the front.' 
    print *, 'This program will run a loop from a minimum detector position' 
    print *, 'To a maximum detector position in an increment of x_step.'

   
    ! In order to ask for a different input in the screen, the new function
    ! we'll define below will take a string as input with the name of the
    ! variable we want the user to give us.
    depth = read_real('depth D')
    width = read_real('width W')
    height = read_real('height h')
    y_zero = read_real("y coordinate of the detector's position y_zero") 
    x_min = read_real("minimum x coordinate of the detector's position x_min") 
    x_max = read_real("maximum x coordinate of the detector's position x_max") 
    x_step = read_real("increment size for the x coordinate of the detector's position x_step") 
    
    
    

    ! The function read_real used above returns a double precision real,
    ! however n_grin and m_samples are integers, that means that we need
    ! another function to get integers. For that we'll define read_integer
    ! as well
    n_grid = read_integer('number of lattice points N')
    m_samples = read_integer('number Monte Carlo samples M')

    print *,
end subroutine read_input

!-----------------------------------------------------------------------
!! Function: read_real
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function takes in input of a parameter from the previous subroutine 
!! and makes sure it is a real number greater than 0. 
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        real        A positive non negative number given by the user
!-----------------------------------------------------------------------
real(dp) function read_real(name) result(x)
    implicit none 
    
    character(len=*), intent(in) :: name
    character(len=120) :: string
    integer :: ierror

    print *, 
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Here we are using a do loop to check if the user entered a real number greater than 0. 
    ! If the user entered something non empty and if it is a real number greater than 0 then the loop 
    ! exits because the input is valid. 
    do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .AND. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a number or is a negative number, please provide a non negative number'
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
    
end function read_real

!-----------------------------------------------------------------------
!! Function: read_integer
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function takes in an integer from the user and makes sure the user 
!! has entered an integer greater than 0.
!!----------------------------------------------------------------------
!! Input:
!!
!! name     character   A string with a brief description of the value being asked for
!!----------------------------------------------------------------------
!! Output:
!!
!! x        integer     A positive non negative number given by the user
!-----------------------------------------------------------------------
integer function read_integer(name) result(x)
    implicit none
    character(len=*), intent(in) :: name
    character(len = 120) :: string
    integer :: ierror

    print *,
    print *, 'Provide a nonzero positive value for the '//trim(name)//':'

    ! Here we are using a do loop to check if the user entered an integer greater than 0. 
    ! If the user entered something non empty and if it is an integer greater than 0 then the loop 
    ! exits because the input is valid.

 

    do
        read(*,'(a)', iostat=ierror) string
        if(string /= '') then
            read(string, *, iostat=ierror) x
            if (ierror == 0 .and. x > 0 ) exit
            print *, "'"//trim(string)//"'"//' is not a number or is a negative number, please provide a non negative number'
        else
            print *, 'that was an empty input, please provide a number'
        endif
    enddo
end function read_integer

!-----------------------------------------------------------------------
!Subroutine: write_neutron_flux
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine takes in the parameters of the box neutron source 
!! and the position of the detector and parameters of monte carlo and booles 
!! integration and writes the resultant flux due to the rectangular neutron 
!! source into a file called results_basic.dat
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! y_zero       real        y coordinate of the detector's position
!! x_min        real        minimum x coordinate of the detector's position
!! x_max        real        maximum x coordinate of the detector's position
!! x_step       real        increment size for the x coordinate of the detector's position
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
    implicit none
    real(dp), intent(in) :: depth, width, height, y_zero, x_min, x_max, x_step
    integer, intent(in) :: n_grid, m_samples

    real(dp) :: x_zero, box_booles, box_mc, box_large_x0, sigma_box
    character(len=*), parameter :: file_name = 'results_basic.dat'
    integer :: unit

    open(newunit=unit,file=file_name)
    write(unit,'(5a28)') 'x_0', 'booles', 'large x_0', 'monte carlo', 'MC uncertainty'
    x_zero = x_min
    do 
        if(x_zero > x_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        box_large_x0 = large_x0_flux(depth, width, height, x_zero, y_zero)
        call box_flux_monte_carlo(depth, width, height, x_zero, y_zero, m_samples, box_mc, sigma_box)
        write(unit,'(5e28.16)') x_zero, box_booles, box_large_x0, box_mc, sigma_box
        x_zero = x_zero + x_step
    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    print *, 'This next part of the program will calculate the flux of a neutron detector' 
    print *, 'with a hollow sphere inside'
end subroutine write_neutron_flux



!-----------------------------------------------------------------------
!! Subroutine: read_advanced_input
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine acts similar to the other subroutine that reads input. 
!! It takes in input from the user for the minimum and max radius of the hollow sphere inside 
!! the neutron source and the increment from the minimum to maximum radius and passes 
!! it to the read_real function that makes sure the input is valid.
!!----------------------------------------------------------------------
!! Output:
!!
!! x_zero       real        x coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!-----------------------------------------------------------------------
  subroutine read_advanced_input(x_zero, r_min, r_max, r_step)
      implicit none
      real(dp), intent(out) :: x_zero, r_min, r_max, r_step
    x_zero = read_real("x coordinate of the detector's position")
    r_min = read_real("minimum radius of the hollow sphere r_min") 
    r_max = read_real("maximum radius of the hollow sphere r_max") 
    r_step = read_real("increment size for the radius of the hollow sphere r_step") 
      ! Base this new subroutine on the read_input one
  end subroutine read_advanced_input

!-----------------------------------------------------------------------
!Subroutine: write_advanced_flux
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine uses the user input of the parameters of the sphere 
!! as well as earlier defined parameters of the neutron detector and box 
!! source to pass it to other subroutines and functions to calculate total flux 
!! and write the resultant flux into a file called results_advanced.dat
!!----------------------------------------------------------------------
!! Input:
!! 
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! r_min        real        minimum radius of the hollow sphere
!! r_max        real        maximum radius of the hollow sphere
!! r_step       real        increment size for the radius of the hollow sphere
!! n_grid       integer     number of grid points in each dimension of Boole's integration
!! m_samples    integer     number of sample points in the Monte Carlo integration
!-----------------------------------------------------------------------
subroutine write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)
      implicit none
      real(dp), intent(in) :: depth, width, height, x_zero, y_zero, r_min, r_max, r_step
      integer, intent(in) :: n_grid, m_samples

      real(dp) :: radius, box_booles, hollow_booles, hollow_mc, sigma_hollow, distance
      character(len=*), parameter :: file_name = 'results_advanced.dat'
      integer :: unit

      open(newunit=unit, file=file_name)
      write(unit,'(5a28)') 'radius', 'box booles', 'hollow booles', 'hollow monte carlo', 'MC uncertainty'
      radius = r_min
     
     ! The goal is to compare Boole's and Monte Carlo integration when there's a hollow 
     ! sphere inside the reactor to the calculation of the solid reactor using Boole's method
      do 
        if(radius > r_max) exit
        box_booles = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid)
        hollow_booles = total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid) 
        call hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, m_samples, hollow_mc, sigma_hollow) 
        write(unit,'(5e28.16)') radius, box_booles, hollow_booles, hollow_mc, sigma_hollow 
        radius = radius + r_step 

    enddo
    close(unit)
    print *, 'The fluxes were written in the '//file_name//' file'
    print *,

end subroutine write_advanced_flux

end module read_write
