!-----------------------------------------------------------------------
!Module: neutron_flux
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! You can think of this as the physics module. This will calculate the neutron 
!! flux due to a solid box and the neutron flux due to a box with a hollow 
!! sphere inside. This module is unaware of the necessary mathematics or 
!! integration and is only concerned about the physics such as the consequences 
!! of considering the hollow sphere vs a solid box. This module also approximates 
!! a solution based off a large x0 which we can compare to our numerical result.
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! box_flux_monte_carlo
!! hollow_box_flux_mc
!!----------------------------------------------------------------------
!! Included functions:
!!
!! box_flux_booles
!! sphere_flux_booles
!! total_flux_booles
!! sphere_flux_kernel
!! flux_kernel
!! flux_kernel_vector
!! hollow_box_flux_kernel
!! large_x0_flux
!-----------------------------------------------------------------------
module neutron_flux
use types
use quadrature, only : booles_quadrature, monte_carlo_quad

implicit none

private
public :: box_flux_booles, large_x0_flux, box_flux_monte_carlo, hollow_box_flux_mc, total_flux_booles, sphere_flux_booles

contains

! Let's do the Boole's integration first

!-----------------------------------------------------------------------
!! Function: box_flux_booles
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_grid
    
    real(dp) :: delta_x, delta_y, delta_z
    real(dp), allocatable :: f_x(:), g_xy(:), h_xyz(:)
    integer :: n_bins, i_x, i_y, i_z
    real(dp) :: x, y, z

    

    ! First we need to determine the distance between
    ! the lattice points at which the function to integrate
    ! will be evaluated $\Delta x$.

    

    ! bins:              1   2   3   4   5   6   7   8
    ! x interval:      |---|---|---|---|---|---|---|---|
    ! grid points:     1   2   3   4   5   6   7   8   9 
    
    ! interval length: |-------------------------------|
    !                  0                               depth
    ! delta x length:  |---|
    !                  0   delta_x
    ! The interval is the number of grid points - 1 
    ! Lets assign that to a variable called n_bins 
      n_bins = n_grid - 1

      delta_x = depth/n_bins
      delta_y = width/n_bins
      delta_z = height/n_bins


    ! Now we need to allocate memory for the arrays that will contain
    ! the evaluated function to integrate
    allocate(  f_x(1:n_grid))
    allocate( g_xy(1:n_grid))
    allocate(h_xyz(1:n_grid))

    ! Now we need to implement the do loop in the README file

    do i_x = 1, n_grid
        
          x = (i_x - 1) * delta_x
        do i_y = 1, n_grid
            
              y = (i_y - 1) * delta_y
            do i_z = 1, n_grid
                  z = (i_z - 1) * delta_z
                
                h_xyz(i_z) = flux_kernel(x, y, z, x_zero, y_zero)
            enddo
            
              g_xy(i_y) = booles_quadrature(h_xyz, delta_z)
        enddo
        
          f_x(i_x) = booles_quadrature(g_xy, delta_y)
    enddo
    
      flux = booles_quadrature(f_x, delta_x)
end function box_flux_booles

!-----------------------------------------------------------------------
!! Function: flux_kernel
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This is the function that we want to integrate. The evaluations of this function 
!! at certain points will fill the h_xyz array defined above.
!!
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x            real        x coordinate of the small integration volume
!! y            real        y coordinate of the small integration volume
!! y            real        z coordinate of the small integration volume
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function flux_kernel(x, y, z, x0, y0) result(k)
    implicit none
    real(dp), intent(in) :: x, y, z, x0, y0

    
    
      k = (1.0_dp)/(4*pi*((x+x0)**2 + (y-y0)**2 + z**2))
end function flux_kernel

!-----------------------------------------------------------------------
!! Function: large_x0_flux
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Imagine if x0 got really large, meaning we got really far away from the detector. 
!! We can approximate this as a point source and in this case we would put all the source 
!! that generates the flux at the center of the detector. So the resultant flux 
!! is the formula down below.
!!----------------------------------------------------------------------
!! Input:
!!
!! d            real        Depth of the rectangular nuclear reactor
!! w            real        Width of the rectangular nuclear reactor
!! h            real        Height of the rectangular nuclear reactor
!! x0           real        x coordinate of the detector's position
!! y0           real        y coordinate of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!!----------------------------------------------------------------------
real(dp) function large_x0_flux(d, w, h, x0, y0) result(flux)
    implicit none
    real(dp), intent(in) :: d, w, h, x0, y0
    
     flux = d * w * h/(4*pi*((x0 + d/2)**2 + (y0 - w/2)**2 + (h/2)**2))
end function large_x0_flux

!-----------------------------------------------------------------------
!! Subroutine: box_flux_monte_carlo
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This evaluates the box flux using monte carlo integration. Same physical 
!! parameters as before, but n_samples is how many times we want to take a random 
!! number in the interval of the integrand.
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine box_flux_monte_carlo(depth, width, height, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: flux, sigma_f
    
    real(dp) :: a(1:3), b(1:3), data(1:2)

    

    ! a is the lower integration limit in the x, y, z coordinates. 
    ! since the origin was placed at the corner of the nuclear reactor the
    ! lower limit is zero in all coordinates
    a = 0._dp

    ! b is the upper integration limit in the x, y, z coordinates.
    b(1) = depth
    b(2) = width
    b(3) = height

    ! This is the 'work array' and contains parameters
    ! (other than the sample point) needed to evaluate the function to
    ! integrate
    data(1) = x_zero
    data(2) = y_zero

    ! Notice that we're sending a function called flux_kernel_vector to
    ! the Monte Carlo subroutine. We need to define that function below. 
    call monte_carlo_quad(flux_kernel_vector, a, b, data, n_samples, flux, sigma_f)
end subroutine box_flux_monte_carlo

!-----------------------------------------------------------------------
!! Function: flux_kernel_vector
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function works with the interface defined in the quadrature module 
!! to make 2 arrays that hold the parameters of the detector and the random 
!! sampling points. k is going to hold the result of sending these parameters over 
!! to the function to be evaluated. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
! Because of the interface defined in the quadrature module the 
! Monte Carlo subroutine expects a kernel function that receives two
! arrays, the first one contains the sampling point, the second one
! contains the parameters needed to calculate the kernel. 
real(dp) function flux_kernel_vector(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0

     x = x_vector(1)
     y = x_vector(2)
     z = x_vector(3)
     x0 = data(1)
     y0 = data(2)

    ! We're going to use the function we already defined for the 
    ! Boole's integration.
     k = flux_kernel(x, y, z, x0, y0)
end function flux_kernel_vector


! As explained in the README, for Boole's method we need to calculate
! the flux of the solid box and subtract the flux from a solid sphere.
! Let's start defining the function that calculates the flux of a 
! solid sphere

!-----------------------------------------------------------------------
!! Function: sphere_flux_booles
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This calculates the flux of the sphere inside the reactor. This is 
!! important because the total flux of the reactor with the sphere inside 
!! is the flux of the box - the flux of the sphere
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! distance     real        Distance from the center of the reactor to the detector
!! radius       real        Radius of the spherical reactor
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function sphere_flux_booles(distance, radius, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: distance, radius
    integer, intent(in) :: n_grid

    real(dp) :: delta_r, delta_theta
    real(dp), allocatable :: f_r(:), g_rtheta(:)
    integer :: n_bins, i_r, i_theta
    real(dp) :: r, theta


    n_bins = n_grid - 1
    ! Base the rest of the function on the one for the solid box.
    ! Here we're integrating only two variables (r and theta).
    ! r is integrated from 0 to radius while theta is integrated
    ! from 0 to pi 

    ! distance is the distance from the center of the sphere 
    ! to the position of the detector (capital R in the README).
    ! It will be given as a input to this function and passed
    ! to the sphere_flux_kernel function defined below 
      delta_r = radius/n_bins
      delta_theta = pi/n_bins
      


    ! Now we need to allocate memory for the arrays that will contain
    ! the evaluated function to integrate
    allocate(  f_r(1:n_grid))
    allocate( g_rtheta(1:n_grid))
    

    ! Same procedure as before but with a double integral and integrating 
    ! With respect to r and theta instead of x, y or z. 
    do i_r = 1, n_grid
        
          r = (i_r - 1) * delta_r
        do i_theta = 1, n_grid
            
              theta = (i_theta - 1) * delta_theta
            
               g_rtheta(i_theta) = sphere_flux_kernel(r, theta, distance)
            enddo
           f_r(i_r) = booles_quadrature(g_rtheta, delta_theta)
    enddo 
    flux = booles_quadrature(f_r, delta_r)
    
      

end function sphere_flux_booles

!-----------------------------------------------------------------------
!! Function: sphere_flux_kernel
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! Contains function to be integrated. The evaluations of this function 
!! will be held inside the g_rtheta array. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! r_prime      real        r coordinate of the small integration volume
!! theta        real        theta coordinate of the small integration volume
!! big_r        integer     distance from the center of the sphere to the detector
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function sphere_flux_kernel(r_prime, theta, big_r) result(k)
    implicit none
    real(dp), intent(in) :: r_prime, theta, big_r
    
    
    k = (r_prime**2) * (sin(theta)/(2*( r_prime**2 + big_r**2 - 2 * r_prime * big_r * cos(theta))))
end function sphere_flux_kernel

!-----------------------------------------------------------------------
!! Function: total_flux_booles
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! The total flux of the reactor is the flux due to the box - the flux 
!! due to the sphere. This function calculates that and holds it in a variable 
!! called flux. 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_grid       integer     number of grid points (in each dimension) used the the quadrature
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the 3 dimensional integral
!-----------------------------------------------------------------------
real(dp) function total_flux_booles(depth, width, height, radius, x_zero, y_zero, n_grid) result(flux)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_grid

    real(dp) :: distance, box_flux, sphere_flux

    ! Now that we have a function to calculate the flux of the solid box and
    ! another one for the solid sphere we just need to use both functions 
    ! and calculate the difference.

    ! distance is the distance between the position of the detector (x_zero, y_zero)
    ! and the center of the sphere (which is also the center of the box) 
    ! The distance between 2 points in xyz plane is the square root of 
    ! (x2- x1)^2 + (y2- y1)^2 + (z2 - z1)^2. 
    ! The detector's position is (x0, y0, 0) while the position of the center 
    ! of the sphere is (-depth/2, width/2, height/2). This is because 
    ! respect to the detector, depth is negative.
     
     distance = sqrt((height/2)**2 + (width/2 - y_zero)**2 + (-1.0 * depth/2 - x_zero)**2)
     box_flux = box_flux_booles(depth, width, height, x_zero, y_zero, n_grid) 
     sphere_flux = sphere_flux_booles(distance, radius, n_grid)
     flux = box_flux - sphere_flux

end function total_flux_booles

! As explained in the README the Monte Carlo approach is simpler.
! We just need to define a new kernel function that is zero if the 
! sampling point is inside the sphere and the original kernel if
! the sampling point is outside of the sphere.

! Again, this new kernel will take to arrays, one with the coordinates
! of the sampling point and one with all the other needed parameters
! this time it's more than just the position of the detector.

!-----------------------------------------------------------------------
!! Function: hollow_box_flux_kernel
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Similar thing as before with the flux kernel vector. We made an interface 
!! that expects 2 arrays. The x vector here will hold the randomly sampled points 
!! while the data array will hold the physical parameters of the sphere, box and 
!! detector this time. We will then run a check. This check will say that 
!! the result is 0 if inside the sphere because there is no flux inside the sphere 
!! or it is equivalent to the kernel function that we defined earlier if outside the sphere. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! x_vector     real        array containing the x, y, z, coordinates of the integration volume
!! data         real        work array containing the sphere's radius and x, y coordinates of the detector's position
!!----------------------------------------------------------------------
!! Output:
!!
!! k            real        kernel to be integrated
!-----------------------------------------------------------------------
real(dp) function hollow_box_flux_kernel(x_vector, data) result(k)
    implicit none
    real(dp), intent(in) :: x_vector(:), data(:)

    real(dp) :: x, y, z, x0, y0, radius, depth, width, height
    real(dp) :: distance_to_center

      x = x_vector(1)
      y = x_vector(2)
      z = x_vector(3)

      x0 = data(1)
      y0 = data(2) 
      radius = data(3) 
      depth = data(4) 
      width = data(5) 
      height = data(6)

    ! We need to determine whether or not the sampling point is inside the 
    ! sphere. For that you can calculate the distance from the sampling point
    ! (THIS IS NOT THE POSITION OF THE DETECTOR) and the center of the sphere 
    ! and compare it with the sphere's radius
    ! Remember that the origin was located at one corner of the box.
    ! Lets use that distance formula again. Lets take the origin to the front left corner or front right corner 
    ! of the box. This means the height of the box is positive, as is the width and the depth. This is because 
    ! unlike the detector case, the center of the sphere and the sampled point lie in the same quadrant
     
     distance_to_center = sqrt((depth/2 - x)**2 + (width/2 - y)**2 + (height/2 - z)**2)

    ! Now we just need to check if that point is inside the sphere.

      if (distance_to_center < radius) then
          k = 0._dp
      else
          k = flux_kernel(x, y, z, x0, y0)
      end if

end function hollow_box_flux_kernel


!-----------------------------------------------------------------------
!! Subroutine: hollow_box_flux_mc
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This finds the result of the monte carlo integration by passing 2 arrays 
!! to another function to check the distance between that randomly sampled point 
!! and the radius and by passing the result of that and the physical parameters 
!! to the monte carlo subroutine in the quadrature module to calculate the result 
!! of the integral and the error. 
!! 
!!----------------------------------------------------------------------
!! Input:
!!
!! depth        real        Depth of the rectangular nuclear reactor
!! width        real        Width of the rectangular nuclear reactor
!! height       real        Height of the rectangular nuclear reactor
!! radius       real        Radius of the hollow sphere
!! x_zero       real        x coordinate of the detector's position
!! y_zero       real        y coordinate of the detector's position
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! flux         real        Result of the Monte Carlo integral
!! sigma_f      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine hollow_box_flux_mc(depth, width, height, radius, x_zero, y_zero, n_samples, flux, sigma_f)
    implicit none
    real(dp), intent(in) :: depth, width, height, radius, x_zero, y_zero
    integer, intent(in) :: n_samples
    real(dp), intent(out) ::  flux, sigma_f

    real(dp) :: a(1:3), b(1:3), data(1:6)

    ! Base the rest of the function on the one from the basic part of the project
    ! Remember that the 'work array' data now contains more information 
    ! 

      a = 0.0_dp
      b(1) = depth
      b(2) = width
      b(3) = height

      data(1) = x_zero 
      data(2) = y_zero 
      data(3) = radius 
      data(4) = depth 
      data(5) = width 
      data(6) = height
      
    
    
      call monte_carlo_quad(hollow_box_flux_kernel, a, b, data, n_samples, flux, sigma_f)
end subroutine hollow_box_flux_mc

end module neutron_flux
