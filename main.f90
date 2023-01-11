! Program: nuclear_reactor
! By: Nathan Crawford
!-----------------------------------------------------------------------------
! The first evaluates the neutron flux of a solid box neutron source, 
! in this case it is a reactor. 
! The specific dimensions of the reactor, the coordinates of the detector 
! That detects the neutrons are all user input as well as the 
! x interval and x increment values. 

! It does this using Boole's method and Monte Carlo integration 
! The second part assumes a hollow sphere inside the reactor 
! and finds the total flux of that configuration. 
! This time the minimum and maximum radius of the sphere 
! are user input as well as the radius increment. 

! You're free to give different values when running the code, the ones 
! bellow are just a suggestion that worked fine for me.

! Basic part of the project
! depth = 40.0_dp
! width = 100.0_dp
! height = 60.0_dp
! y_zero = 30._dp
! x_min = 5_dp
! x_max = 200._dp
! x_step = 5._dp
! n_grid = 25
! m_samples = 10000

! Advanced part of the project
! x_zero = 80._dp
! r_min = 1.0_dp
! r_max = 19._dp
! r_step = 1.0_dp
!-----------------------------------------------------------------------------
program nuclear_reactor
use types 

use read_write, only : read_input, write_neutron_flux, read_advanced_input, write_advanced_flux

implicit none

real(dp) :: depth, width, height, y_zero, x_min, x_max, x_step
integer :: n_grid, m_samples
real(dp) :: x_zero, r_min, r_max, r_step


! Basic part of the project
call read_input(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)
call write_neutron_flux(depth, width, height, y_zero, x_min, x_max, x_step, n_grid, m_samples)

! Advanced part of the project
  call read_advanced_input(x_zero, r_min, r_max, r_step)
  call write_advanced_flux(depth, width, height, x_zero, y_zero, r_min, r_max, r_step, n_grid, m_samples)

end program nuclear_reactor
