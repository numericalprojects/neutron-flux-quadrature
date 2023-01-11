!-----------------------------------------------------------------------
!Module: quadrature
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! You can think of this as the math module. The subroutines and functions 
!! in this module will take care of the mathematics without even knowing 
!! about or considering the physics. The booles quadrature and booles rule 
!! functions evaluate the integral using booles rule. 
!! Likewise the monte_carlo_quad subroutine evaluates the integral using 
!! monte carlo integration.
!! If you look closely at the monte carlo subroutine you will notice 
!! it is generalized for any number of dimensional integrals. 
!!----------------------------------------------------------------------
!! Included subroutines:
!!
!! 
!! monte_carlo_quad 
!!----------------------------------------------------------------------
!! Included functions:
!!
!! booles_quadrature
!! booles_rule
!-----------------------------------------------------------------------
module quadrature
use types

implicit none

private
public :: booles_quadrature, monte_carlo_quad

!-----------------------------------------------------------------------
!Interface: func
!-----------------------------------------------------------------------
!! This defines a new type of procedure in order to allow callbacks
!! in the Monte Carlo quadrature subroutine of an arbitrary function that is given
!! as input and declared as a procedure
!!
!! The arbitrary function receives two rank 1 arrays of arbitrary size.
!! The first array contains an n-dimensional vector representing the
!! point sampled by the Monte Carlo method. The second is a "work array"
!! that contains parameters  necessary to calculate the function to be
!! integrated.
!!----------------------------------------------------------------------
interface
    real(dp) function func(x, data)
        use types, only : dp
        implicit none
        real(dp), intent(in) :: x(:), data(:)
        
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Function: booles_quadrature
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function takes in an array filled with evaluations of a function 
!! at the points from 1 to the number of lattice points and sends slices 
!! of that array to a different booles function to perform the actual 
!! integration. 
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_quadrature(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size, i

    fx_size = size(fx)

    ! As the diagram below shows, only certain number of grid points
    ! fit the scheme of Boole's quadrature. 

    ! |--interval 1---|--interval 2---|--interval 3---|
    ! 1   2   3   4   5   6   7   8   9   10  11  12  13
    ! |---|---|---|---|---|---|---|---|---|---|---|---|
    ! x0  x1  x2  x3  x4
    !                 x0  x1  x2  x3  x4
    !                                 x0  x1  x2  x3  x4

    ! Using the modulo operator we can see if the size of the array is correct 
    ! by checking the remainder left over by dividing fx_size - 1, which according 
    ! to the above diagram is always a multiple of 4,by 4. If there is any remainder 
    ! then the array is an incorrect size.
      
      if (modulo((fx_size-1),4) /= 0) then 
          print *, 'fx array size in booles_quadrature has to be a multiple of 4'
          stop
      endif

    ! We could implement the full integration here, however to make a cleaner,
    ! easy to read (and debug or maintain) code we will define a smaller
    ! function that returns Boole's five point rule and pass slices (1:5), (5:9),
    ! (9:13), ... of fx to such function to then add all the results. 

    s = 0._dp
    
    ! In terms of the actual integration, we want to send the slices as described above 
    ! up to fx_size. To do that we can use the syntax i:i+4 where i is an element of an array. 
    ! This sends the current element to the element 4 steps ahead of it, guaranteeing the correct 
    ! section of the array being sent for evaluation.
      
      do i = 1, fx_size - 1, 4
          s = s + booles_rule(fx(i:i+4), delta_x)
      enddo
end function booles_quadrature

!-----------------------------------------------------------------------
!! Function: booles_rule
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function performs the actual integration. It takes in the delta x 
!! or dx and slices of that fx array with the evaluations and uses it to 
!! evaluate the boole's rule formula. 
!! ----------------------------------------------------------------------
!! Input:
!!
!! fx           real        Array containing the evaluated function
!! delta_x      real        Distance between the evaluation points
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Boole's quadrature
!-----------------------------------------------------------------------
real(dp) function booles_rule(fx, delta_x) result(s)
    implicit none
    real(dp), intent(in) :: fx(1:), delta_x

    integer :: fx_size
    real(dp) :: fx0, fx1, fx2, fx3, fx4

    fx_size = size(fx)

    ! Let's make an additional test to make sure that the array
    ! received has 5 and only 5 points 

      if (fx_size /= 5) then
          print *, "Did not receive 5 points. Please adjust quadrature algorithm." 
          stop
      endif
    
    ! To improve readability lets assign the values of the function at specific points 
    ! to real variables.
      
      fx0 = fx(1)
      fx1 = fx(2)
      fx2 = fx(3) 
      fx3 = fx(4) 
      fx4 = fx(5)
   
   ! Formula for Booles rule being held in variable s.  
      
      s = (delta_x * 2.0_dp/45.0_dp)*(7 * fx4 + 32 * fx3 + 12 * fx2 + 32 * fx1 + 7 * fx0)
end function booles_rule

!-----------------------------------------------------------------------
!! Subroutine: monte_carlo_quad
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This performs monte carlo integration. This is much more straightforward. 
!! It samples random points using the pseudorandom number generator which 
!! has a uniform sampling distribution from 0 to 1, scales that to be between 
!! the proper interval, and sums it up and finds the integral value using the 
!! Monte Carlo formula. 
!! ----------------------------------------------------------------------
!! Input:
!!
!! f            procedure   function to be integrated
!! a            real        array containing the lower limits of the integral
!! b            real        array containing the upper limits of the integral
!! data         real        array containing parameters necessary to calculate the function f
!! n_samples    integer     number of sample points in the Monte Carlo integration
!!----------------------------------------------------------------------
!! Output:
!!
!! s            real        Result of the Monte Carlo integral
!! sigma_s      real        Estimate of the uncertainty in the Monte Carlo integral
!-----------------------------------------------------------------------
subroutine monte_carlo_quad(f, a, b, data, n_samples, s, sigma_s)
    implicit none
    procedure(func) :: f
    real(dp), intent(in) :: a(:), b(:), data(:)
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: s, sigma_s

    integer :: i, vector_size
    real(dp), allocatable :: x_vector(:), fx(:), f_squared(:) ! ...you might need to declare other arrays here
    real(dp) :: avg, summ_squared, avg_squared, volume, variance

    vector_size = size(a)

    ! We're defining a Monte Carlo routine that works for an arbitrary number of 
    ! dimensions in the integral (Remember, that's the advantage of Monte Carlo integration,
    ! it's very efficient for high dimensional integrals)

    ! Since a and b give the lower and upper limits they need to have the same size.
    ! Make a check to see if they do have the same size

      if (size(a) /= size(b)) then 
           print *, 'a and b arrays in monte_carlo_quad have to be the same size'
           stop       
      endif

    ! Here we allocate memory for the vector containing the sample points and 
    ! for a vector that contains the evaluated function
    allocate(x_vector(1:vector_size))
    allocate(fx(1:n_samples)) 
    allocate(f_squared(1:n_samples))
    
    
    s = 0.0_dp 
    
    summ_squared = 0.0_dp
    ! This do loop evaluates the sum of all the random numbers
    do i=1,n_samples
        

        call random_number(x_vector) !generates an array with random numbers in the [0,1) interval
        x_vector = a + x_vector*(b-a) !rescaling to the integration volume [a,b)
        fx(i) = f(x_vector,data) 
        s = s + fx(i)
    enddo 
    
    avg = s/n_samples
    volume = 1.0_dp 
   
   ! This do loop evaluates the volume for any number of dimensions. It is generalized for any 
   ! dimensional integral
   
   do i = 1, size(a) 
        volume = (b(i) - a(i)) * volume
   end do
   
   ! This evaluates the f^2 term in the Monte Carlo formula. 
   ! It is the evaluations squared at every evaluation and summed up. 
   
   do i = 1, n_samples 
        f_squared(i) = fx(i)**2 
        summ_squared = summ_squared + f_squared(i) 
    end do     
   
   ! Average of the squared terms and the intrinsic variance term for the uncertainty
    
    avg_squared = summ_squared/n_samples
    variance = avg_squared - (avg**2) 
   
   ! Final results
   
   sigma_s = sqrt(variance/n_samples) * (volume)
    s = volume * avg  
    
    
end subroutine monte_carlo_quad

end module quadrature
