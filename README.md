## Goal of the Program

Consider a reactor of depth D along the $x$ axis, width W along the $y$ axis and height H along the $z$ axis. Assume the reactor emits neutrons uniformly throughout its volume and that the emission rate is 1 neutron per second per cubic meter. Now place a neutron detector a distance $x_0$ in front of the reactor, and over a distance $y_0$ along the front from one corner (either corner, it does not matter), as seen in the image below. 

![reactor](https://user-images.githubusercontent.com/89489977/210156596-133cc5f8-d3ef-4521-a14a-e7723ad9447e.png)

The goal of this part of the project is to plot the total flux as a function of the distance $x_0$ using 2 different methods and then compare. 

The flux from a small volume $dV = dxdydz$ at a distance r is $\frac{dxdydz}{4πr^2(x, y, z, x_0, y_0)}$ where $r(x, y, z, x_0, y_0)$ is the distance 
from each small volume $(x, y, z)$ to the reactor $(x_0, y_0)$. 

Now you can place the origin anywhere but it is convenient to place the origin at one of the front corners. So the integral then becomes 

![integral](https://user-images.githubusercontent.com/89489977/210156753-3a17f425-81bd-414e-a552-4d2c3c796b6e.PNG)

For this integral I am going to use Boole's rule and Monte Carlo and compare the results. Monte Carlo is obviously more 
accurate for multiple integrals such as this.  

For the next part let us assume a modified reactor configuration in which a sphere of radius $r$ is inside the center of the reactor. 

![sphere](https://user-images.githubusercontent.com/89489977/210156870-a94d8237-a18a-4ab5-b05c-b9924f850cb9.png) 

Taking advantage of the spherical symmetry, the integral for Boole's now becomes: 

![integral2](https://user-images.githubusercontent.com/89489977/210156893-be812c03-7a5e-468c-a759-afa32d142018.PNG)

However, with Monte Carlo we don't need this since we can just sample points and set it to 0 if the point is in the sphere. 

## Formulas 
 Boole's: 
  
  ![image](https://user-images.githubusercontent.com/89489977/210156928-1e9800e5-4ac8-4272-baa1-170196af9541.png) 
  
 Monte Carlo: 
 
  ![image](https://user-images.githubusercontent.com/89489977/210156944-81fa9a9c-f526-4ecf-998f-20ef4c96e51c.png)

Error for Monte Carlo: 

 $σ_A = \frac{V}{\sqrt{M}}σ_f$ 
 where 
 
 ![image](https://user-images.githubusercontent.com/89489977/210156977-36c55f5a-a4df-4495-8595-8a0538f72774.png)

## Running the program 

It would be quite beneficial to you if you had a Linux system because it would enable you to use the makefile included. 

If this is the case then what you do is open a terminal, use the cd command to change to this directory. 

Then type make. 

You'll see some gfortran commands being executed. All of this has created an exectuable file called nuclear-reactor. 

In order to run this executable you type ./nuclear-reactor. Make sure not to put any spaces! 

Next comes the most important part, your user input! 

For the physical parameters, please enter positive real numbers. There are checks in the code itself 
to stop the program and let you know if your input is valid but know that the program expects positive 
real values. 

For the number of lattice points in Boole's rule, remember that the number of points _minus 1_ has to be 
a multiple of 4 and that the number of Monte Carlo integration points has to be positive and an integer(whole number)

Next open the plots.ipynb to look at the analysis of the results. 
