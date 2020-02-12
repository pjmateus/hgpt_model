# hgpt model
**An ERA5-based hourly global temperature and pressure model**

Hourly global temperature and pressure (HGTP) model is based on the full spatial and temporal resolution of the new ERA5 reanalysis produced by the ECMWF. The HGPT is based on the time-segmentation concept and uses three periodicities for surface air temperature and two for surface pressure. The weighted mean temperature is determined using 20-years of monthly data to contemplate its seasonality and geographic variability. We also introduced a linear trend to account for the global climate change scenario.

The model was developed at the Dom Luiz Institute (IDL), Faculty of Sciences of the University of Lisbon (FCUL), by Pedro Mateus, João Catalão and Virgílio Mendes. Also by Giovanni Nico from the Istituto per le Applicazioni del Calcolo (IAC), Consiglio Nazionale delle Ricerche (CNR), 70126 Bari, Italy. 

<img src="https://github.com/pjmateus/hgpt/blob/master/logos.png" width="450">

The same code is available in three programming languages, Fortran, Matlab and Python. The header contains guidelines for running each of these codes.

Simple Fortran code to call the hgpt subroutine (in file hgpt.f90) 
```Fortran
program call_hgpt_model
	implicit None
	real :: x0, y0, z0, P, T, Tm, ZHD
	real*8, dimension(1) :: dt
	y0 = 38.5519    ! Latitude, degrees
	x0 = -9.0147    ! Longitude, degrees
	z0 = 25         ! Orthometric height, m
	dt(1) = 58119.5 ! MJD
	call hgpt(dt, size(dt), x0, y0, z0, 'orth', P, T, Tm, ZHD)
	write(*,*) P, T, Tm, ZHD
end program call_hgpt_model             
```
Save this file with the name call_hgpt_model.f90 and compile with gfortran hgpt.f90 call_hgpt_model.f90 -o call_hgpt_model.exe
Module constructs can also be easily implemented.

Matlab code to call the hgpt function 
```Matlab
y0 = 38.5519;  % Latitude, degrees
x0 = -9.0147;  % Longitude, degrees
z0 = 25;       % Orthometric height, m
dt = 58119.5;  % MJD
[P, T, Tm, ZHD] = hgpt(dt, x0, y0, z0, 'orth')
```
Locate the function hgpt.m and the binary grid files (or set the coeffiles variable in hgpt.m file) and run this code.

Python code to call the hgpt function 
```Python
y0 = 38.5519  # Latitude, degrees
x0 = -9.0147  # Longitude, degrees
z0 = 25       # Orthometric height, m
dt = 58119.5  # MJD
P, T, Tm, ZHD = hgpt(dt, x0, y0, z0, 'orth')
```
