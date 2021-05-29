# hgpt model
**An ERA5-based hourly global pressure and temperature model**

Hourly global pressure and temperature and pressure (HGTP) model is based on the full spatial and temporal resolution of the new ERA5 reanalysis produced by the ECMWF. The HGPT is based on the time-segmentation concept and uses three periodicities for surface air temperature and two for surface pressure. The weighted mean temperature is determined using 20-years of monthly data to contemplate its seasonality and geographic variability. We also introduced a linear trend to account for the global climate change scenario.

**Version 2 is also based on the time-segmentation concept. We introduced the relative humidity (RH, %), the zenith wet delay (ZWD, m), and the precipitable water vapor (PWV, m). The zenith total delay (ZTD) is easily calculated by ZTD = ZHD + ZWD, and the water vapor (e) or saturated water vapor (es) can be calculated using the temperature and relative humidity, e.g., using the Wexler formulation (available in the code).** 

The model was developed at the Dom Luiz Institute (IDL), Faculty of Sciences of the University of Lisbon (FCUL), by Pedro Mateus, João Catalão and Virgílio Mendes. Also by Giovanni Nico from the Istituto per le Applicazioni del Calcolo (IAC), Consiglio Nazionale delle Ricerche (CNR), 70126 Bari, Italy.
In the development of the second version we also count on the help of Sandra M. Plecha from FCUL/IDL.

<img src="https://github.com/pjmateus/hgpt/blob/master/logos.png" width="450">

The same code is available in three programming languages, Fortran, Matlab and Python. The header contains guidelines for running each of these codes. 
**We advise using the version 2 code as it offers significant improvements, mainly in the interpolation process, eliminating intrinsic function as griddedInterpolant or RegularGridInterpolator. We also eliminate the use of mjuliandate in Matlab once it requires Aerospace Toolbox.**

**To download the binary grid files use the [release section (v1.0)](https://github.com/pjmateus/hgpt_model/releases)**

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
Save this file with the name call_hgpt_model.f90 and compile with "gfortran hgpt.f90 call_hgpt_model.f90 -o call_hgpt_model.exe". Module constructs can also be easily implemented.

Matlab code to call the hgpt function [![View hgpt_model on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/74247-hgpt_model)
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
from hgpt import hgpt
y0 = 38.5519  # Latitude, degrees
x0 = -9.0147  # Longitude, degrees
z0 = 25       # Orthometric height, m
dt = 58119.5  # MJD
P, T, Tm, ZHD = hgpt(dt, x0, y0, z0, 'orth')
```
Locate the function hgpt.py and the binary grid files (or set the coeffiles variable in hgpt.py file) and run this code.
Requirements:

You need Python 3.2 or later and julian 0.14. You can install julian like this:
```
$ pip install julian
```
**If you have any questions do not hesitate to contact me by email pjmateus@fc.ul.pt**

Don't forget to cite:

Mateus, P.; Catalão, J.; Mendes, V.B.; Nico, G. An ERA5-Based Hourly Global Pressure and Temperature (HGPT) Model. Remote Sens. 2020, 12, 1098; https://doi.org/10.3390/rs12071098 
