#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This routine determines the surface pressure (P), surface air temperature (T), 
weighed mean temperature (Tm), and zenith hydrostatic delay (ZHD) from binary coefficient files
As available from:
https://github.com/pjmateus/hgpt_model (release v1.0)
press_grid.bin; temp_grid.bin; and tm_grid.bin

It is admitted that the binary files with the coefficients are in the same directory as this script.
In alternative you can define the "coeffiles" variable

The epoch can be an array of size 1, and in this case is the Modified Julian Date (MJD) 
or can be an array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)

All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude

Reference for HGPT:
An ERA5-based hourly global temperature and pressure (HGTP) model (submitted to Remote Sensing, MDPI)

INPUT:
             dt : if size(dt)=1 => modified julian date
                  if size(dt)=6 => year, month, day, hour, min, sec
             x0 : ellipsoidal longitude (degrees)   
             y0 : ellipsoidal latitude (degrees)         
             z0 : height (m)
        z0_type : ‘orth’ for orthometric height or ‘elli’ for ellipsoidal height
OUTPUT:
              P : surface pressure valid at (x0, y0, z0), in hPa
              T : surface air temperature valid at (x0, y0, z0), in Kelvins
             Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
            ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
--------------------------------------------------------------------------
 Example:       
	y0 = 38.5519
	x0 = -9.0147
	z0 = 25
	dt = 58119.5   or   dt = np.array([2018, 1, 1, 12, 0, 0])
    P, T, Tm, ZHD = hgpt(dt, x0, y0, z0, 'orth')
--------------------------------------------------------------------------
written by Pedro Mateus (2020/01/15)
Instituto Dom Luiz (IDL), Faculdade de Ciências, Universidade de Lisboa, 1749-016 Lisboa, Portugal
pjmateus@fc.ul.pt

Dependencies:
https://pypi.org/project/julian/
pip install julian
"""

import numpy as np
import julian
from datetime import datetime
from scipy.interpolate import RegularGridInterpolator

def hgpt(dt, x0, y0, z0, z0_type):
       
       # Grid files location
       coeffiles='' # put '/' or '\' at the end

       # Constants
       row = 721
       col = 1440
       p1  = 365.250
       p2  = 182.625
       p3  = 91.3125
       deg2rad = np.pi/180.0

       # Geographic coordinates ( equal to ERA5 )
       lon = np.linspace(-179.75, 180, col)
       lat = np.linspace(-90, 90, row)     

       # Modified Julian date 
       if np.size(dt) == 6:
              # Input: Gregorian calendar              
              mjd = julian.to_jd(datetime(np.int(dt[0]),np.int(dt[1]),np.int(dt[2]), \
                                          np.int(dt[3]),np.int(dt[4]),np.int(dt[5])), fmt='mjd')                       
              hour = np.int(dt[3]) 
       elif np.size(dt) == 1:
              # Input: Modified Julian date 
              gre = julian.from_jd(dt, fmt='mjd')
              mjd = dt
              hour = np.int(np.around(gre.hour))
       else:
              raise NameError('Use 1) Modified Julian Date (MJD) or 2) Gregorian date (y,m,d,HH,MM,SS).')
       
       # Open and read the surface air temperature coefficients file       
       fid = open(coeffiles+'temp_grid.bin', 'rb')
       fid.seek((row*col*26)*hour, 0)       
       y_intercept = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       slope = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       a1 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       f1 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a2 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       f2 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a3 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       f3 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       fid.close()    

       # Bilinear interpolation       
       F = RegularGridInterpolator((lat, lon), y_intercept, method='linear'); a = F(np.array([y0, x0]))[0] 
       F = RegularGridInterpolator((lat, lon), slope, method='linear'); b = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), a1, method='linear'); amp1 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), f1, method='linear'); pha1 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), a2, method='linear'); amp2 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), f2, method='linear'); pha2 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), a3, method='linear'); amp3 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), f3, method='linear'); pha3 = F(np.array([y0, x0]))[0]

       # Surface air temperature model
       T = a + b*(mjd - 51178) + amp1*np.cos(2*np.pi*(mjd - 51178)/p1+pha1) + \
                                 amp2*np.cos(2*np.pi*(mjd - 51178)/p2+pha2) + \
                                 amp3*np.cos(2*np.pi*(mjd - 51178)/p3+pha3)

       # Open and read the surface pressure coefficients file
       fid = open(coeffiles+'press_grid.bin', 'rb')
       fid.seek((row*col*20)*hour, 0)       
       y_intercept = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       slope = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       a1 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       f1 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a2 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       f2 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       fid.close()    
       
       # Bilinear interpolation
       F = RegularGridInterpolator((lat, lon), y_intercept, method='linear'); a = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), slope, method='linear'); b = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), a1, method='linear'); amp1 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), f1, method='linear'); pha1 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), a2, method='linear'); amp2 = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), f2, method='linear'); pha2 = F(np.array([y0, x0]))[0]
       
       # Surface pressure model
       P = a + b*(mjd - 51178) + amp1*np.cos(2*np.pi*(mjd - 51178)/p1+pha1) + \
                                 amp2*np.cos(2*np.pi*(mjd - 51178)/p2+pha2)
                       
       # Open and read the Tm coefficients and undulation file
       fid = open(coeffiles+'tm_grid.bin', 'rb')   
       y_intercept = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       slope = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       orography = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       undu = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F') 
       fid.close()   
       
       # Bilinear interpolation
       F = RegularGridInterpolator((lat, lon), y_intercept, method='linear'); a = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), slope, method='linear'); b = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), orography, method='linear'); geo_height = F(np.array([y0, x0]))[0]
       F = RegularGridInterpolator((lat, lon), undu, method='linear'); N = F(np.array([y0, x0]))[0]

       # Zenith hydrostatic delay (ZHD), Saastamoinen model
       if z0_type=='orth':
              H_orth = z0
       elif z0_type=='elli':
              H_orth = z0 - N
       else:
              raise NameError('Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).')

       # Correction to P and T (see Guochanf Xu, GPS Theory, Algorithms and Applications, 2nd Edition, page 56)
       P = (P*100.0 * (1.0 - 0.0065/T * (H_orth - geo_height))**5.2559)/100.0
       T = T - 0.0065*(H_orth - geo_height)
       
       # Weight mean temperature, Tm
       Tm = a + b*T
       
       # ZHD using the Saastamoinen Model (see Saastamoinen, 1973)
       ZHD = (0.0022768 * P)/(1 - 0.0026*np.cos(2*deg2rad*y0)-0.00000028*H_orth)
              
       return P, T, Tm, ZHD
