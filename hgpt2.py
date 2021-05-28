#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This routine determines the surface pressure (P, hPa), surface air temperature (T, K), relative humidity (RH, %)
weighed mean temperature (Tm, K), zenith hydrostatic delay (ZHD, m), zenith wet delay (ZWD, m), and
precipitable water vapor (PWV, m) from binary coefficient files
As available from:
https://github.com/pjmateus/hgpt_model (release v2.0)
press_grid.bin; temp_grid.bin; tm_grid.bin; and rh_grid.bin

It is admitted that the binary files with the coefficients are in the same directory as this script.
In alternative you can define the "coeffiles" variable

The epoch can be an array of size 1, and in this case is the Modified Julian Date (MJD)
or can be an array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)

All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude

Reference for HGPT and HGPT2:
Mateus, P.; Catalão, J.; Mendes, V.B.; Nico, G. An ERA5-Based Hourly Global Pressure and Temperature (HGPT) Model.
Remote Sens. 2020, 12, 1098. https://doi.org/10.3390/rs12071098
HGPT2: an ERA5-based global model to estimate relative humidity (Remote Sensing, MDPI)

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
             RH : relative humidity valid at (x0, y0, z0), in %
             Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
            ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
            ZWD : zenith wet delay, valid at (x0, y0, z0), in meters
            PWV : precipitable water vapor, valid at (x0, y0, z0), in meters
--------------------------------------------------------------------------
 Example:
	y0 = 38.5519
	x0 = -9.0147
	z0 = 25
	dt = 58119.5   or   dt = np.array([2018, 1, 1, 12, 0, 0])
    P, T, RH, Tm, ZHD, ZWD, PWV = hgpt2(dt, x0, y0, z0, 'orth')
    (1020.262,
      286.589,
       74.908,
      277.236,
        2.324,
        0.117,
        0.018)
--------------------------------------------------------------------------
written by Pedro Mateus (2021/05/15)
Instituto Dom Luiz (IDL), Faculdade de Ciências, Universidade de Lisboa, 1749-016 Lisboa, Portugal
pjmateus@fc.ul.pt

Dependencies:
https://pypi.org/project/julian/
pip install julian
"""
import numpy as np
import julian
from datetime import datetime

def es_wexler(T, P):
       ''' Saturation Water Vapor Pressure (es) using Wexler formulation with new coefficients (adjusted for ITS-90)
       INPUT : T = Temperature (in Kelvins)
               P = atmospheric Pressure (in hPa)
       OUTPUT: es= saturation water vapor pressure (in hPa)
       References:
       1)
       Wexler, A. Vapor Pressure Formulation for Water in Range 0 to 100 Degrees C. A Revision.
       J. Res. Natl. Bur. Stand. 1976, 80A, 775–785.
       2)
       Wexler, A. Vapor Pressure Formulation for Ice. J. Res. Natl. Bur. Stand. 1977, 81A, 5–20.
       3)
       Hardy, B. ITS-90 Formulations for Water Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in Range -100 to +100 C.
       In Proceedings of the Third International Symposium on Humidity and Moisture;
       UK National Physical Laboratory (NPL): Teddington, UK, April 6 1998; pp. 1–8.'''

       if T >= 273.15:
              # Saturation Vapor Pressure over Water
              g0 =-2.8365744*10**3
              g1 =-6.028076559*10**3
              g2 = 1.954263612*10**1
              g3 =-2.737830188*10**-2
              g4 = 1.6261698*10**-5
              g5 = 7.0229056*10**-10
              g6 =-1.8680009*10**-13
              g7 = 2.7150305
              es = 0.01 * np.exp(g0*T**-2 + g1*T**-1 + g2 + g3*T + g4*T**2 + g5*T**3 + g6*T**4 + g7*np.log(T))

              # Enhancement Factor coefficients for Water 0 to 100°C
              A0 =-1.6302041*10**-1
              A1 = 1.8071570*10**-3
              A2 =-6.7703064*10**-6
              A3 = 8.5813609*10**-9
              B0 =-5.9890467*10**1
              B1 = 3.4378043*10**-1
              B2 =-7.7326396*10**-4
              B3 = 6.3405286*10**-7
       else:
              # Saturation Vapor Pressure over Ice
              k0 =-5.8666426*10*3
              k1 = 2.232870244*10**1
              k2 = 1.39387003*10**-2
              k3 =-3.4262402*10**-5
              k4 = 2.7040955*10**-8
              k5 = 6.7063522*10**-1
              es = 0.01 * np.exp(k0*T**-1 + k1 + k2*T + k3*T**2 + k4*T**3 + k5*np.log(T))

              # Enhancement Factor coefficients for Ice –100 to 0°C
              A0 =-6.0190570*10**-2
              A1 = 7.3984060*10**-4
              A2 =-3.0897838*10**-6
              A3 = 4.3669918*10**-9
              B0 =-9.4868712*10**1
              B1 = 7.2392075*10**-1
              B2 =-2.1963437*10**-3
              B3 = 2.4668279*10**-6

       # Enhancement Factor
       alpha = A0 + A1*T + A2*T**2 + A3*T**3
       beta = np.exp(B0 + B1*T + B2*T**2 + B3*T**3)
       f = np.exp( alpha*(1-es/P) + beta*(P/es-1) )
       return es * f

def hgpt2(dt, x0, y0, z0, z0_type):

       # Grid files location
       coeffiles='' # put '/' or '\' at the end

       # Constants
       row = 721
       col = 1440
       p1  = 365.250
       p2  = 182.625
       p3  = 91.3125

       # Geographic coordinates ( equal to ERA5 )
       lon = np.linspace(-180, 179.75, col)
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

       # Finding indexes for bilinear interpolation
       # x-location
       indx = np.argsort( np.sqrt( (lon-x0)**2 ) )
       ix1 = indx[ 0 ]; ix2 = indx[ 1 ]
       x1 = lon[ ix1 ]; x2 = lon[ ix2 ]
       x = [ix1, ix1, ix2, ix2]
       # y-location
       indy = np.argsort( np.sqrt( (lat-y0)**2 ) )
       jy1 = indy[ 0 ]; jy2 = indy[ 1 ]
       y1 = lat[ jy1 ]; y2 = lat[ jy2 ]
       y = [jy1, jy2, jy1, jy2]
       # xy-distances (weights)
       dx1y1= 1/np.sqrt( (x1 - x0)**2 + (y1 - y0)**2 )
       dx1y2= 1/np.sqrt( (x1 - x0)**2 + (y2 - y0)**2 )
       dx2y1= 1/np.sqrt( (x2 - x0)**2 + (y1 - y0)**2 )
       dx2y2= 1/np.sqrt( (x2 - x0)**2 + (y2 - y0)**2 )
       dxy = np.array([dx1y1, dx1y2, dx2y1, dx2y2], dtype=np.float64)
       if np.any(np.isinf(dxy)) == True:
              # Exact point grid
              dxy = np.array([1,0,0,0], dtype=np.float64)

       # ******************************************************************
       # Open and read the surface air temperature coefficients file
       # ******************************************************************
       fid = open(coeffiles+'temp_grid.bin', 'rb')
       fid.seek((row*col*26)*hour, 0)
       a0 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       b0 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       a1 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f1 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a2 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f2 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a3 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f3 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       fid.close()

       # Surface air temperature model
       fun_t = lambda a0, b0, a1, f1, a2, f2, a3, f3: a0 + b0*(mjd - 51178) + a1*np.cos(2*np.pi*(mjd - 51178)/p1+f1) + \
       a2*np.cos(2*np.pi*(mjd - 51178)/p2+f2) + a3*np.cos(2*np.pi*(mjd - 51178)/p3+f3)

       # Applying the bilinear interpolation
       tij = np.array([0,0,0,0], dtype=np.float64)
       for j in range(0, len(x)):
              tij[j] = fun_t(a0[y[j],x[j]], b0[y[j],x[j]], a1[y[j],x[j]], f1[y[j],x[j]], \
                             a2[y[j],x[j]], f2[y[j],x[j]], a3[y[j],x[j]], f3[y[j],x[j]])
       T = np.sum(tij*dxy)/np.sum(dxy)

       # ******************************************************************
       # Open and read the surface pressure coefficients file
       # ******************************************************************
       fid = open(coeffiles+'press_grid.bin', 'rb')
       fid.seek((row*col*20)*hour, 0)
       a0 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       b0 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       a1 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f1 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a2 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f2 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       fid.close()

       # Surface pressure model
       fun_p = lambda a0, b0, a1, f1, a2, f2: a0 + b0*(mjd - 51178) + a1*np.cos(2*np.pi*(mjd - 51178)/p1+f1) + \
       a2*np.cos(2*np.pi*(mjd - 51178)/p2+f2)

       # Applying the bilinear interpolation
       pij = np.array([0,0,0,0], dtype=np.float64)
       for j in range(0, len(x)):
              pij[j] = fun_p(a0[y[j],x[j]], b0[y[j],x[j]], a1[y[j],x[j]], f1[y[j],x[j]], a2[y[j],x[j]], f2[y[j],x[j]])
       P = np.sum(pij*dxy)/np.sum(dxy)

       # ******************************************************************
       # Open and read the surface relative humidity coefficients file
       # ******************************************************************
       fid = open(coeffiles+'rh_grid.bin', 'rb')
       fid.seek((row*col*22)*hour, 0)
       a0 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       a1 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f1 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a2 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f2 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       a3 = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       f3 = np.fromfile(fid, dtype=np.int16, count=row*col).reshape((row, col), order='F')/10000.0
       fid.close()

       # Surface relative humidity model
       fun_rh = lambda a0, a1, f1, a2, f2, a3, f3: a0 + a1*np.cos(2*np.pi*(mjd - 51178)/p1+f1) + \
       a2*np.cos(2*np.pi*(mjd - 51178)/p2+f2) + a3*np.cos(2*np.pi*(mjd - 51178)/p3+f3)

       # Applying the bilinear interpolation
       rhij = np.array([0,0,0,0], dtype=np.float64)
       for j in range(0, len(x)):
              rhij[j] = fun_rh(a0[y[j],x[j]], a1[y[j],x[j]], f1[y[j],x[j]], \
                        a2[y[j],x[j]], f2[y[j],x[j]], a3[y[j],x[j]], f3[y[j],x[j]])
       RH = np.sum(rhij*dxy)/np.sum(dxy)

       # ************************************************************
       # Open and read the Tm coefficients and undulation file
       # ************************************************************
       fid = open(coeffiles+'tm_grid.bin', 'rb')
       a = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       b = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       orography = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       undu = np.fromfile(fid, dtype=np.float32, count=row*col).reshape((row, col), order='F')
       fid.close()

       # Applying the bilinear interpolation
       aij = np.array([0,0,0,0], dtype=np.float64)
       bij = np.copy(aij)
       hij = np.copy(aij)
       nij = np.copy(aij)
       for j in range(0, len(x)):
              aij[j] = a[y[j],x[j]]
              bij[j] = b[y[j],x[j]]
              hij[j] = orography[y[j],x[j]]
              nij[j] = undu[y[j],x[j]]
       a = np.sum(aij*dxy)/np.sum(dxy)
       b = np.sum(bij*dxy)/np.sum(dxy)
       geo_height = np.sum(hij*dxy)/np.sum(dxy)
       N = np.sum(nij*dxy)/np.sum(dxy)

       # Zenith hydrostatic delay (ZHD), Saastamoinen model
       if z0_type=='orth':
              H_orth = z0
       elif z0_type=='elli':
              H_orth = z0 - N
       else:
              raise NameError('Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).')

       # Correction to P, T, and RH (see Guochanf Xu, GPS Theory, Algorithms and Applications, 3nd Edition, page 82)
       dh = H_orth - geo_height

       # Temperature lapse rate by latitude
       # rate = (6.25 - 2*sin( np.deg2rad(y0) )^4)/1000
       rate = (10.3 + 0.03182*(T-273.15) - 0.00436*P)/1000

       # Altimetric corrections
       P = (P*100 * ( 1 - (rate*dh)/T )**5.6004)/100
       T = T - rate * ( dh )

       # Ensure that we eliminate any noise outside this range
       # Gamit has an unidentified error when assuming rh = 100% (met-files)
       if RH >= 100:
              RH = 99.9
       if RH < 0:
              RH = 0

       # Call external function
       es = es_wexler(T, P)
       e = es * (RH/100)

       # Weight mean temperature, Tm
       Tm = a + b*T

       # ZHD using the Saastamoinen Model (see Saastamoinen, 1973)
       ZHD = (0.0022768 * P)/(1 - 0.0026*np.cos(2*np.deg2rad(y0))-0.00000028*H_orth)

       # ZWD using the Saastamoinen Tropospheric Model (see Saastamoinen, 1972, 1973)
       ZWD = 0.002277 * ((1255/T + 0.05) * e)

       # PWV
       # Dimensionless constant
       const = 10**8 / (461525 * (22.9744 + 375463/Tm))
       PWV = ZWD * const

       return P, T, RH, Tm, ZHD, ZWD, PWV