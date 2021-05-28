! hgpt2.f90
!
! This routine determines the surface pressure (P, hPa), surface air temperature (T, K), relative humidity (RH, %)
! weighed mean temperature (Tm, K), zenith hydrostatic delay (ZHD, m), zenith wet delay (ZWD, m), and
! precipitable water vapor (PWV, m) from binary coefficient files
! As available from:
! https://github.com/pjmateus/hgpt_model (release v2.0)
! press_grid.bin; temp_grid.bin; tm_grid.bin; and rh_grid.bin 
!
! It is admitted that the binary files with the coefficients are in the same directory as this script
! If you don't want it you can change it in the open statement
!
! The epoch can be an real*8 array of size 1, and in this case is the Modified Julian Date (MJD) 
! or can be an integer array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)
!  
! All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude
!
! Reference for HGPT and HGPT2:
! Mateus, P.; Catalão, J.; Mendes, V.B.; Nico, G. An ERA5-Based Hourly Global Pressure and Temperature (HGPT) Model. 
! Remote Sens. 2020, 12, 1098. https://doi.org/10.3390/rs12071098
! HGPT2: an ERA5-based global model to estimate relative humidity (Remote Sensing, MDPI)
!
! INPUT:
!              dt : if size(dt)=1 => modified julian date
!                   if size(dt)=6 => year, month, day, hour, min, sec
!             ndt : size of dt
!              x0 : ellipsoidal longitude (degrees)   
!              y0 : ellipsoidal latitude (degrees)         
!              z0 : height (m)
!         z0_type : ‘orth’ for orthometric height or ‘elli’ for ellipsoidal height
!
! OUTPUT:
!               P : surface pressure valid at (x0, y0, z0), in hPa
!               T : surface air temperature valid at (x0, y0, z0), in Kelvins
!              RH : relative humidity valid at (x0, y0, z0), in %
!              Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
!             ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
!             ZWD : zenith wet delay, valid at (x0, y0, z0), in meters
!             PWV : precipitable water vapor, valid at (x0, y0, z0), in meters
!
! Set the following variables in the script from which you call hgpt2.f90:
!--------------------------------------------------------------------------!
!       real :: x0, y0, z0, P, T, RH, Tm, ZHD, ZWD, PWV                    !
!	real*8, dimension(6) :: dt           ! Gregorian Calendar          !
!                                                                          !
! Define variables values and call hgtp:                                   !
! 	call hgpt2(dt, 6, x0, y0, z0, 'orth', P, T, RH, Tm, ZHD, ZWD, PWV) !
! or                                                                       !
!       real :: x0, y0, z0, P, T, RH, Tm, ZHD, ZWD, PWV                    ! 
!	real*8, dimension(1) :: dt           ! Modified Julian Date (MJD)  !
!                                                                          !
! Define variables values and call hgtp in the same way:                   !
! 	call hgpt2(dt, 1, x0, y0, z0, 'orth', P, T, RH, Tm, ZHD, ZWD, PWV) !
!--------------------------------------------------------------------------!
!Example:                                                                  !
!program call_hgpt2                                                        !
!	implicit None                                                      !
!       real :: x0, y0, z0, P, T, RH, Tm, ZHD, ZWD, PWV                    !
!	real*8, dimension(1) :: dt                                         !
!	y0 = 38.5519                                                       !
!	x0 = -9.0147                                                       !
!	z0 = 25                                                            !
!	dt(1) = 58119.5                                                    !
!	call hgpt2(dt,size(dt),x0,y0,z0,'orth',P,T,RH,Tm,ZHD,ZWD,PWV)      !
!	write(*,*) P, T, RH, Tm, ZHD, ZWD, PWV                             !
!end program call_hgpt2                                                    !
!--------------------------------------------------------------------------!
! written by Pedro Mateus (2020/01/15)
! Instituto Dom Luiz (IDL), Faculdade de Ciências, Universidade de Lisboa, 1749-016 Lisboa, Portugal
! pjmateus@fc.ul.pt
! 
!
subroutine hgpt2(dt, ndt, x0, y0, z0, z0_type, P, T, RH, Tm, ZHD, ZWD, PWV)
	implicit None
	real, intent(in) :: x0, y0, z0
	integer, intent(in) :: ndt
	real*8, intent(in) :: dt(ndt)
	character(len=4), intent(in) :: z0_type
	real, intent(out) :: P, T, RH, Tm, ZHD, ZWD, PWV
	real*8  :: mjd
	integer :: ios, i, j, year, month, day, hour, minute, sec 
	integer, parameter :: row = 721, col = 1440     
	real, parameter :: pi = 4*atan(1.0), deg2rad = pi/180.0 
	real, dimension(:)   :: lon(col), lat(row), indx(col), indy(row)
	real, dimension(:,:) :: a0(row, col), b0(row, col), a1(row, col), a2(row, col), a3(row, col), & 
	                        orography(row, col), undu(row, col)
	integer*2, dimension(:,:) :: f1(row, col), f2(row, col), f3(row, col) 
	real, dimension(:,:)      :: pha1(row, col), pha2(row, col), pha3(row, col)
	real :: a, b, geo_height, N, H_orth, H_ellip, x1, x2, y1, y2, es_wexler, const, &
                dx1y1, dx1y2, dx2y1, dx2y2, dxy(4), fun_t, fun_p, fun_rh, tij(4), pij(4), & 
                rhij(4), aij(4), bij(4), hij(4), nij(4), dh, rate, es, e
        integer :: ix1(1), ix2(1), x(4), iy1(1), iy2(1), y(4)

	! Input datetime format
	if ( ndt == 6 )  then
		year = idint(dt(1)); month = idint(dt(2)); day = idint(dt(3))
		hour = idint(dt(4)); minute = idint(dt(5)); sec = idint(dt(6))		
		call greg2modjul(year, month, day, hour, minute, sec, mjd)
	else if ( ndt == 1 ) then
		mjd = dt(1)
		call modjul2greg(mjd, year, month, day, hour, minute, sec)
	else
		stop 'Use 1) Modified Julian Date (MJD) or 2) Gregorian date (y,m,d,HH,MM,SS).'
	end if
	
	! Geographic latitude and longitude (equal to ERA5 grid)  
	lon = (/((-180.25+0.25*i), i=1, col)/)
	lat = (/(( -90+0.25*i), i=0, row-1)/)	

	! Finding indexes for bilinear interpolation
	! x-location
        indx = sqrt( (lon-x0)**2 ) 
	ix1 = minloc( indx ); ix2 = minloc( indx, MASK = (indx > indx(ix1(1))) )	
	x1 = lon( ix1(1) ); x2 = lon( ix2(1) )
	x = (/ ix1(1), ix1(1), ix2(1), ix2(1) /)	
	! y-location
	indy = sqrt( (lat-y0)**2 ) 
	iy1 = minloc( indy ); iy2 = minloc( indy, MASK = (indy > indy(iy1(1))) )
	y1 = lat( iy1(1) ); y2 = lat( iy2(1) )
	y = (/ iy1(1), iy2(1), iy1(1), iy2(1) /)
	! xy-distances (weights)
        dx1y1 = sqrt( (x1 - x0)**2 + (y1 - y0)**2 )	
	dx1y2 = sqrt( (x1 - x0)**2 + (y2 - y0)**2 )
	dx2y1 = sqrt( (x2 - x0)**2 + (y1 - y0)**2 )
	dx2y2 = sqrt( (x2 - x0)**2 + (y2 - y0)**2 )
	if ( dx1y1 == 0 .or. dx1y2 == 0 .or. dx2y1 == 0 .or. dx2y2 == 0 ) then    
		dxy = (/ 1, 0, 0, 0 /)
        else  
		dxy = (/ 1/dx1y1, 1/dx1y2, 1/dx2y1, 1/dx2y2 /)
	end if

        ! ***********************************************************************
	! Open and read the surface air temperature coefficients file
        ! ***********************************************************************
	open(10, file='temp_grid.bin', status='old', action='read', &
	     form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 )	stop 'Error opening file temp_grid.bin'		 

	! hour varies from 0 to 23 UTC
	call fseek(10, (row*col*26)*hour, 0)
	read(10) a0	
	read(10) b0			
	read(10) a1
	read(10) f1; pha1 = float(f1)/10000.0
	read(10) a2
	read(10) f2; pha2 = float(f2)/10000.0
	read(10) a3
	read(10) f3; pha3 = float(f3)/10000.0
	close(10)			
	! Applying the bilinear interpolation
	do j = 1, 4
    		tij(j) = fun_t( a0(y(j),x(j)), b0(y(j),x(j)), a1(y(j),x(j)), pha1(y(j),x(j)), &
        		        a2(y(j),x(j)), pha2(y(j),x(j)), a3(y(j),x(j)), pha3(y(j),x(j)), &
                                mjd)
	end do
	T = sum(tij*dxy)/sum(dxy)

        ! ***********************************************************************
	! Open and read the surface pressure coefficients file
        ! ***********************************************************************
	open(11, file='press_grid.bin', status='old', action='read', &
	     form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 ) stop 'Error opening file press_grid.bin'
        		
	! hour varies from 0 to 23 UTC
	call fseek(11, (row*col*20)*hour, 0)
	read(11) a0	
	read(11) b0			
	read(11) a1
	read(11) f1; pha1 = float(f1)/10000.0
	read(11) a2
	read(11) f2; pha2 = float(f2)/10000.0
	close(11)	
	! Applying the bilinear interpolation
	do j = 1, 4
    		pij(j) = fun_p(a0(y(j),x(j)), b0(y(j),x(j)), a1(y(j),x(j)), pha1(y(j),x(j)), &
                               a2(y(j),x(j)), pha2(y(j),x(j)), mjd)
	end do
	P = sum(pij*dxy)/sum(dxy)
	
	! ***********************************************************************
	! Open and read the surface relative humidity coefficients file
	! ***********************************************************************
	open(12, file='rh_grid.bin', status='old', action='read', &
	     form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 ) stop 'Error opening file rh_grid.bin'

	! hour varies from 0 to 23 UTC
	call fseek(12, (row*col*22)*hour, 0)
	read(12) a0				
	read(12) a1
	read(12) f1; pha1 = float(f1)/10000.0
	read(12) a2
	read(12) f2; pha2 = float(f2)/10000.0
	read(12) a3
	read(12) f3; pha3 = float(f3)/10000.0
	close(12)
	! Applying the bilinear interpolation
	do j = 1, 4
    		rhij(j) = fun_rh(a0(y(j),x(j)), a1(y(j),x(j)), pha1(y(j),x(j)), &
                                 a2(y(j),x(j)), pha2(y(j),x(j)), a3(y(j),x(j)), pha3(y(j),x(j)), &
			         mjd)
	end do
	RH = sum(rhij*dxy)/sum(dxy)

        ! ***********************************************************************
	! Open and read the weight mean temperature and undulation coefficients file
        ! ***********************************************************************
	open(13, file='tm_grid.bin', status='old', action='read', &
             form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 ) stop 'Error opening file tm_grid.bin'		 
        
	read(13) a0	
	read(13) b0			
	read(13) orography
	read(13) undu
	close(13)	
	! Applying the bilinear interpolation
	do j = 1, 4
	    aij(j) = a0(y(j),x(j))
	    bij(j) = b0(y(j),x(j))
	    hij(j) = orography(y(j),x(j))
	    nij(j) = undu(y(j),x(j))
	end do
	a = sum(aij*dxy)/sum(dxy)
	b = sum(bij*dxy)/sum(dxy)
	geo_height = sum(hij*dxy)/sum(dxy)
	N = sum(nij*dxy)/sum(dxy)

	! Zenith hydrostatic delay (ZHD), Saastamoinen model
	if ( z0_type == 'orth' ) then		
		H_orth = z0
	else if ( z0_type == 'elli' ) then 
		H_orth = z0 - N		
	else
		stop 'Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).'
	end if
	
	! Correction to P, T, and RH (see Guochanf Xu, GPS Theory, Algorithms and Applications, 3nd Edition, page 82)
	dh = H_orth - geo_height

	! Temperature lapse rate by latitude 
	! rate = (6.25 - 2*sin( deg2rad(y0) )^4)/1000;
	rate = (10.3 + 0.03182*(T-273.15) - 0.00436*P)/1000

	! Altimetric corrections
	P = (P*100 * ( 1 - (rate*dh)/T )**5.6004)/100
	T = T - rate * ( dh )

	! Ensure that we eliminate any noise outside this range
	! Gamit has an unidentified error when assuming rh = 100% (met-files)
	if ( RH >= 100 ) RH = 99.9 
	if ( RH < 0 ) RH = 0.0 

	! Call external function 
	es = es_wexler(T, P)
	e = es * (RH/100)

	! Weight mean temperature, Tm
	Tm = a + b*T

	! ZHD using the Saastamoinen Model (see Saastamoinen, 1973)
	ZHD = (0.0022768 * P)/(1 - 0.0026*cos(2*deg2rad*y0)-0.00000028*H_orth)

	! ZWD using the Saastamoinen Tropospheric Model (see Saastamoinen, 1972, 1973)
	ZWD = 0.002277 * ((1255/T + 0.05) * e)

	! PWV 
	! Dimensionless constant
	const = 10**8 / (461525 * (22.9744 + 375463/Tm))
	PWV = ZWD * const

	return
end subroutine hgpt2 

function fun_t(a0, b0, a1, f1, a2, f2, a3, f3, mjd)
	! Surface air temperature model
	real, intent(in)   :: a0, b0, a1, f1, a2, f2, a3, f3
	real*8, intent(in) :: mjd
  	real               :: fun_t
        real, parameter    :: pi = 4*atan(1.0), p1 = 365.250, p2 = 182.6250, p3 = 91.3125
	fun_t = a0 + b0*(mjd - 51178) + a1*cos(2*pi*(mjd - 51178)/p1+f1) + &
    		a2*cos(2*pi*(mjd - 51178)/p2+f2) + a3*cos(2*pi*(mjd - 51178)/p3+f3)
	return
end function

function fun_p(a0, b0, a1, f1, a2, f2, mjd)
	! Surface pressure model
	real, intent(in)   :: a0, b0, a1, f1, a2, f2
	real*8, intent(in) :: mjd
  	real               :: fun_p
        real, parameter    :: pi = 4*atan(1.0), p1 = 365.250, p2 = 182.6250
	fun_p = a0 + b0*(mjd - 51178) + a1*cos(2*pi*(mjd - 51178)/p1+f1) + &
                a2*cos(2*pi*(mjd - 51178)/p2+f2)
	return
end function

function fun_rh(a0, a1, f1, a2, f2, a3, f3, mjd)
	! Surface relative humidity model
	real, intent(in)   :: a0, a1, f1, a2, f2, a3, f3
	real*8, intent(in) :: mjd
  	real               :: fun_rh
        real, parameter    :: pi = 4*atan(1.0), p1 = 365.250, p2 = 182.6250, p3 = 91.3125
	fun_rh = a0 + a1*cos(2*pi*(mjd - 51178)/p1+f1) + &
                 a2*cos(2*pi*(mjd - 51178)/p2+f2) + a3*cos(2*pi*(mjd - 51178)/p3+f3)
	return
end function

function es_wexler(T, P)
        !Saturation Water Vapor Pressure (es) using Wexler formulation with new coefficients (adjusted for ITS-90)
        !INPUT : T = Temperature (in Kelvins)
        !        P = atmospheric Pressure (in hPa)
        !OUTPUT: es= saturation water vapor pressure (in hPa)
        !References:
        !1)
        !Wexler, A. Vapor Pressure Formulation for Water in Range 0 to 100 Degrees C. A Revision.
        !J. Res. Natl. Bur. Stand. 1976, 80A, 775–785.
        !2)
        !Wexler, A. Vapor Pressure Formulation for Ice. J. Res. Natl. Bur. Stand. 1977, 81A, 5–20.
        !3)
        !Hardy, B. ITS-90 Formulations for Water Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in Range -100 to +100 C.
        !In Proceedings of the Third International Symposium on Humidity and Moisture;
        !UK National Physical Laboratory (NPL): Teddington, UK, April 6 1998; pp. 1–8.
	real, intent(in) :: T, P
  	real             :: es_wexler   
        real :: g0, g1, g2, g3, g4, g5, g6, g7, A0, A1, A2, A3, B0, B1, B2, B3, &
                k0, k1, k2, k3, k4, k5, es, alpha, beta, f                

        if ( T >= 273.15 ) then
		! Saturation Vapor Pressure over Water
                g0 =-2.8365744d3
                g1 =-6.028076559d3
                g2 = 1.954263612d1
                g3 =-2.737830188d-2
                g4 = 1.6261698d-5
                g5 = 7.0229056d-10
                g6 =-1.8680009d-13
                g7 = 2.7150305
                es = 0.01 * exp(g0*T**(-2) + g1*T**(-1) + g2 + g3*T + g4*T**2 + g5*T**3 + g6*T**4 + g7*log(T))

                ! Enhancement Factor coefficients for Water 0 to 100°C
                A0 =-1.6302041d-1
                A1 = 1.8071570d-3
                A2 =-6.7703064d-6
                A3 = 8.5813609d-9
                B0 =-5.9890467d1
                B1 = 3.4378043d-1
                B2 =-7.7326396d-4
                B3 = 6.3405286d-7
	else
                ! Saturation Vapor Pressure over Ice
                k0 =-5.8666426d3
                k1 = 2.232870244d1
                k2 = 1.39387003d-2
                k3 =-3.4262402d-5
                k4 = 2.7040955d-8
                k5 = 6.7063522d-1
                es = 0.01 * exp(k0*T**(-1) + k1 + k2*T + k3*T**2 + k4*T**3 + k5*log(T))

                ! Enhancement Factor coefficients for Ice –100 to 0°C
                A0 =-6.0190570d-2
                A1 = 7.3984060d-4
                A2 =-3.0897838d-6
                A3 = 4.3669918d-9
                B0 =-9.4868712d1
                B1 = 7.2392075d-1
                B2 =-2.1963437d-3
                B3 = 2.4668279d-6
	end if 
        ! Enhancement Factor
        alpha = A0 + A1*T + A2*T**2 + A3*T**3
        beta = exp(B0 + B1*T + B2*T**2 + B3*T**3)
        f = exp( alpha*(1-es/P) + beta*(P/es-1) )
        es = es * f
	return
end function

subroutine greg2modjul(year, month, day, hour, minute, sec, mjd) 
	! greg2modjul(..., mjd) returns the
	! Modified Julian Date (MJD) corresponding to the 
	! Gregorian calendar date (year, month, day, hour, minute, and second).
	implicit none
        integer, intent(in) :: year, month, day, hour, minute, sec
	real*8, intent(out) :: mjd 
        real*8              :: yyyy, mm, jd   

        ! January & February
	if ( month <= 2 ) then 
        	yyyy = dfloat(year) - 1_8		
        	mm = dfloat(month) + 12_8
	else 
        	yyyy = dfloat(year)		
        	mm = dfloat(month)
        end if

	jd = floor( 365.25_8*(yyyy + 4716.0_8)) + floor( 30.6001_8*( mm + 1.0_8)) + 2.0_8 - &
	     floor( yyyy/100.0_8 ) + floor( floor( yyyy/100.0_8 )/4.0_8 ) + dfloat(day) - 1524.5_8 + &
             (dfloat(hour) + dfloat(minute)/60_8 + dfloat(sec)/3600_8)/24_8
        mjd = jd - 2400000.5_8

	return
end subroutine greg2modjul

subroutine modjul2greg(mjd, year, month, day, hour, minute, sec) 
	! modjul2greg(mjd, ...) returns the
	! Gregorian calendar date (year, month, day, hour, minute, and second)
	! corresponding to the Modified Julian Date (MJD).
	implicit none
        real*8, intent(in)  :: mjd
	real*8              :: jd, ijd, fjd, a, b, c, d, e, m
	integer, intent(out):: year, month, day, hour, minute, sec            

	jd = mjd + 2400000.5_8
	ijd = floor(jd + 0.5_8)
        fjd = jd - ijd + 0.5_8

        sec = 86400_8 * fjd
        hour = ifix( sec / 3600.0 )
        sec = sec - 3600_8 * hour
        minute = ifix( sec / 60.0 )
        sec = sec - 60_8 * minute

        a = ijd + 32044_8
        b = floor((4_8 * a + 3_8) / 146097_8)
        c = a - floor((b * 146097_8) / 4_8)
	d = floor((4_8 * c + 3_8) / 1461_8)
        e = c - floor((1461_8 * d) / 4_8)
        m = floor((5_8 * e + 2_8) / 153_8)

        day   = e - floor((153_8 * m + 2_8) / 5_8) + 1_8
        month = m + 3_8 - 12_8 * floor(m / 10_8)
        year  = b * 100_8 + d - 4800_8 + floor(m / 10_8)

	return
end subroutine modjul2greg

