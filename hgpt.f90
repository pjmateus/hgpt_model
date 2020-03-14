! hgpt.f90
!
! This routine determines the surface pressure (P), surface air temperature (T), 
! weighed mean temperature (Tm), and zenith hydrostatic delay (ZHD) from binary coefficient files
! As available from:
! https://github.com/pjmateus/hgpt_model (release v1.0)
! press_grid.bin; temp_grid.bin; and tm_grid.bin
!
! It is admitted that the binary files with the coefficients are in the same directory as this script
! If you don't want it you can change it in the open statement
!
! The epoch can be an real*8 array of size 1, and in this case is the Modified Julian Date (MJD) 
! or can be an integer array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)
!  
! All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude
!
! Reference for HGPT:
! An ERA5-based hourly global temperature and pressure (HGTP) model (submitted to Remote Sensing, MDPI)
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
!              Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
!             ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
!
! Set the following variables in the script from which you call hgpt.f90:
!--------------------------------------------------------------------------!
!       real :: x0, y0, z0, P, T, Tm, ZHD                                  ! 
!	real*8, dimension(6) :: dt           ! Gregorian Calendar          !
!                                                                          !
! Define variables values and call hgtp:                                   !
! 	call hgpt(dt, 6, x0, y0, z0, 'orth', P, T, Tm, ZHD)                !
! or                                                                       !
!       real :: x0, y0, z0, P, T, Tm, ZHD                                  ! 
!	real*8, dimension(1) :: dt           ! Modified Julian Date (MJD)  !
!                                                                          !
! Define variables values and call hgtp in the same way:                   !
! 	call hgpt(dt, 1, x0, y0, z0, 'orth', P, T, Tm, ZHD)                !
!--------------------------------------------------------------------------!
!Example:                                                                  !
!program call_hgpt                                                         !
!	implicit None                                                      !
!       real :: x0, y0, z0, P, T, Tm, ZHD                                  !
!	real*8, dimension(1) :: dt                                         !
!	y0 = 38.5519                                                       !
!	x0 = -9.0147                                                       !
!	z0 = 25                                                            !
!	dt(1) = 58119.5                                                    !
!	call hgpt(dt, size(dt), x0, y0, z0, 'orth', P, T, Tm, ZHD)         !
!	write(*,*) P, T, Tm, ZHD                                           !
!end program call_hgpt                                                     !
!--------------------------------------------------------------------------!
! written by Pedro Mateus (2020/01/15)
! Instituto Dom Luiz (IDL), Faculdade de Ciências, Universidade de Lisboa, 1749-016 Lisboa, Portugal
! pjmateus@fc.ul.pt
! 
!
subroutine hgpt(dt, ndt, x0, y0, z0, z0_type, P, T, Tm, ZHD)
	implicit None
	real, intent(in) :: x0, y0, z0
	integer, intent(in) :: ndt
	real*8, intent(in) :: dt(ndt)
	character(len=4), intent(in) :: z0_type
	real, intent(out) :: T, P, Tm, ZHD
	real*8  :: mjd
	integer :: ios, i, j, year, month, day, hour, minute, sec 
	integer, parameter :: row = 721, col = 1440     
	real, parameter :: pi = 4*atan(1.0), p1 = 365.25, p2 = 182.6250, p3 = 91.3125, deg2rad = pi/180.0 
	real, dimension(:) :: lon(col), lat(row)
	real, dimension(:,:) :: y_intercept(row, col), slope(row, col), a1(row, col), a2(row, col), a3(row, col), & 
	                        orography(row, col), undu(row, col)
	integer*2, dimension(:,:) :: f1(row, col), f2(row, col), f3(row, col) 
	real :: a, b, amp1, pha1, amp2, pha2, amp3, pha3, alt, N, H_orth, H_ellip

	! Input datetime format
	if ( ndt == 6)  then
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
	lon = (/((  -180+0.25*i), i=1, col)/)
	lat = (/((-90.25+0.25*i), i=1, row)/)

	! Locate i, j indices 
	call locate_ij(row, lat, y0, i)
	call locate_ij(col, lon, x0, j)

	! Open and read the surface air temperature coefficients file
	open(10, file='temp_grid.bin', status='old', action='read', &
	     form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 )	stop 'Error opening file temp_grid.bin'		 

	! hour varies from 0 to 23 UTC
	call fseek(10, (row*col*26)*hour, 0)
	read(10) y_intercept	
	read(10) slope			
	read(10) a1
	read(10) f1
	read(10) a2
	read(10) f2
	read(10) a3
	read(10) f3
	close(10)

	! Bilinear interpolation
	call interpolate(lon, lat, y_intercept, x0, y0, i, j, a)
	call interpolate(lon, lat, slope, x0, y0, i, j, b)	
	call interpolate(lon, lat, a1, x0, y0, i, j, amp1)
	call interpolate(lon, lat, float(f1)/10000.0, x0, y0, i, j, pha1)
	call interpolate(lon, lat, a2, x0, y0, i, j, amp2)
	call interpolate(lon, lat, float(f2)/10000.0, x0, y0, i, j, pha2)
	call interpolate(lon, lat, a3, x0, y0, i, j, amp3)
	call interpolate(lon, lat, float(f3)/10000.0, x0, y0, i, j, pha3)

	! Surface air temperature model
	T = a + b*(mjd - 51178) + amp1*cos(2*pi*(mjd - 51178)/p1+pha1) + & 
             		          amp2*cos(2*pi*(mjd - 51178)/p2+pha2) + & 
             	   	          amp3*cos(2*pi*(mjd - 51178)/p3+pha3)
	
	! Open and read the surface pressure coefficients file
	open(11, file='press_grid.bin', status='old', action='read', &
	     form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 ) stop 'Error opening file press_grid.bin'
        		
	! hour varies from 0 to 23 UTC
	call fseek(11, (row*col*20)*hour, 0)
	read(11) y_intercept	
	read(11) slope			
	read(11) a1
	read(11) f1
	read(11) a2
	read(11) f2
	close(11)	

	! Bilinear interpolation
	call interpolate(lon, lat, y_intercept, x0, y0, i, j, a)
	call interpolate(lon, lat, slope, x0, y0, i, j, b)	
	call interpolate(lon, lat, a1, x0, y0, i, j, amp1)
	call interpolate(lon, lat, float(f1)/10000.0, x0, y0, i, j, pha1)
	call interpolate(lon, lat, a2, x0, y0, i, j, amp2)
	call interpolate(lon, lat, float(f2)/10000.0, x0, y0, i, j, pha2)

	! Surface pressure model
	P = a + b*(mjd - 51178) + amp1*cos(2*pi*(mjd - 51178)/p1+pha1) + & 
                	          amp2*cos(2*pi*(mjd - 51178)/p2+pha2) 

	! Open and read the weight mean temperature and undulation coefficients file
	open(12, file='tm_grid.bin', status='old', action='read', &
             form='unformatted', access='stream', iostat=ios)
	if ( ios /= 0 ) stop 'Error opening file tm_grid.bin'		 
        
	read(12) y_intercept	
	read(12) slope			
	read(12) orography
	read(12) undu
	close(12)	

	! Bilinear interpolation
	call interpolate(lon, lat, y_intercept, x0, y0, i, j, a)
	call interpolate(lon, lat, slope, x0, y0, i, j, b)	
	call interpolate(lon, lat, orography, x0, y0, i, j, alt)
	call interpolate(lon, lat, undu, x0, y0, i, j, N)

	! Zenith hydrostatic delay (ZHD), Saastamoinen model
	if ( z0_type == 'orth' ) then		
		H_orth = z0
	else if ( z0_type == 'elli' ) then 
		H_orth = z0 - N		
	else
		stop 'Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).'
	end if

	! Correction to P and T (see Guochanf Xu, GPS Theory, Algorithms and Applications, 2nd Edition, page 56)
	P = P * (1.0 - 0.000226*(H_orth - alt))**5.225
	T = T - 0.0065*(H_orth - alt)

	! Weight mean temperature (Tm) linear model 
	Tm = a + b * T

	! ZHD using the Saastamoinen Model (see Saastamoinen, 1973) 
	ZHD = (0.0022768 * P)/(1 - 0.0026*cos(2*deg2rad*y0-0.00000028)*H_orth)

	return
end subroutine hgpt 

subroutine locate_ij(length, array, value, idx)
        ! Given an array and a value, returns the index of the element that
        ! is closest to, but less than, the given value.
        ! Uses a binary search algorithm.
        implicit none
        integer, intent(in) :: length
        real, dimension(length), intent(in) :: array      
        real, intent(in) :: value
	integer, intent(out) :: idx
        integer :: left, middle, right
       
        left = 1
        right = length
        do
            if (left > right) then
                exit
            end if
            middle = nint((left+right) / 2.0)
            if ( abs(array(middle) - value) <= 1e-9) then
                idx = middle
                return
            else if (array(middle) > value) then
                right = middle - 1
            else
                left = middle + 1
            end if
        end do
        idx = right

	return
end subroutine locate_ij

subroutine interpolate(lon, lat, f, x0, y0, i, j, value)
        ! This function uses bilinear interpolation to estimate the value
        ! of a function f at point (x0,y0) using the point grid (i, j), (i+1,j), (i,j+1) 
        ! and (i+1,j+1), see locate_ij subroutine
	! In the limits (borders) the nearest-neighbor interpolation is used
        ! f is assumed to be sampled on a regular grid, with the grid x values specified
        ! by lon array and the grid y values specified by lat array
        implicit none   
	integer, parameter :: col = 1440, row = 721     
        real, dimension(:), intent(in) :: lon(col), lat(row)	
        real, dimension(row, col), intent(in) :: f
        integer, intent(in) :: i, j
	real, intent(in) :: x0, y0
	integer :: ix, jy
	real, intent(out) :: value
        real :: x1, x2, y1, y2, xy1, xy2

	ix = i
	jy = j
	if ( ix == 0 .or. ix >= row .or. jy == 0 .or. jy >= col ) then
		! Nearest-neighbor interpolation
		if ( ix == 0 ) ix = ix + 1
		if ( ix >= row ) ix = row
		if ( jy == 0 ) jy = jy + 1
		if ( jy >= col ) jy = col
		value = f(ix, jy)
	else
		! Bilinear interpolation
		x1 = lon( jy )
		x2 = lon( jy + 1 )
		y1 = lat( ix )
		y2 = lat( ix + 1 )		
		xy1 = (x2 - x0)/(x2 - x1)*f(ix, jy) + (x0 - x1)/(x2 - x1)*f(ix, jy+1)
		xy2 = (x2 - x0)/(x2 - x1)*f(ix+1, jy) + (x0 - x1)/(x2 - x1)*f(ix+1, jy+1)
                value = (y2 - y0)/(y2 - y1)*xy1 + (y0 - y1)/(y2 - y1)*xy2
	end if  

	return
end subroutine interpolate

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

