function [P, T, RH, Tm, ZHD, ZWD, PWV] = hgpt2(dt, x0, y0, z0, z0_type)
% hgpt2.m
%
% This routine determines the surface pressure (P, hPa), surface air temperature (T, K), relative humidity (RH, %)
% weighed mean temperature (Tm, K), zenith hydrostatic delay (ZHD, m), zenith wet delay (ZWD, m), and precipitable water vapor (PWV, m)
% from binary coefficient files
% As available from:
% https://github.com/pjmateus/hgpt_model (release v2.0)
% press_grid.bin; temp_grid.bin; tm_grid.bin; and rh_grid.bin
%
% It is admitted that the binary files with the coefficients are in the same directory as this script.
% In alternative you can define the "coeffiles" variable
%
% The epoch can be an array of size 1, and in this case is the Modified Julian Date (MJD) 
% or can be an array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)
%
% All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude
%
% Reference for HGPT and HGPT2:
% Mateus, P.; Catalão, J.; Mendes, V.B.; Nico, G. An ERA5-Based Hourly Global Pressure and Temperature (HGPT) Model. 
% Remote Sens. 2020, 12, 1098. https://doi.org/10.3390/rs12071098
% HGPT2: an ERA5-based global model to estimate relative humidity (Remote Sensing, MDPI)
%
% INPUT:
%              dt : if size(dt)=1 => modified julian date
%                   if size(dt)=6 => year, month, day, hour, min, sec
%              x0 : ellipsoidal longitude (degrees)   
%              y0 : ellipsoidal latitude (degrees)         
%              z0 : height (m)
%         z0_type : ‘orth’ for orthometric height or ‘elli’ for ellipsoidal height
%
% OUTPUT:
%               P : surface pressure valid at (x0, y0, z0), in hPa
%               T : surface air temperature valid at (x0, y0, z0), in Kelvins
%              RH : relative humidity valid at (x0, y0, z0), in %
%              Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
%             ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
%             ZWD : zenith wet delay, valid at (x0, y0, z0), in meters
%             PWV : precipitable water vapor, valid at (x0, y0, z0), in meters
%
%--------------------------------------------------------------------------
% Example:       
%	y0 = 38.5519;
%	x0 = -9.0147;
%	z0 = 25;     
%	dt = 58119.5; or dt = [2018, 1, 1, 12, 0, 0]
%   [P, T, RH, Tm, ZHD, ZWD, PWV] = hgpt2(dt, x0, y0, z0, 'orth')
%   [1020.262, 
%     286.589, 
%      74.908,
%     277.236,
%       2.324,
%       0.117,
%       0.018]
%--------------------------------------------------------------------------
% written by Pedro Mateus (2021/05/15)
% Instituto Dom Luiz (IDL), Faculdade de Ciências, Universidade de Lisboa, 1749-016 Lisboa, Portugal
% pjmateus@fc.ul.pt
%

% Location of coefficient files
coeffiles = '';

% Constants
row = 721;
col = 1440;
p1  = 365.250;
p2  = 182.625;
p3  = 91.3125;

% Geographic coordinates ( equal to ERA5 )
lon = linspace(-180, 179.75, col);
lat = linspace(-90, 90, row); 

% Input datetime format
if length(dt) == 6
    mjd = juliandate( datetime(dt(1), dt(2), dt(3), dt(4), dt(5), dt(6)) ) - 2400000.5;    
    hour = dt(4);
elseif length(dt) == 1
    mjd = dt;
    mjdconv = datevec(mjd+678942);
    hour = mjdconv(4);
else
    error('Use 1) Modified Julian Date (MJD) or 2) Gregorian date (y,m,d,HH,MM,SS).')
end

% Finding indexes for bilinear interpolation
% x-location
[~, indx] = sort( sqrt( (lon-x0).^2 ) );
ix1 = indx( 1 ); ix2 = indx( 2 );
x1 = lon( ix1 ); x2 = lon( ix2 );
x = [ix1, ix1, ix2, ix2];
% y-location
[~, indy] = sort( sqrt( (lat-y0).^2 ) );
jy1 = indy( 1 ); jy2 = indy( 2 );
y1 = lat( jy1 ); y2 = lat( jy2 );
y = [jy1, jy2, jy1, jy2];
% xy-distances (weights)
dx1y1= 1/sqrt( (x1 - x0)^2 + (y1 - y0)^2 );
dx1y2= 1/sqrt( (x1 - x0)^2 + (y2 - y0)^2 );
dx2y1= 1/sqrt( (x2 - x0)^2 + (y1 - y0)^2 );
dx2y2= 1/sqrt( (x2 - x0)^2 + (y2 - y0)^2 );
dxy = [dx1y1, dx1y2, dx2y1, dx2y2];
if sum(isinf(dxy)) > 0
    % Exact point grid
    dxy = [1,0,0,0];    
end

% ******************************************************************
% Open and read the surface air temperature coefficients file
% ******************************************************************
if ~isempty(coeffiles)
    if isunix && ~strcmp(coeffiles(end),'/'), coeffiles(end+1) = '/'; end
    if ~isunix && ~strcmp(coeffiles(end),'\'), coeffiles(end+1) = '\'; end
end
[fid, errmsg] = fopen([coeffiles,'temp_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    fseek(fid, (row*col*26)*hour, -1);
    a0 = fread(fid, [row,col], 'single');
    b0 = fread(fid, [row,col], 'single');
    a1 = fread(fid, [row,col], 'single');
    f1 = (fread(fid, [row,col], 'int16'))./10000;
    a2 = fread(fid, [row,col], 'single');
    f2 = (fread(fid, [row,col], 'int16'))./10000;
    a3 = fread(fid, [row,col], 'single');
    f3 = (fread(fid, [row,col], 'int16'))./10000;
    fclose(fid);
end

% Surface air temperature model
fun_t = @(a0, b0, a1, f1, a2, f2, a3, f3) a0 + b0*(mjd - 51178) + a1*cos(2*pi*(mjd - 51178)/p1+f1) + ...
    a2*cos(2*pi*(mjd - 51178)/p2+f2) + a3*cos(2*pi*(mjd - 51178)/p3+f3);

% Applying the bilinear interpolation
for j = 1 : length(x)
    tij(j) = fun_t(a0(y(j),x(j)), b0(y(j),x(j)), a1(y(j),x(j)), f1(y(j),x(j)), ...
        a2(y(j),x(j)), f2(y(j),x(j)), a3(y(j),x(j)), f3(y(j),x(j)));
end
T = sum(tij.*dxy)/sum(dxy);


% ******************************************************************
% Open and read the surface pressure coefficients file
% ******************************************************************    
[fid, errmsg] = fopen([coeffiles,'press_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    fseek(fid, (row*col*20)*hour, -1);
    a0 = fread(fid, [row,col], 'single');
    b0 = fread(fid, [row,col], 'single');    
    a1 = fread(fid, [row,col], 'single');
    f1 = (fread(fid, [row,col], 'int16'))./10000;
    a2 = fread(fid, [row,col], 'single');
    f2 = (fread(fid, [row,col], 'int16'))./10000;
    fclose(fid);
end

% Surface pressure model
fun_p = @(a0, b0, a1, f1, a2, f2) a0 + b0*(mjd - 51178) + a1*cos(2*pi*(mjd - 51178)/p1+f1) + ...
        a2*cos(2*pi*(mjd - 51178)/p2+f2);

% Applying the bilinear interpolation
for j = 1 : length(x)
    pij(j) = fun_p(a0(y(j),x(j)), b0(y(j),x(j)), a1(y(j),x(j)), f1(y(j),x(j)), a2(y(j),x(j)), f2(y(j),x(j)));
end
P = sum(pij.*dxy)/sum(dxy);


% ******************************************************************
% Open and read the surface relative humidity coefficients file
% ******************************************************************
[fid, errmsg] = fopen([coeffiles,'rh_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    fseek(fid, (row*col*22)*hour, -1);
    a0 = fread(fid, [row,col], 'single');
    a1 = fread(fid, [row,col], 'single');
    f1 = (fread(fid, [row,col], 'int16'))./10000;
    a2 = fread(fid, [row,col], 'single');
    f2 = (fread(fid, [row,col], 'int16'))./10000;
    a3 = fread(fid, [row,col], 'single');
    f3 = (fread(fid, [row,col], 'int16'))./10000;
    fclose(fid);
end

% Surface relative humidity model
fun_rh = @(a0, a1, f1, a2, f2, a3, f3) a0 + a1*cos(2*pi*(mjd - 51178)/p1+f1) + ...
    a2*cos(2*pi*(mjd - 51178)/p2+f2) + a3*cos(2*pi*(mjd - 51178)/p3+f3);

% Applying the bilinear interpolation
for j = 1 : length(x)
    rhij(j) = fun_rh(a0(y(j),x(j)), a1(y(j),x(j)), f1(y(j),x(j)), ...
        a2(y(j),x(j)), f2(y(j),x(j)), a3(y(j),x(j)), f3(y(j),x(j)));
end
RH = sum(rhij.*dxy)/sum(dxy);


% ******************************************************************
% Open and read the Tm coefficients and undulation file
% ******************************************************************
[fid, errmsg] = fopen([coeffiles,'tm_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    a = fread(fid, [row,col], 'single');
    b = fread(fid, [row,col], 'single');
    orography = fread(fid, [row,col], 'single');
    undu = fread(fid, [row,col], 'single');
    fclose(fid);
end

% Applying the bilinear interpolation
for j = 1 : length(x)
    aij(j) = a(y(j),x(j));
    bij(j) = b(y(j),x(j));
    hij(j) = orography(y(j),x(j));
    nij(j) = undu(y(j),x(j));
end
a = sum(aij.*dxy)/sum(dxy);
b = sum(bij.*dxy)/sum(dxy);
geo_height = sum(hij.*dxy)/sum(dxy);
N = sum(nij.*dxy)/sum(dxy);

% Zenith hydrostatic delay (ZHD), Saastamoinen model
if strcmp(z0_type, 'orth')
    H_orth = z0;
elseif strcmp(z0_type, 'elli')
    H_orth = z0 - N;
else
    error('Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).')
end

% Correction to P, T, and RH (see Guochanf Xu, GPS Theory, Algorithms and Applications, 3nd Edition, page 82)
dh = H_orth - geo_height;

% Temperature lapse rate by latitude 
% rate = (6.25 - 2*sin( deg2rad(y0) )^4)/1000;
rate = (10.3 + 0.03182*(T-273.15) - 0.00436*P)/1000;

% Altimetric corrections
P = (P*100 * ( 1 - (rate*dh)/T )^5.6004)/100;
T = T - rate * ( dh ); 

% Ensure that we eliminate any noise outside this range
% Gamit has an unidentified error when assuming rh = 100% (met-files)
if RH >= 100, RH = 99.9; end
if RH < 0, RH = 0.0; end

% Call external function 
es = es_wexler(T, P);
e = es * (RH/100);

% Weight mean temperature, Tm
Tm = a + b*T;

% ZHD using the Saastamoinen Model (see Saastamoinen, 1973)
ZHD = (0.0022768 * P)/(1 - 0.0026*cos(2*deg2rad(y0))-0.00000028*H_orth);

% ZWD using the Saastamoinen Tropospheric Model (see Saastamoinen, 1972, 1973)
ZWD = 0.002277 * ((1255/T + 0.05) * e);

% PWV 
% Dimensionless constant
const = 10^8 / (461525 * (22.9744 + 375463/Tm));
PWV = ZWD * const;

return