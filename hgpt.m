function [P, T, Tm, ZHD] = hgpt(dt, x0, y0, z0, z0_type)
% hgpt.m
%
% This routine determines the surface pressure (P), surface air temperature (T), 
% weighed mean temperature (Tm), and zenith hydrostatic delay (ZHD) from binary coefficient files
% As available from:
% https://github.com/pjmateus/hgpt_model (release v1.0)
% press_grid.bin; temp_grid.bin; and tm_grid.bin
%
% It is admitted that the binary files with the coefficients are in the same directory as this script.
% In alternative you can define the "coeffiles" variable
%
% The epoch can be an array of size 1, and in this case is the Modified Julian Date (MJD) 
% or can be an array of size 6, with the Gregorian Calendar in the following format (year, month, day, hour, min, sec)
%
% All parameters are bilinear interpolated to the input ellipsoidal longitude and latitude
%
% Reference for HGPT:
% An ERA5-based hourly global temperature and pressure (HGTP) model (submitted to Remote Sensing, MDPI)
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
%              Tm : weighed mean temperature valid at (x0, y0, z0), in Kelvins
%             ZHD : zenith hydrostatic delay, valid at (x0, y0, z0), in meters
%
%--------------------------------------------------------------------------
% Example:       
%	y0 = 38.5519;
%	x0 = -9.0147;
%	z0 = 25;     
%	dt = 58119.5;
%   [P, T, Tm, ZHD] = hgpt(dt, x0, y0, z0, 'orth')
%--------------------------------------------------------------------------
% written by Pedro Mateus (2020/01/15)
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
lon = linspace(-179.75, 180, col);
lat = linspace(-90, 90, row);

% Input datetime format
if length(dt) == 6
    mjd = mjuliandate(dt(1), dt(2), dt(3), dt(4), dt(5), dt(6));
    hour = dt(4);
elseif length(dt) == 1
    mjd = dt;
    mjdconv = datevec(mjd+678942);
    hour = mjdconv(4);
else
    error('Use 1) Modified Julian Date (MJD) or 2) Gregorian date (y,m,d,HH,MM,SS).')
end

% Open and read the surface air temperature coefficients file
if ~isempty(coeffiles)
    if isunix && ~strcmp(coeffiles(end),'/'), coeffiles(end+1) = '/'; end
    if ~isunix && ~strcmp(coeffiles(end),'\'), coeffiles(end+1) = '\'; end
end
[fid, errmsg] = fopen([coeffiles,'temp_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    fseek(fid, (row*col*26)*hour, -1);
    y_intercept = fread(fid, [row,col], 'single');
    slope = fread(fid, [row,col], 'single');
    a1 = fread(fid, [row,col], 'single');
    f1 = (fread(fid, [row,col], 'int16'))./10000;
    a2 = fread(fid, [row,col], 'single');
    f2 = (fread(fid, [row,col], 'int16'))./10000;
    a3 = fread(fid, [row,col], 'single');
    f3 = (fread(fid, [row,col], 'int16'))./10000;
    fclose(fid);
end

% Bilinear interpolation
F = griddedInterpolant({lon,lat}, y_intercept', 'linear', 'linear'); a = F(x0, y0);
F = griddedInterpolant({lon,lat}, slope', 'linear', 'linear'); b = F(x0, y0);
F = griddedInterpolant({lon,lat}, a1', 'linear', 'linear'); amp1 = F(x0, y0);
F = griddedInterpolant({lon,lat}, f1', 'linear', 'linear'); pha1 = F(x0, y0);
F = griddedInterpolant({lon,lat}, a2', 'linear', 'linear'); amp2 = F(x0, y0);
F = griddedInterpolant({lon,lat}, f2', 'linear', 'linear'); pha2 = F(x0, y0);
F = griddedInterpolant({lon,lat}, a3', 'linear', 'linear'); amp3 = F(x0, y0);
F = griddedInterpolant({lon,lat}, f3', 'linear', 'linear'); pha3 = F(x0, y0);

% Surface air temperature model
T = a + b*(mjd - 51178) + amp1*cos(2*pi*(mjd - 51178)/p1+pha1) + ...
    amp2*cos(2*pi*(mjd - 51178)/p2+pha2) + ...
    amp3*cos(2*pi*(mjd - 51178)/p3+pha3);

% Open and read the surface pressure coefficients file
[fid, errmsg] = fopen([coeffiles,'press_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    fseek(fid, (row*col*20)*hour, -1);
    y_intercept = fread(fid, [row,col], 'single');
    slope = fread(fid, [row,col], 'single');
    a1 = fread(fid, [row,col], 'single');
    f1 = (fread(fid, [row,col], 'int16'))./10000;
    a2 = fread(fid, [row,col], 'single');
    f2 = (fread(fid, [row,col], 'int16'))./10000;
    fclose(fid);
end

% Bilinear interpolation
F = griddedInterpolant({lon,lat}, y_intercept', 'linear', 'linear'); a = F(x0, y0);
F = griddedInterpolant({lon,lat}, slope', 'linear', 'linear'); b = F(x0, y0);
F = griddedInterpolant({lon,lat}, a1', 'linear', 'linear'); amp1 = F(x0, y0);
F = griddedInterpolant({lon,lat}, f1', 'linear', 'linear'); pha1 = F(x0, y0);
F = griddedInterpolant({lon,lat}, a2', 'linear', 'linear'); amp2 = F(x0, y0);
F = griddedInterpolant({lon,lat}, f2', 'linear', 'linear'); pha2 = F(x0, y0);

% Surface pressure model
P = a + b*(mjd - 51178) + amp1*cos(2*pi*(mjd - 51178)/p1+pha1) + ...
    amp2*cos(2*pi*(mjd - 51178)/p2+pha2);

% Open and read the Tm coefficients and undulation file
[fid, errmsg] = fopen([coeffiles,'tm_grid.bin'], 'r');
if fid == -1
    error(errmsg)
else
    y_intercept = fread(fid, [row,col], 'single');
    slope = fread(fid, [row,col], 'single');
    orography = fread(fid, [row,col], 'single');
    undu = fread(fid, [row,col], 'single');
    fclose(fid);
end

% Bilinear interpolation
F = griddedInterpolant({lon,lat}, y_intercept', 'linear', 'linear'); a = F(x0, y0);
F = griddedInterpolant({lon,lat}, slope', 'linear', 'linear'); b = F(x0, y0);
F = griddedInterpolant({lon,lat}, orography', 'linear', 'linear'); geo_height = F(x0, y0);
F = griddedInterpolant({lon,lat}, undu', 'linear', 'linear'); N = F(x0, y0);

% Zenith hydrostatic delay (ZHD), Saastamoinen model
if strcmp(z0_type, 'orth')
    H_orth = z0;
elseif strcmp(z0_type, 'elli')
    H_orth = z0 - N;
else
    error('Use 1) <<orth>> for Orthometric height or 2) <<elli>> for Ellipsoidal height (in m).')
end

% Correction to P and T (see Guochanf Xu, GPS Theory, Algorithms and Applications, 2nd Edition, page 56)
% P = P * (1.0 - 0.000226*(H_orth - geo_height))^5.225;
% T = T - 0.0065*(H_orth - geo_height);

% Weight mean temperature, Tm
Tm = a + b*T;

% ZHD using the Saastamoinen Model (see Saastamoinen, 1973)
ZHD = (0.0022768 * P)/(1 - 0.0026*cos(2*deg2rad(y0))-0.00000028*H_orth);

return
