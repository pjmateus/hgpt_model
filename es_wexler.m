function es = es_wexler(T, P)
% Saturation Water Vapor Pressure (es) using Wexler formulation with new coefficients (adjusted for ITS-90)
% INPUT : T = Temperature (in Kelvins)
%         P = atmospheric Pressure (in hPa)
% OUTPUT: es= saturation water vapor pressure (in hPa)
%
% References:
% 1)
% Wexler, A. Vapor Pressure Formulation for Water in Range 0 to 100 Degrees C. A Revision. 
%   J. Res. Natl. Bur. Stand. 1976, 80A, 775–785.
% 2)
% Wexler, A. Vapor Pressure Formulation for Ice. 
%   J. Res. Natl. Bur. Stand. 1977, 81A, 5–20.
% 3)
% Hardy, B. ITS-90 Formulations for Water Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature, and Enhancement Factors in Range -100 to +100 C. 
%   In Proceedings of the Third International Symposium on Humidity and Moisture; 
%   UK National Physical Laboratory (NPL): Teddington, UK, April 6 1998; pp. 1–8.
%
if T >= 273.15    
    % Saturation Vapor Pressure over Water
    g0 =-2.8365744*10^3;
    g1 =-6.028076559*10^3;
    g2 = 1.954263612*10^1;
    g3 =-2.737830188*10^-2;
    g4 = 1.6261698*10^-5;
    g5 = 7.0229056*10^-10;
    g6 =-1.8680009*10^-13;
    g7 = 2.7150305;
    
    es = 0.01 * exp(g0*T^-2 + g1*T^-1 + g2 + g3*T + g4*T^2 + g5*T^3 + g6*T^4 + g7*log(T));
    
    % Enhancement Factor coefficients for Water 0 to 100°C
    A0 =-1.6302041*10^-1;
    A1 = 1.8071570*10^-3;
    A2 =-6.7703064*10^-6;
    A3 = 8.5813609*10^-9;
    B0 =-5.9890467*10^1;
    B1 = 3.4378043*10^-1;
    B2 =-7.7326396*10^-4;
    B3 = 6.3405286*10^-7;    
else    
    % Saturation Vapor Pressure over Ice
    k0 =-5.8666426*10^3;
    k1 = 2.232870244*10^1;
    k2 = 1.39387003*10^-2;
    k3 =-3.4262402*10^-5;
    k4 = 2.7040955*10^-8;
    k5 = 6.7063522*10^-1;
    
    es = 0.01 * exp(k0*T^-1 + k1 + k2*T + k3*T^2 + k4*T^3 + k5*log(T));
    
    % Enhancement Factor coefficients for Ice –100 to 0°C
    A0 =-6.0190570*10^-2;
    A1 = 7.3984060*10^-4;
    A2 =-3.0897838*10^-6;
    A3 = 4.3669918*10^-9;
    B0 =-9.4868712*10^1;
    B1 = 7.2392075*10^-1;
    B2 =-2.1963437*10^-3;
    B3 = 2.4668279*10^-6;
end

% Enhancement Factor
alpha = A0 + A1*T + A2*T^2 + A3*T^3;
beta = exp(B0 + B1*T + B2*T^2 + B3*T^3);
f = exp( alpha*(1-es/P) + beta*(P/es-1) );

es = es * f;
return
