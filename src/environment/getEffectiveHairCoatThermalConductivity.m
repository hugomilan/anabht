function keff = getEffectiveHairCoatThermalConductivity(Ta, N, D, HL, Lh, kh, z = 595, Lt = 21.2583)
%
% ANABHT - ANAlytical solver for steady-state BioHeat Transfer problems in 1D
% 
% Copyright (C) 2018 by Cornell University. All Rights Reserved.
% 
% Written by Hugo Fernando Maia Milan.
% 
% Free for educational, research and non-profit purposes.
% Refer to the license file for details.
%
% 
% File:   getEffectiveHairCoatThermalConductivity.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on February 9, 2018.
%
%
% Function description:
% Calculates the effective hair coat thermal conductivity given environmental
% parameters and hair coat characteristics.
%
%
% Usage:
% 1) Minimum usage. Uses altitude and latitude from Jaboticabal, Brazil
% keff = getHairThermalConductivity(Ta, N, D, HL, Lh, kh)
%
% 2) Maximum usage. User-defined altitude and latitude.
% keff = getHairThermalConductivity(Ta, N, D, HL, Lh, kh, z, Lt)
%
% Input:
% Ta: air temperature (oC)
% N: number of hair per m2
% D: hair diameter (m)
% HL: hair length (m)
% Lh: hair coat thickness (m)
% z: elevation of the location (m). Default input parameter = 595 m (for Jaboticabal, Brazil)
% Lt: latitute of the location (decimal degrees). Default input parameter = 21.2583 o (for Jaboticabal, Brazil)
% 
% Output:
% keff: effective hair coat thermal conductivity (W/(m oC))
% 

% Calculations:

% obtaining air parameters
[ka, ~, ~, ~, ~, ~, ~] = getAirParameters(Ta, z, Lt);

% AF/AT: Cross section area of hairs per unit area. Cannot be greater than 1
AF_AT = N*HL/Lh*pi*D^2/4;
if (AF_AT > 1)
  AF_AT = 1; % correcting for overflows
end

% lc: average fiber spacing. Cannot be lower than D
lc = 1/sqrt(N*HL/Lh);
if (lc < D)
  lc = D; % correcting for overflows
end

% keff
keff = 0.5*( AF_AT*(kh - ka) + ka ) + ... % thermal conductivity in the direction perpendicular to the skin
       0.5*( ka*(lc - D)/lc + D*ka*kh/( D*ka +(lc - D)*kh ) ); % thermal conductivity in the direction parallel to the skin

end
