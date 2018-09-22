function TMR = getTMR(Tg, Ta, ua = 0, dg = 0.15, epsilong = 0.95, z = 595, Lt = 21.2583)
%
% ANABHT - ANAlytical solver for steady-state BioHeat Transfer problems in 1D
% 
% Copyright (C) 2018 by Cornell University. All Rights Reserved.
% 
% Written by Hugo Fernando Maia Milan.
% 
% Free for educational, research and non-profit purposes.=
% Refer to the license file for details.
%
% 
% File:   getTMR.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on February 9, 2018.
%
%
% Function description:
% Calculates the mean randiant temperature (TMR, oC) given the black globe
% temperature (Tg, oC), air temperature (Ta, oC), air velocity (ua, m/s),
% black globe diameter (dg, m), radiant thermal emissivity of the globe (eg),
% elevation (z, m), and latitude (Lt, decimal degrees).
%
%
% Usage:
% 1) Minimum usage. Uses default values
% TMR = getTMR(Tg, Ta)
%
% 2) Maximum usage. User-defined input values
% TMR = getTMR(Tg, Ta, ua, dg, epsilong, z, Lt);
%
% Input:
% Tg: black globe temperature (oC)
% Ta: air temperature (oC)
% ua: air velocity (m/s). Default input parameter = 0 m/s
% dg: black globe diameter (m). Default input parameter = 0.15 m
% epsilong: radiant thermal emissivity of the globe. Default input parameter = 0.95
% z: elevation of the location (m). Default input parameter = 595 m (for Jaboticabal, Brazil)
% Lt: latitute of the location (decimal degrees). Default input parameter = 21.2583 o (for Jaboticabal, Brazil)
% 
% Output:
% TMR: Mean radiant temperature (oC)
% 

% Calculations:

% obtaining air parameters
[ka, rhoa, ca, Da, Pa, upsilona, g] = getAirParameters(Ta, z, Lt);

% Prandtl number
Pr = rhoa*ca*upsilona/ka;

% Reynolds number
Re = ua*dg/upsilona;

% Grashof number
Gr = g*dg^3*abs( Tg - Ta )/( upsilona^2*( Ta + 273.15 ) );

% relation between Grashof number and Reynolds number
if Re == 0 % avoiding division by zero
  xi = inf;
else
  xi = Gr/Re^2;
end



if xi <= 0.08
  % forced convection
  
  % air dynamic viscosity at air temperature
  mua = upsilona*rhoa;
  
  % air dynamic viscosity at black globe temperature
  [~, rhob, ~, ~, ~, upsilonb, ~] = getAirParameters(Tg, z, Lt);
  mub = upsilonb*rhob;
  
  % Nusselt number
  Nu = 2 + ( 0.4*Re^(1/2) + 0.06*Re^(2/3) )*Pr^0.4*( mua/mub )^(1/4);
elseif xi > 0.08 & xi < 3
  % mixed convection
  
  % Nusselt number for natural convection
  NuN = 2 + 0.589*( Gr*Pr )^(1/4)/( ( 1 + ( 0.469/Pr )^(9/16) )^(4/9) );
  
  % air dynamic viscosity at air temperature
  mua = upsilona*rhoa;
  
  % air dynamic viscosity at black globe temperature
  [~, rhob, ~, ~, ~, upsilonb, ~] = getAirParameters(Tg, z, Lt);
  mub = upsilonb*rhob;
  
  % Nusselt number for forced convection
  NuF = 2 + ( 0.4*Re^(1/2) + 0.06*Re^(2/3) )*Pr^0.4*( mua/mub )^(1/4);
  
  % Nusselt number
  Nu = ( NuN^4 + NuF^4 )^(1/4);
else
  % natural convection
  
  % Nusselt number
  Nu = 2 + 0.589*( Gr*Pr )^(1/4)/( ( 1 + ( 0.469/Pr )^(9/16) )^(4/9) );
  
end

% convection heat transfer coefficient for the black globe
hb = ka*Nu/dg;

% Stefan-Boltzmann constant
sigma = 5.670373e-8;

TMR = ( hb/( epsilong*sigma )*( Tg - Ta ) + ( Tg + 273.15 )^4 )^(1/4) - 273.15;

end