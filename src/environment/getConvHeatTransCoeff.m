function [h] = getConvHeatTransCoeff(Ta, Ts, ua = 0, d = 0.15, z = 595, Lt = 21.2583)
%
% ANABHT - ANAlytical solver for steady-state BioHeat Transfer problems in 1D
% 
% Copyright (C) 2018 by Cornell University. All Rights Reserved.
% % Written by Hugo Fernando Maia Milan.
% 
% Free for educational, research and non-profit purposes.
% Refer to the license file for details.
%
% 
% File:   getConvHeatTransCoeff.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on February 9, 2018.
%
%
% Function description:
% Calculates the convection heat transfer coefficient given animal's
% and environmental parameters.
%
% Usage:
% 1) Minimum usage. Use default input values
% h = getConvHeatTransCoeff(Ta, Ts)
%
% 1) Maximum usage. Specify all input parameters
% h = getConvHeatTransCoeff(Ta, Ts, ua, d, z, Lt)
%
% Input:
% Ta: air temperature (oC)
% Ts: surface temperature (oC)
% ua: air velocity (m/s). Default input parameter = 0 m/s
% d: diameter of the animal (m). Default input parameter = 0.15 m
% z: elevation of the location (m). Default input parameter = 595 m (for Jaboticabal, Brazil)
% Lt: latitute of the location (decimal degrees). Default input parameter = 21.2583 o (for Jaboticabal, Brazil)
% 
% Output:
% h: convection heat transfer coefficient (W/m2)
% 
% Assumptions:
% Animal is a horizontal cylinder and airflow is perpendicular to the cylinder.
%

% Calculations:

% obtaining air parameters
[ka, rhoa, ca, Da, Pa, upsilona, g] = getAirParameters(Ta, z, Lt);

% Prandtl number
Pr = rhoa*ca*upsilona/ka;

% Reynolds number
Re = ua*d/upsilona;

% Grashof number
Gr = g*d^3*abs( Ts - Ta )/( upsilona^2*( Ta + 273.15 ) );

% relation between Grashof number and Reynolds number
if Re == 0 % avoiding division by zero
  xi = inf;
else
  xi = Gr/Re^2;
end


if xi <= 0.08
  % forced convection
  
  % Nusselt number
  Nu = 0.3 + 0.62*Re^(1/2)*Pr^(1/3)/( ( 1 + ( 0.4/Pr )^(2/3) )^(1/4) )*( 1 + ( Re/282000 )^(5/8) )^(4/5);
elseif xi > 0.08 & xi < 3
  % mixed convection
  
  % Nusselt number for natural convection
  NuN = ( 0.6 + 0.387*( Gr*Pr )^(1/6)/( ( 1 + ( 0.559/Pr )^(9/16) )^(8/27) ) )^2;
  
  % Nusselt number for forced convection
  NuF =  0.3 + 0.62*Re^(1/2)*Pr^(1/3)/( ( 1 + ( 0.4/Pr )^(2/3) )^(1/4) )*( 1 + ( Re/282000 )^(5/8) )^(4/5);
  
  % Nusselt number
  Nu = ( NuN^3.5 + NuF^3.5 )^(1/3.5);
else
  % natural convection
  
  % Nusselt number
  Nu = ( 0.6 + 0.387*( Gr*Pr )^(1/6)/( ( 1 + ( 0.559/Pr )^(9/16) )^(8/27) ) )^2;
end

% coefficient of convective heat transfer
h = ka*Nu/d;

end