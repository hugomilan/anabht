function [ka, rhoa, ca, Da, Pa, upsilona, g] = getAirParameters(Ta, z = 595, Lt = 21.2583)
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
% File:   getAirParameters.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on February 9, 2018.
%
%
% Function description:
% Calculates air parameters given air temperature (Ta, oC),
% elevation (z, m), and latitude (Lt, decimal degrees).

%
% Usage:
% 1) Minimum usage. Uses altitude and latitude from Jaboticabal, Brazil
% [ka, rhoa, ca, Da, Pa, upsilona, g] = getAirParameters(Ta)
%
% 2) Maximum usage. User-defined altitude and latitude.
% [ka, rhoa, ca, Da, Pa, upsilona, g] = getAirParameters(Ta, z, Lt)
%
% Input:
% Ta: air temperature (oC)
% z: elevation of the location (m). Default input parameter = 595 m (for Jaboticabal, Brazil)
% Lt: latitute of the location (decimal degrees). Default input parameter = 21.2583 o (for Jaboticabal, Brazil)
% 
% Output:
% ka: air thermal conductivity (W/(m oC))
% rhoa: air density (kg/(m3))
% ca: air specific heat (J/(kg oC))
% Da: air thermal diffusivity (m2/s)
% Pa: air pressure (kPa)
% upsilona: air kinematic viscosity (m2/s)
% g: gravity (m/s2)
% 

% Calculations:

% gravity
g = 9.78013 + 8.18e-5*Lt + 1.168e-5*Lt^2 - 3.1e-6*z;

% air pressure
Pa = 101.325*exp( -z*g/( 287.04*( Ta + 273.15 ) ) );

% air density
rhoa = 3481.965*Pa/( Ta + 273.15 )*1e-3;

% air specific heat
ca = (1.0052 + 4.577e-4*exp( Ta/32.07733 ))*1e3;

% air thermal diffusivity
Da = 1.888e-5 + 1.324e-7*Ta;

% air thermal conductivity
ka = rhoa*ca*Da;

% air kinematic viscosity
upsilona = 1.32743e-5 + 9.22286e-8*Ta;

end