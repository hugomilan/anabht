function grad = gradientObjectiveFunctionSkin(beta, matrix_temp, temp_offset, ...
      matrix_heat, heat_offset, Ta, TMR, ua, d, epsilon, qo_end, z, Lt, h_m)
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
% File:   gradientOjectiveFunctionSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 2, 2018.
%
%
% Function description:
% Calculates the gradient of the objective funciton.
%
%
% Usage:
% grad = gradientObjectiveFunctionSkin(beta, matrix_temp, temp_offset, ...
%        matrix_heat, heat_offset, Ta, TMR, ua, d, epsilon, qo_end, z, Lt, h_m)
%
%
% Input:
% beta: vector of beta used in the model
% matrix_temp: matrix containing the relationship between betas and the temperatures
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
% Ta: air temperature (oC)
% TMR: Mean radiant temperature (oC)
% ua: air velocity (m/s). Default input parameter = 0 m/s
% d: diameter of the animal (m). Default input parameter = 0.15 m
% epsilon: emissivity of the animal surface
% qo_end: offset of heat flux at the skin
% z: elevation of the location (m).
% Lt: latitute of the location (decimal degrees).
% h_m: multiplicative factor for h to consider for possible measurement errors
% 
% Output:
% grad: Vector of the gradients
% 
% 

% Calculate the gradient of the objective function
grad = zeros(length(beta), 1);

Ts = matrix_temp(end,:)*beta + temp_offset(end);
qs = matrix_heat(end,:)*beta + heat_offset(end);

h = getConvHeatTransCoeff(Ta, Ts, ua, d, z, Lt);
h *= h_m;
qconv = h*(Ts - Ta);

sigma = 5.670373e-8;
qrad = sigma*epsilon*( (Ts + 273.15)^4 - (TMR + 273.15)^4 );

diff2 = 2*(qs - qo_end - qconv - qrad);
const = diff2*(h + 4*sigma*epsilon*(Ts + 273.15)^3);
grad(end-1) = diff2*matrix_heat(end,end-1) - const*matrix_temp(end,end-1);
grad(end) = diff2*matrix_heat(end,end) - const*matrix_temp(end,end);

end
