function [A, b, matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
        getProblemIndependentMatricesSkin(k, L, qtotal, wb, rhob, cb, ...
        To, qo)
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
% File:   getProblemIndependentMatricesSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 29, 2018.
%
%
% Function description:
% Calculates the matrices componentes that are independent of the problem. That is,
% each problem will change Ta, Tg, etc., and the matrices that this function outputs
% do not depend on those.
%
%
% Usage:
% [A, b, matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
%         getProblemIndependentMatricesSkin(k, L, qtotal, wb, rhob, cb, ...
%         To, qo)
%
% Input:
% k: vector column with thermal conductivity values for layers (W/(m oC))
% L: vector column with legth of layers (m)
% qtotal: vector column with total heat source values for layers (W/m3)
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% rhob: vector column with blood density values for layers (kg/m2)
% cb: vector column with blood specific heat values for layers (J/(kg oC))
% To: vector of offsets temperatures (oC)
% qo: vector of offsets heat fluxes (W/m2)
%
%
% Output:
% A: matrix containing the relationship between betas
% b: vector containing offset for the relationship between betas
% matrix_temp: matrix containing the relationship between betas and the temperatures
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
% 

% Matrix and vector to solve for beta
[A b] = getBetaParametersProblemIndependentSkin(k, L, qtotal, wb, rhob, cb, To, qo);

% Matrix and vector to calculate temperature given betas
[matrix_temp temp_offset] = getTempParametersProblemIndependent(k, L, qtotal, ...
                                                                wb, rhob, cb, To);

% Matrix and vector to calculate heat flux given betas        
[matrix_heat heat_offset] = getHeatFluxParametersProblemIndependent(k, L, qtotal, ...
                                                                wb, rhob, cb, qo);
        
end
