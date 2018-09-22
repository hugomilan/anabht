function [A, b, A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin, ...
  matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
  getProblemIndependentMatricesKM(k, L, qtotal, wb, rhob, cb, ...
  To, qo, NTissueLayers, alpha_eff_KM)
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
% File:   getProblemIndependentMatricesKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrices componentes that are independent of the problem. That is,
% each problem will change Ta, Tg, etc., and the matrices that this function outputs
% do not depend on those.
%
%
% Usage:
% [A, b, A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin, ...
%  matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
%  getProblemIndependentMatricesKM(k, L, qtotal, wb, rhob, cb, ...
%  To, qo, NTissueLayers, alpha_eff_KM)
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
% NTissueLayers: number of tissue layers
% alpha_eff_KM: effective absorption coefficient of the hair coat for the KM model
%
%
% Output:
% A: matrix containing the relationship between betas
% b: vector containing offset for the relationship between betas
% A_FD: matrix that multiplies T for the KM model
% E: matrix that multiplies T^4 for the KM model
% b_FD1: matrix that contains two vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% E_q_skin: vector that multiplies T^4 to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% matrix_temp: matrix containing the relationship between betas and the temperatures
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
% 

% Matrix and vector to solve for beta
[A b] = getBetaParametersProblemIndependentKM(k(1:NTissueLayers), L(1:NTissueLayers), ...
        qtotal(1:NTissueLayers), wb(1:NTissueLayers), rhob(1:NTissueLayers), ...
        cb(1:NTissueLayers), To(1:(NTissueLayers + 1)), ...
        qo(1:(NTissueLayers + 1)) );

% Matrices and vectors to solve for T inside the hair coat with radiation (KM model)
deltax = L(end);
[A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin] = ...
       getProblemIndependentMatricesRadiationKM(alpha_eff_KM, deltax, qtotal((NTissueLayers + 1):end));

% Matrix and vector to calculate temperature given betas
[matrix_temp temp_offset] = getTempParametersProblemIndependent(k(1:NTissueLayers), L(1:NTissueLayers), ...
        qtotal(1:NTissueLayers), wb(1:NTissueLayers), rhob(1:NTissueLayers), ...
        cb(1:NTissueLayers), To(1:(NTissueLayers + 1)));

% Matrix and vector to calculate heat flux given betas        
[matrix_heat heat_offset] = getHeatFluxParametersProblemIndependent(k(1:NTissueLayers), L(1:NTissueLayers), ...
        qtotal(1:NTissueLayers), wb(1:NTissueLayers), rhob(1:NTissueLayers), ...
        cb(1:NTissueLayers), qo(1:(NTissueLayers + 1)) );
        
end
