function [b, A_FD, b_FD, A_q_skin, b_q_skin, temp_offset] = ...
  updateProblemDependentMatricesKM(b, A_FD, b_FD1, A_q_skin, ...
  b_q_skin, temp_offset, k, wb, Tb, Tr, TMR, NTissueLayers)
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
% File:   updateProblemDependentMatricesKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Updates matrices for the specific parameters of a problem.
%
%
% Usage:
% [b, A_FD, b_FD, A_q_skin, b_q_skin, temp_offset] = ...
%  updateProblemDependentMatricesKM(b, A_FD, b_FD1, A_q_skin, ...
%  b_q_skin, temp_offset, k, wb, Tb, Tr, TMR, NTissueLayers)
%
% Input:
% b: vector containing offset for the relationship between betas
% A_FD: matrix that multiplies T for the KM model
% b_FD1: matrix that contains two vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% k: vector column with thermal conductivity values for layers (W/(m oC))
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% Tb: vector column with blood temperature values for layers (oC)
% Tr: Rectal temperature (oC)
% TMR: Mean radiant temperature (oC)
% NTissueLayers: number of tissue layers
%
% Output:
% b: vector containing offset for the relationship between betas
% A_FD: matrix that multiplies T for the KM model
% b_FD: matrix that contains the vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% 

% Matrix and vector to solve for beta
b = updateBetaParametersProblemDependent(b, wb(1:NTissueLayers), ...
        Tb(1:NTissueLayers), Tr);

% Matrices and vectors to solve for T inside the hair coat with radiation (KM model)
[A_FD, b_FD, A_q_skin, b_q_skin] = updateProblemDependentMatricesRadiationKM(TMR, ...
        k((NTissueLayers + 1):end), A_FD, b_FD1, A_q_skin, b_q_skin);

% Matrix and vector to calculate temperature given betas
temp_offset = updateTempParametersProblemDependent(temp_offset, ...
        wb(1:NTissueLayers), Tb(1:NTissueLayers));
        
end
