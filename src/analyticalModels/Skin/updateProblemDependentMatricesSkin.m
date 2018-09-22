function [b, temp_offset] = updateProblemDependentMatricesSkin(b, temp_offset, wb, Tb, Tr)
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
% File:   updateProblemDependentMatricesSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 27, 2018.
%
%
% Function description:
% Updates matrices for the specific parameters of a problem.
%
%
% Usage:
% [b, temp_offset] = updateProblemDependentMatricesSkin(b, temp_offset, wb, Tb, Tr)
%
% Input:
% b: vector containing offset for the relationship between betas
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% Tb: vector column with blood temperature values for layers (oC)
% Tr: Rectal temperature (oC)
%
% Output:
% b: vector containing offset for the relationship between betas
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% 

% Matrix and vector to solve for beta
b = updateBetaParametersProblemDependent(b, wb, Tb, Tr);

% Matrix and vector to calculate temperature given betas
temp_offset = updateTempParametersProblemDependent(temp_offset, wb, Tb);
        
end
