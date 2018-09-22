function grad = gradientEqualityConstraintOptimumWeightsMSE(weights)
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
% File:   gradientEqualityConstraintOptimumWeightsMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Calculates the gradient for the constraints
%
%
% Usage:
% grad = gradientEqualityConstraintOptimumWeightsMSE(weights)
%
% Input:
% weights: Vector that contains the weights for the setsPredictions.
%
% Output:
% grad: Gradient for weight_j
% 
  
% error in the constraint
grad = ones(1, length(weights));
