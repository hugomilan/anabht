function const = equalityConstraintOptimumWeightMSE(weights)
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
% File:   equalityConstraintOptimumWeightMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Errors in satisfying the constraints
%
%
% Usage:
% const = equalityConstraintOptimumWeightMSE(weights)
%
% Input:
% weights: Vector that contains the weights for the setsPredictions.
%
% Output:
% const: Constraint that the sum of weights should be 1
% 
  
% error in the constraint
const = sum(weights) - 1;
