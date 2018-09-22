function [errors] = equalityConstraintFunctionSkin(beta, A, b)
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
% File:   equalityConstraintFunctionSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 1, 2018.
%
%
% Function description:
% Applies equality constraints to be solved by the sqp algorithm.
%
% Usage:
% [errors] = equalityConstraintFunctionSkin(beta, A, b)
%
% Input:
% beta: vector of beta used in the model
% A: Matrix that multiplies beta
% b: Vector that should be equal to A*beta
% 
% Output:
% errors: Vector of the errors in the equality
% 

% calculate the errors in the constrains
errors = A*beta - b;

end
