function obj = objectiveFunctionKM(T, A, E, b)
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
% File:   objectiveFunctionKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 07, 2018.
%
%
% Function description:
% Calculates the error in the KM model
%
%
% Usage:
% obj = objectiveFunctionKM(T, A, E, b)
%
% Input:
% T: Temperatures (oC)
% A: matrix that multiplies T for the KM model
% E: matrix that multiplies T^4 for the KM model
% b_FD: matrix that contains the vectors that will add to the KM model equation
%
% Output:
% obj: squared error of the KM model.
% 
  
errors = A*T + E*T.^4 + b;
% squared error
obj = errors'*errors;
