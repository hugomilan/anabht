function obj = objectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
                                                setsPredictions, lambda, length_MSE)
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
% File:   objectiveFunctionOptimumWeightMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Calculates the error in using a given weight
%
%
%
% Usage:
% obj = objectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
%                                         setsPredictions, lambda, length_MSE)
%
% Input:
% weights: Vector that contains the weights for the setsPredictions.
% measuredValues: a vector that contains the measured data. This is organized as [Th; Ts].
% setsPredictions: matrix that contains the predictions from the sets. The lines
%                  represent the predictions and the columns represent the 
%                  different sets.
% lambda: balance between minimizing MSE and maximazing log likelihood
% length_MSE: number of samples
%
% Output:
% obj: squared error.
%      

mean_sets = setsPredictions*weights;        

errorsMSE = measuredValues - mean_sets; % == measuredMinusMean
MSE = errorsMSE'*errorsMSE/length_MSE;

variance_sets = (setsPredictions.^2)*weights - mean_sets.^2;
LogLikelihood = length_MSE*0.5*log(2*pi) - 0.5*sum( log(variance_sets) + ...
                                                  errorsMSE.^2./variance_sets );
      
obj = lambda*MSE + (lambda - 1)*LogLikelihood;
