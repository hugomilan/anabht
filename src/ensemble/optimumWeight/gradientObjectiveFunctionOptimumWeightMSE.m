function grad = gradientObjectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
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
% File:   gradientObjectiveFunctionOptimumWeightMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Calculates the gradient for the model
%
%
% Usage:
% grad = gradientObjectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
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
% grad: Gradient for weight_j T_i
% 

mean_sets = setsPredictions*weights;

errorsMSE = measuredValues - mean_sets; % == measuredMinusMean
gradMSE = -2/length_MSE*setsPredictions'*errorsMSE;

%% These computations are correct
setsPredictionsSquared = setsPredictions.^2;
oneOverVariance = 1./( setsPredictionsSquared*weights - mean_sets.^2 );
firstDerivativeOfVariance = setsPredictionsSquared - 2*mean_sets.*setsPredictions;
errorDividedByVariance = errorsMSE.*oneOverVariance;
errorSquaredDividedByVarianceSquared = errorDividedByVariance.^2;

gradLL = 0.5*( errorSquaredDividedByVarianceSquared -  oneOverVariance)'*firstDerivativeOfVariance + ...
         errorDividedByVariance'*setsPredictions;

grad = lambda*gradMSE + (lambda - 1)*gradLL';

%% Old way to compute gradient for log-likelihood
%variance_sets = (setsPredictions.^2)*weights - mean_sets.^2;
%setsMinusTwoMean = setsPredictions - 2*mean_sets;
%% measuredMinusMean = measuredValues - mean_sets;
%oneMinusMeasuredMinusMeanSquaredDividedByVariance = 1 - errorsMSE.^2./variance_sets;
% 
%gradLL = zeros(length(weights), 1);
%for kk = 1:length(weights)
%  gradLL(kk) = -0.5*sum( setsPredictions(:, kk)./variance_sets.*( ...
%    setsMinusTwoMean(:, kk).*oneMinusMeasuredMinusMeanSquaredDividedByVariance ...
%    -2*errorsMSE ) );
%end
%
%grad = lambda*gradMSE + (lambda - 1)*gradLL;

end
