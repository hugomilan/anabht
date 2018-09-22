function Hessian = HessianObjectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
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
% File:   HessianObjectiveFunctionOptimumWeightMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Calculates the Hessian for the model
%
%
% Usage:
% Hessian = HessianObjectiveFunctionOptimumWeightMSE(weights, measuredValues, ...
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
% Hessian: Hessian matrix
% 

HessianMSE = 2/length_MSE*setsPredictions'*setsPredictions;

%% These computations are correct
mean_sets = setsPredictions*weights;
errorsMSE = measuredValues - mean_sets;
setsPredictionsSquared = setsPredictions.^2;
oneOverVariance = 1./( setsPredictionsSquared*weights - mean_sets.^2 );
oneOverVarianceSquared = oneOverVariance.^2;
firstDerivativeOfVariance = setsPredictionsSquared - 2*mean_sets.*setsPredictions;
errorDividedByVarianceSquared = errorsMSE.*oneOverVarianceSquared;
errorSquaredDividedByVarianceSquared = errorsMSE.*errorDividedByVarianceSquared;
errorSquaredDividedByVarianceCubic = errorSquaredDividedByVarianceSquared.*oneOverVariance;

HessianLL = ( -2*errorDividedByVarianceSquared.*setsPredictions + ...
              (0.5*oneOverVarianceSquared - errorSquaredDividedByVarianceCubic).*firstDerivativeOfVariance )'...
              *firstDerivativeOfVariance ...
            - ( errorSquaredDividedByVarianceSquared.*setsPredictions )'*setsPredictions;

%%% Old way to compute Hessian for log-likelihood
%mean_sets = setsPredictions*weights;  
%variance_sets = (setsPredictions.^2)*weights - mean_sets.^2;
%setsMinusTwoMean = setsPredictions - 2*mean_sets;
%measuredMinusMean = measuredValues - mean_sets;
%twoMeasuredMinusMeanSquaredDividedByVarianceMinusOne = 2*measuredMinusMean.^2./variance_sets - 1;
%
%HessianLL = zeros(length(weights), length(weights));
%for ll = 1:length(weights)
%  for kk = ll:length(weights)
%    HessianLL(ll, kk) = -0.5*sum( setsPredictions(:, kk).*setsPredictions(:, ll)./variance_sets.^2.*( ...
%                             setsMinusTwoMean(:, ll).*setsMinusTwoMean(:, kk).*(...
%                             twoMeasuredMinusMeanSquaredDividedByVarianceMinusOne ) + ...
%                             2*measuredMinusMean.*( measuredMinusMean + ...
%                             setsMinusTwoMean(:, ll) + setsMinusTwoMean(:, kk) ) ) );
%    
%    
%    % taking advantage of the symmetry
%    HessianLL(kk, ll) = HessianLL(ll, kk);
%  end
%end

Hessian = lambda*HessianMSE + (lambda - 1)*HessianLL;

end
