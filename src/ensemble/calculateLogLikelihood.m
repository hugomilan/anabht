function l = calculateLogLikelihood(predictions, setsPositions, positions, solveFor, data, ...
                            weights = ones(length(setsPositions),1)/length(setsPositions) )
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
% File:   calculateLogLikelihood.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on May 10, 2018.
%
%
% Function description:
% Calculates the log likelihood given the sets and their weights
%
%
%
% Usage:
% MSE = calculateLogLikelihood(predictions, setsPositions, positions, solveFor, data, weights)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
% setsPositions: the positions of the sets in the predictions variables.
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% positions = positions of the samples to use to calculate the MSE
% solveFor: Specify if should solve for all models or for a particular model.
%           Valid inputs are:
%            (1) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (2) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
% data: data that contains the ground truth.
% weights: weights to use in the calculations
%
% Output:
% l: log likelihood

n_sets = size(predictions, 3);

if ( strcmp(solveFor, "Skin") )
  data_columns = 9;
  predictions_columns = 2;
  length_MSE = length(positions);

elseif ( strcmp(solveFor, "KM") )
  data_columns = [14 9];
  predictions_columns = [1 2];
  length_MSE = 2*length(positions);
else
  error(["I don't know what to do with model " solveFor])
end

measuredValues = reshape(data(positions, data_columns), [], 1);
setsPredictions = reshape( predictions(positions, predictions_columns, ...
                                       setsPositions), [], length(setsPositions) );
                                       
mean_sets = setsPredictions*weights;
variance_sets = (setsPredictions.^2)*weights - mean_sets.^2;

l = length_MSE/2*log(2*pi) - 1/2*sum( log(variance_sets) + ...
      1./variance_sets.*(measuredValues - mean_sets).^2 );

end
