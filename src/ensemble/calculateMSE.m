function MSE = calculateMSE(sets, predictions, positions, solveFor, data)
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
% File:   calculateMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Calculates the MSE of the given positions of the specified predictions
%
%
%
% Usage:
% MSE = calculateMSE(sets, predictions, positions, solveFor data)
%
% Input:
% sets: a vector with the #sets.
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
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
%
% Output:
% MSE: matrix containing the MSE, where the first column is the MSE, the 
%      second column is the number of the set, and the third column is the position
%      of the set in the MSE.

n_sets = size(predictions, 3);
MSE = NA*zeros(n_sets, 2);
MSE(:, 3) = 1:n_sets;

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

for ii = 1:n_sets

  MSE(ii, 1) = sum( sum ( ( data(positions, data_columns) - predictions(positions, predictions_columns, ii) ).^2 ) )/ ...
                          ( length_MSE );
  MSE(ii, 2) = sets(ii);
end

end
