function MSE_bagging = calculateBaggingMSE(predictions, setsPositions, samplePositions, solveFor, data)
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
% File:   calculateBaggingMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Calculates the MSE for the bagging composed of the input sets that are defined
% by the positions of the variable "setsPositions". Bagging predictions are the
% average prediction of the sets' predictions.
%
%
%
% Usage:
% MSE_bagging = calculateBaggingMSE(predictions, setsPositions, samplePositions, solveFor, data)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% setsPositions: the positions of the sets in the predictions variables.
% samplePositions: positions of the samples to use the MSE calculation
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
% MSE_selected: sets selected

if ( strcmp(solveFor, "Skin") )
  data_columns = 9;
  predictions_columns = 2;
  length_MSE = length(samplePositions);

elseif ( strcmp(solveFor, "KM") )
  data_columns = [14 9];
  predictions_columns = [1 2];
  length_MSE = 2*length(samplePositions);

else
  error(["I don't know what to do with model " solveFor])
end

predictions_bagging = mean(predictions(samplePositions, predictions_columns, setsPositions), 3);

MSE_bagging = sum( sum ( ( data(samplePositions, data_columns) - predictions_bagging ).^2 ) )/ ...
                         ( length_MSE );

end
