function [MSE_randomOptimumWeights randomWeights ...
          LL_randomOptimumWeights Obj_randomOptimumWeights] = ...
                        randomOptimumWeight(predictions, ...
                                            samplePositions, testingPositions, ...
                                            setsCellPositions, sets, solveFor, data, ...
                                            msg = "", lambda = 1)
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
% File:   randomOptimumWeight.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on May 12, 2018.
%
%
% Function description:
% Given the sets in the setsCells, this code calculates the optimum weights
% for the given sets.
%
%
%
% Usage:
% [MSE_randomOptimumWeights randomWeights LL_randomOptimumWeights Obj_randomOptimumWeights] = ...
%                         randomOptimumWeight(predictions, ...
%                                             samplePositions, testingPositions, ...
%                                             setsCellPositions, sets, solveFor, data, msg, lambda, numberOfInitialSets)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% samplePositions: positions of the samples to use to add new best set
% testingPositions: positions to calculate the MSE_testing
% setsCellPositions: cell that contains the sets to do the random search
% sets: vector that has the best sets (in minimum MSE training sense) as first sets
% solveFor: Specify if should solve for all models or for a particular model.
%           Valid inputs are:
%            (1) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (2) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
% data: data that contains the ground truth.
% msg: mensage to show at the beginning of the text
% lambda: balance between minimizing MSE and maximazing log likelihood
%
% Output:
% MSE_randomOptimumWeights: Matrix where the lines are for the different random 
%                           combinations. First column is MSE calculated using samplePositions and
%                           second column is MSE calculated using testingPositions
% randomWeights: cell that contains the weights for each random search.
% LL_randomOptimumWeights: Vector that contains the log likelihood.
% Obj_randomOptimumWeights: Vector that contains the value of the objective function.
%     

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
         
totalNumberOfSets = size(setsCellPositions,1);
MSE_randomOptimumWeights = ones(totalNumberOfSets,2)*Inf;
LL_randomOptimumWeights = ones(totalNumberOfSets,1)*Inf;
Obj_randomOptimumWeights = ones(totalNumberOfSets,1)*Inf;
randomWeights = cell(1,totalNumberOfSets);
                     
for ii = 1:totalNumberOfSets
    numberOfSetsTesting = setsCellPositions{ii, 1};
    [MSE_randomOptimumWeights(ii, :) randomWeights{ii} LL_randomOptimumWeights(ii) ...
        Obj_randomOptimumWeights(ii)] = optimumWeightsMSE(predictions, ...
              sets(setsCellPositions{ii, 2}), samplePositions, testingPositions, ...
              solveFor, data, ones(numberOfSetsTesting, 1)/numberOfSetsTesting, lambda);
    
   
    disp([msg "Random search with optimum weight (lambda = " num2str( lambda, "%2.4f" ) ") " num2str(ii) ...
                  ": MSE training: " num2str( MSE_randomOptimumWeights(ii, 1), "%2.4f" ) ...
                  "; MSE testing: "  num2str( MSE_randomOptimumWeights(ii, 2), "%2.4f" ) ...
                  "; Log likelihood: "  num2str( LL_randomOptimumWeights(ii), "%2.4f" ) ...
                  "; Objective function: "  num2str( Obj_randomOptimumWeights(ii), "%2.4f" ) ...
                  "; non-zero weights: " num2str( sum( randomWeights{ii} > 0) ) "/" ...
                  num2str( length( randomWeights{ii} ) ) ]);
end

end
