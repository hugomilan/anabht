function [MSE_bagging setsBagging LL_bagging obj_bagging] = greedyBagging(predictions, totalNumberOfSets, ...
                                     samplePositions, testingPositions, sets, solveFor, data, msg = "", lambda = 1)
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
% File:   greedyBagging.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Adds one set at a time. Search for all sets in the given list of sets, which is
% in the variable predictions, and adds one set that minimize the MSE using samplePositions.
% Returns the MSE using samplePosition and the MSE for the testing positions.
%
%
%
% Usage:
% [MSE_bagging setsBagging LL_bagging obj_bagging] = greedyBagging(predictions, totalNumberOfSets, ...
%                                     samplePositions, testingPositions, sets, solveFor, data, lambda)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% totalNumberOfSets: the total number of sets to add during the greedy search.
% samplePositions: positions of the samples to use to add new best set
% testingPositions: positions to calculate the MSE_testing
% sets: set IDs to do the greedy search.
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
% MSE_bagging: Matrix where the lines represents the quantity of sets in the greedy
%              search. First column is MSE calculated using samplePositions and
%              second column is MSE calculated using testingPositions
% setsBagging: Vector that contains the position of the sets selected during the
%              greedy search.
% LL_bagging: vector that contains the log likelihood.
% obj_bagging: vector that contains the objective function values.
%

if ( strcmp(solveFor, "Skin") )
  data_columns = 9;
  predictions_columns = 2;
  length_MSE = length(samplePositions);
  length_MSE_testing = length(testingPositions);

elseif ( strcmp(solveFor, "KM") )
  data_columns = [14 9];
  predictions_columns = [1 2];
  length_MSE = 2*length(samplePositions);
  length_MSE_testing = 2*length(testingPositions);

else
  error(["I don't know what to do with model " solveFor])
end

setsBagging = [];
MSE_bagging = ones(totalNumberOfSets,2)*Inf;
LL_bagging = ones(totalNumberOfSets,1)*-Inf;
obj_bagging = ones(totalNumberOfSets,1)*Inf;
predictions_greedy = 0;
predictions_greedy_squared = 0;
measuredValues = data(samplePositions, data_columns);

% initial set
ii = 1;
predictions_greedy = predictions(samplePositions, predictions_columns, sets(ii));
predictions_greedy_squared = predictions_greedy.^2;
setsBagging(end + 1,1) = sets(ii);
MSE_bagging(ii, 1) = sum( sum ( ( measuredValues - predictions_greedy ).^2 ) )/ ...
                                ( length_MSE );
LL_bagging(ii) =  NA;
obj_bagging(ii) = NA;
% calculating MSE_bagging on the testingPositions
predictions_greedy_testing = mean(predictions(testingPositions, predictions_columns, sets(ii)), 3);
    
MSE_bagging(ii, 2) = sum( sum ( ( data(testingPositions, data_columns) - predictions_greedy_testing ).^2 ) )/ ...
                                ( length_MSE_testing );
                            
disp([msg "greedy search with bagging (lambda = " num2str( lambda, "%2.4f" ) ") " num2str(ii) ...
                 ": Objective: " num2str( obj_bagging(ii), "%2.4f" ) ...
                 "; MSE training: " num2str( MSE_bagging(ii, 1), "%2.4f" ) ...
                 "; MSE testing: "  num2str( MSE_bagging(ii, 2), "%2.4f" ) ...
                 "; Log likelihood: "  num2str( LL_bagging(ii), "%2.4f" )]);

for ii = 2:totalNumberOfSets
  MSEToAdd = Inf;
  LLToAdd = Inf;
  objToAdd = Inf;
  setToAdd = NA;
  predictions_greedyToAdd = 0;
  predictions_greedy_squaredToAdd = 0;
  length_setsPositions = length(setsBagging) + 1;
  
  % I'll define two ways to calculate the objective function. I need to define them
  % because log likelihood is -Inf for ii == 1
  if (ii == 1)
    objFunction = @(MSE, LL) MSE;
  else
    objFunction = @(MSE, LL) lambda*MSE + (lambda - 1)*LL;
  end
  
  for kk = 1:length(sets)
    set_j = sets(kk);
    % if set_j is already in the selected set, we skip it
    if ( sum(set_j == setsBagging) )
        continue
    end
    setsPositions = [setsBagging; set_j];
        
    predictions_greedy_test = ( predictions_greedy*( length_setsPositions - 1) + ...
                                predictions(samplePositions, predictions_columns, set_j) ) ...
                                /( length_setsPositions );
                                
    predictions_greedy_squared_test = ( predictions_greedy_squared*( length_setsPositions - 1) + ...
                                        predictions(samplePositions, predictions_columns, set_j).^2 ) ...
                                       /( length_setsPositions );
        
    MSE_newSet = sum( sum ( ( measuredValues - predictions_greedy_test ).^2 ) )/ ...
                            ( length_MSE );
                            
    % calculating log likelihood
    variance_sets = predictions_greedy_squared_test - predictions_greedy_test.^2;
    LL_newSet = length_MSE/2*log(2*pi) - 1/2*sum( sum( log(variance_sets) +  ...
                   1./variance_sets.*(measuredValues - predictions_greedy_test).^2 ) );
                            
    obj_newSet = objFunction(MSE_newSet, LL_newSet);
                             
    if (obj_newSet < objToAdd)
      objToAdd = obj_newSet;
      MSEToAdd = MSE_newSet;
      LLToAdd = LL_newSet;
      setToAdd = set_j;
      predictions_greedyToAdd = predictions_greedy_test;
      predictions_greedy_squaredToAdd = predictions_greedy_squared_test;
    end
  end
    
  % now we have a new position to add
  setsBagging(end + 1,1) = setToAdd;
  MSE_bagging(ii, 1) = MSEToAdd;
  LL_bagging(ii) = LLToAdd;
  obj_bagging(ii) = objToAdd;
  predictions_greedy = predictions_greedyToAdd;
  predictions_greedy_squared = predictions_greedy_squaredToAdd;
  % calculating MSE_bagging on the testingPositions
  predictions_greedy_testing = mean(predictions(testingPositions, predictions_columns, setsBagging), 3);
    
  MSE_bagging(ii, 2) = sum( sum ( ( data(testingPositions, data_columns) - predictions_greedy_testing ).^2 ) )/ ...
                                  ( length_MSE_testing );
                            
 disp([msg "greedy search with bagging (lambda = " num2str( lambda, "%2.4f" ) ") " num2str(ii) ...
                   ": Objective: " num2str( obj_bagging(ii), "%2.4f" ) ...
                   "; MSE training: " num2str( MSE_bagging(ii, 1), "%2.4f" ) ...
                   "; MSE testing: "  num2str( MSE_bagging(ii, 2), "%2.4f" ) ...
                   "; Log likelihood: "  num2str( LL_bagging(ii), "%2.4f" )]);
end


end
