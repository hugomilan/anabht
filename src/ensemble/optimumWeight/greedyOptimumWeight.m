function [MSE_optimumWeights setsOptimumWeights greedyWeights ...
          LL_optimumWeights Obj_optimumWeights] = ...
                        greedyOptimumWeight(predictions, ...
                                            samplePositions, testingPositions, ...
                                            sets, solveFor, data, ...
                                            msg = "", lambda = 1, numberOfInitialSets = 20)
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
% File:   greedyOptimumWeight.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Starts combining the best (in minimum training MSE sense) numberOfInitialSets set 
% and checks if the next set improves upon the objective function. If
% it does, it is added to the ensemble. Keep examining next sets, adding the sets
% that improve MSE and removing the sets that end up with weight 0.
%
%
%
% Usage:
% [MSE_optimumWeights setsOptimumWeights greedyWeights LL_optimumWeights Obj_optimumWeights] = ...
%                         greedyOptimumWeight(predictions, ...
%                                             samplePositions, testingPositions, ...
%                                             sets, solveFor, data, msg, lambda, numberOfInitialSets)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% samplePositions: positions of the samples to use to add new best set
% testingPositions: positions to calculate the MSE_testing
% sets: sets to do the greedy search
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
% numberOfInitialSets: number of initial sets to do the optimization. This is 
%                      necessary to get started in the calculations because optimizing
%                      with log-likelihood and with a low number of samples leads to
%                      unstability.
%
% Output:
% MSE_optimumWeights: Matrix where the lines represents the quantity of sets in the greedy
%                     search. First column is MSE calculated using samplePositions and
%                     second column is MSE calculated using testingPositions
% setsOptimumWeights: Cell that contains the position of the sets selected during the
%                     greedy search.
% greedyWeights: cell that contains the weights for each greedy search.
% LL_optimumWeights: Vector that contains the log likelihood.
% Obj_optimumWeights: Vector that contains the value of the objective function.
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
         
totalNumberOfSets = length(sets);
MSE_optimumWeights = ones(totalNumberOfSets,2)*Inf;
LL_optimumWeights = ones(totalNumberOfSets,1)*Inf;
Obj_optimumWeights = ones(totalNumberOfSets,1)*Inf;
setsOptimumWeights = cell(1,totalNumberOfSets);
greedyWeights = cell(1,totalNumberOfSets);

% initial set. Will start with the first numberOfInitialSets
ii = numberOfInitialSets;
setToTest = sets(1:numberOfInitialSets);

[MSEToTest weightsToTest LLToTest ObjToTest] = optimumWeightsMSE(predictions, ...
              setToTest, samplePositions, testingPositions, ...
              solveFor, data, ones(numberOfInitialSets,1)/numberOfInitialSets, lambda);
              
MSE_optimumWeights(ii, :) = MSEToTest;
LL_optimumWeights(ii, :)  = LLToTest;
Obj_optimumWeights(ii, :) = ObjToTest;
% only adds the sets that had weight greater than zero
setsOptimumWeights{ii} = setToTest( find (weightsToTest > 0) );
greedyWeights{ii} = weightsToTest(  find (weightsToTest > 0) );

disp([msg "Greedy search with optimum weight (lambda = " num2str( lambda, "%2.4f" ) ") " num2str(ii) ...
                  ": MSE training: " num2str( MSE_optimumWeights(ii, 1), "%2.4f" ) ...
                  "; MSE testing: "  num2str( MSE_optimumWeights(ii, 2), "%2.4f" ) ...
                  "; Log likelihood: "  num2str( LL_optimumWeights(ii), "%2.4f" ) ...
                  "; Objective function: "  num2str( Obj_optimumWeights(ii), "%2.4f" ) ...
                  "; non-zero weights: " num2str( length(greedyWeights{ii}) )]);
                     
for ii = (numberOfInitialSets + 1):totalNumberOfSets
    setToTest = [setsOptimumWeights{ii - 1}; sets(ii)];
    [MSEToTest weightsToTest LLToTest ObjToTest] = optimumWeightsMSE(predictions, ...
              setToTest, samplePositions, testingPositions, ...
              solveFor, data, [greedyWeights{ii - 1}; 0], lambda);
                  
    if ( ObjToTest(1) < Obj_optimumWeights(ii - 1, 1) )
        MSE_optimumWeights(ii, :) = MSEToTest;
        LL_optimumWeights(ii, :)  = LLToTest;
        Obj_optimumWeights(ii, :) = ObjToTest;
        % only adds the sets that had weight greater than zero
        setsOptimumWeights{ii} = setToTest( find (weightsToTest > 0) );
        greedyWeights{ii} = weightsToTest(  find (weightsToTest > 0) );
    else
        MSE_optimumWeights(ii, :) = MSE_optimumWeights(ii - 1, :);
        LL_optimumWeights(ii, :)  = LL_optimumWeights( ii - 1, :);
        Obj_optimumWeights(ii, :) = Obj_optimumWeights(ii - 1, :);
        setsOptimumWeights{ii} = setsOptimumWeights{ii - 1};
        greedyWeights{ii} = greedyWeights{ii - 1};
    end
    
   
    disp([msg "Greedy search with optimum weight (lambda = " num2str( lambda, "%2.4f" ) ") " num2str(ii) ...
                  ": MSE training: " num2str( MSE_optimumWeights(ii, 1), "%2.4f" ) ...
                  "; MSE testing: "  num2str( MSE_optimumWeights(ii, 2), "%2.4f" ) ...
                  "; Log likelihood: "  num2str( LL_optimumWeights(ii), "%2.4f" ) ...
                  "; Objective function: "  num2str( Obj_optimumWeights(ii), "%2.4f" ) ...
                  "; non-zero weights: " num2str( length(greedyWeights{ii}) )]);
end

end
