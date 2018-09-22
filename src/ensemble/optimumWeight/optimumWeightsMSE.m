function [MSE_optimumWeights weights LL_optimumWeights Obj_optimumWeights] = ...
                               optimumWeightsMSE(predictions, setsPositions, ...
                                     samplePositions, testingPositions, solveFor, data, ...
                                     weights_initial = ones(length(setsPositions),1)/length(setsPositions), ...
                                     lambda = 1, tol = sqrt(eps), maxiter = length(setsPositions)*4)
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
% File:   optimumWeightsMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 20, 2018.
%
%
% Function description:
% Combine the given sets by selecting the optimum weights, the ones that minimize
% the MSE_training
%
%
%
% Usage:
% [MSE_optimumWeights weights] = optimumWeightsMSE(predictions, setsPositions, ...
%                                      samplePositions, testingPositions, solveFor, ...
%                                      data, weights_initial, lambda, tol, maxiter)
%
% Input:
% predictions: a matrix with three-dimensions.The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".
%              The MSE is calculated for Th and Ts.
% setsPositions: the positions of the sets in the predictions variables.
% samplePositions: positions of the samples to use to add new best set
% testingPositions: positions to calculate the MSE_testing
% solveFor: Specify if should solve for all models or for a particular model.
%           Valid inputs are:
%            (1) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (2) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
% data: data that contains the ground truth.
% weights_initial: initial guess for weights
% lambda: balance between minimizing MSE and maximazing log likelihood
% tol: minimum tolerance value to assume convergence of the optimization functions.
% maxiter: number of maximum iterations to run the solver
%
% Output:
% MSE_optimumWeights: Matrix where the lines represents the quantity of sets in the greedy
%                     search. First column is MSE calculated using samplePositions and
%                     second column is MSE calculated using testingPositions
% weights: Vector that contains the weights for the setsPositions.
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

MSE_optimumWeights = ones(1,2)*Inf;

measuredValues = reshape(data(samplePositions, data_columns), [], 1);
setsPredictions = reshape( predictions(samplePositions, predictions_columns, ...
                                       setsPositions), [], length(setsPositions) );

if (lambda == 1)
  H = setsPredictions'*setsPredictions;
  q = -setsPredictions'*measuredValues;
  A = ones(1, length(weights_initial)) ;
  b = 1;
  options.MaxIter = maxiter;
  try
    [weights, Obj_optimumWeights, info, lambdaOptimization] = qp(weights_initial, H, q, A, b, ...
              zeros(length(weights_initial), 1), ones(length(weights_initial), 1), ...
              options );
          
    % Recasting the problem in case the convergence fails
    while (info.info == 3)
      warning(["lambda = 1: qp did not converge with " num2str(maxiter) " iterations. Recasting"])
      [weights, Obj_optimumWeights, info, lambdaOptimization] = qp(weights, H, q, A, b, ...
                zeros(length(weights_initial), 1), ones(length(weights_initial), 1), ...
                options );
    end
  catch
    warning("We couldn't find a solution for the optimization problem. Returning the initial weights");
    weights = weights_initial;
  end
else
  objFunction = @(weights)objectiveFunctionOptimumWeightMSE(weights, measuredValues, setsPredictions, lambda, length_MSE);
  gradObjFunction = @(weights)gradientObjectiveFunctionOptimumWeightMSE(weights, measuredValues, setsPredictions, lambda, length_MSE);
  HessObjFunction = @(weights)HessianObjectiveFunctionOptimumWeightMSE(weights, measuredValues, setsPredictions, lambda, length_MSE);
  
  eqConstFunction = @(weights)equalityConstraintOptimumWeightMSE(weights);
  gradEqConstFunction = @(weights)gradientEqualityConstraintOptimumWeightsMSE(weights);

  % since this function is very unstable, I'll try and if I fail I'll return the initial weights
  try
    [weights, Obj_optimumWeights, info, iter, nf, lambdaOptimization] = sqp(weights_initial, ...
      {objFunction, gradObjFunction, HessObjFunction}, {eqConstFunction, gradEqConstFunction}, [], 0, 1, maxiter, tol);%    {objFunction, gradObjFunction, HessObjFunction}, {eqConstFunction, gradEqConstFunction}, [], 0, 1, maxiter, tol);
  catch
      warning("We couldn't find a solution for the optimization problem. Returning the initial weights");
      weights = weights_initial;
   end
end

% Removing the weights that are below a tolerance
weights( find(weights < tol) ) = 0;
% renormalizing the weigths
weights = weights./sum(weights);

% getting the MSEs and log likelihood
mean_sets = setsPredictions*weights;
MSE_optimumWeights(1) = sum( ( measuredValues - ...
                                mean_sets ).^2 )/ ...
                             ( length_MSE );

variance_sets = (setsPredictions.^2)*weights - mean_sets.^2;
LL_optimumWeights = length_MSE/2*log(2*pi) - 1/2*sum( log(variance_sets) + ...
      1./variance_sets.*(measuredValues - mean_sets).^2 );
      

Obj_optimumWeights = lambda*MSE_optimumWeights(1) + (lambda - 1)*LL_optimumWeights;

setsPredictionsTesting = reshape( predictions(testingPositions, predictions_columns, ...
                                       setsPositions), [], length(setsPositions) );

measuredValuesTesting = reshape(data(testingPositions, data_columns), [], 1);

MSE_optimumWeights(2) = sum( ( measuredValuesTesting - ...
                                setsPredictionsTesting*weights ).^2)/ ...
                             ( length(measuredValuesTesting) );
