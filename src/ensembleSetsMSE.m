function ensembleSetsMSE(sets, solveFor = "all", training = 1, maxSetsToCombine = 1000, ...
                         numberOfInitialSetsOptimumWeights = 20, ...
                         minSetsToRandomlyCombineOptimumWeight = 5, ...
                         maxSetsToRandomlyCombineOptimumWeight = 20, ...
                         lambdaVector = [0 10.^(-4:0.25:-1) 0.2:0.3:0.8 (1 - 10.^(-1:-0.25:-4)) 1])
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
% File:   ensembleSetsMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on May 5, 2018.
%
%
% Function description:
% Given a group of #set, calls the functions that creates ensamble out of them
%
%
% Usage:
% ensembleSetsMSE(sets, solveFor, training, maxSetsToCombine, numberOfInitialSetsOptimumWeights, ...
%                 minSetsToRandomlyCombineOptimumWeight, maxSetsToRandomlyCombineOptimumWeight, ...
%                 lambdaVector);
% 
%
% Input:
% sets: vector with #set to combine
% solveFor: Specify if should solve for all models or for a particular model.
%           Valid inputs are:
%            (1) "all": loads solutions from all models.
%            (2) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (3) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
% training: Should calculate MSE for training, testing, or both? If 0, calculates 
%           MSE for testing. If 1, calculates MSE for training. If 2, calculates MSE for
%           both, testing and training. WARNING: it assumes that the requested dataset
%           (e.g., training or testing) were already used to calculate the predictions.
% maxSetsToCombine: maximum number of sets to combine.
% numberOfInitialSetsOptimumWeights: Number of initial sets to combine when doing
%                                    optimum weights calculation for naive and
%                                    greedy searchers. This speeds up the calculations.
% minSetsToRandomlyCombineOptimumWeight: Minimum numbers of sets to randomly combine
%                                        when calculating optimum weights.
% maxSetsToRandomlyCombineOptimumWeight: Maximum numbers of sets to randomly combine
%                                        when calculating optimum weights.
% lambdaVector: vector with the lambda values to test. lambda is a balance 
%               between minimizing MSE and maximazing log likelihood.
% 
% Output:
% None
%


disp("Loading data")
timei1 = time;

if ( strcmp(solveFor, "all") )
  ensembleSetsMSE(sets, maxSetsToCombine, "KM", training);
  ensembleSetsMSE(sets, maxSetsToCombine, "Skin", training);
  return
end

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  pkg load statistics
  pkg load parallel
  more off
else
  % we are executing in MATLAB
end

predictions = loadPredictions(sets, solveFor);

% loading positions
folder_datasets = locateFolderPath("datasets");
load([folder_datasets "/crossValidationPosition.mat"]) % cross-validation 
load([folder_datasets "/trainingPosition.mat"])

addpath(folder_datasets);
data_piglets %loads ground truth data

%% Adding the necessary paths
folder_ensemble = locateFolderPath("ensemble");
folder_ensembleData = locateFolderPath("ensembleData");
addpath(folder_ensemble);
folder_ensembleBagging = locateFolderPath("ensemble/bagging");
addpath(folder_ensembleBagging);
folder_ensembleOptimumWeight = locateFolderPath("ensemble/optimumWeight");
addpath(folder_ensembleOptimumWeight);
% get a group of sets
MSE_training = calculateMSE(sets, predictions, trainingPosition, solveFor, data);
MSE_cv1 = calculateMSE(sets, predictions, cross_validation_1, solveFor, data);
MSE_cv2 = calculateMSE(sets, predictions, cross_validation_2, solveFor, data);
MSE_cv3 = calculateMSE(sets, predictions, cross_validation_3, solveFor, data);
MSE_cv4 = calculateMSE(sets, predictions, cross_validation_4, solveFor, data);
MSE_cv5 = calculateMSE(sets, predictions, cross_validation_5, solveFor, data);
% removing NAs
NAsTraining = find( isnan(MSE_training(:, 1) ) );
numberOfNonNanSets = size(MSE_training, 1) - length(NAsTraining);
MSE_cv1(NAsTraining, 1) = nan;
MSE_cv2(NAsTraining, 1) = nan;
MSE_cv3(NAsTraining, 1) = nan;
MSE_cv4(NAsTraining, 1) = nan;
MSE_cv5(NAsTraining, 1) = nan;

% sorting the MSEs
MSE_cv1s = sortrows(MSE_cv1, 1);
MSE_cv2s = sortrows(MSE_cv2, 1);
MSE_cv3s = sortrows(MSE_cv3, 1);
MSE_cv4s = sortrows(MSE_cv4, 1);
MSE_cv5s = sortrows(MSE_cv5, 1);
MSE_trainings = sortrows(MSE_training, 1);


timeToLoad = time - timei1;
save([folder_ensembleData "/MSEs" solveFor ".mat"], "-V7", "MSE_training", "MSE_trainings", ...
     "MSE_cv1", "MSE_cv1s", "MSE_cv2", "MSE_cv2s", "MSE_cv3", "MSE_cv3s", ...
     "MSE_cv4", "MSE_cv4s", "MSE_cv5", "MSE_cv5s", "predictions", "NAsTraining", "numberOfNonNanSets", ...
     "maxSetsToCombine", "lambdaVector", "timeToLoad");
% load([folder_ensembleData "/MSEs" solveFor ".mat"])


% Naive search and bagging
timei2 = time;
MSEBagging = NA*zeros(maxSetsToCombine, 7);
LLBagging = NA*zeros(maxSetsToCombine, 7); % log likelihood

for ii = 1:maxSetsToCombine
  MSEBagging(ii, 1) = calculateBaggingMSE(predictions, MSE_cv1s(1:ii, 3), cross_validation_1c, solveFor, data);
  LLBagging(ii, 1)  = calculateLogLikelihood(predictions, MSE_cv1s(1:ii, 3), cross_validation_1, solveFor, data);
  
  MSEBagging(ii, 2) = calculateBaggingMSE(predictions, MSE_cv2s(1:ii, 3), cross_validation_2c, solveFor, data);
  LLBagging(ii, 2)  = calculateLogLikelihood(predictions, MSE_cv2s(1:ii, 3), cross_validation_2, solveFor, data);
  
  MSEBagging(ii, 3) = calculateBaggingMSE(predictions, MSE_cv3s(1:ii, 3), cross_validation_3c, solveFor, data);
  LLBagging(ii, 3)  = calculateLogLikelihood(predictions, MSE_cv3s(1:ii, 3), cross_validation_3, solveFor, data);
  
  MSEBagging(ii, 4) = calculateBaggingMSE(predictions, MSE_cv4s(1:ii, 3), cross_validation_4c, solveFor, data);
  LLBagging(ii, 4)  = calculateLogLikelihood(predictions, MSE_cv4s(1:ii, 3), cross_validation_4, solveFor, data);
  
  MSEBagging(ii, 5) = calculateBaggingMSE(predictions, MSE_cv5s(1:ii, 3), cross_validation_5c, solveFor, data);
  LLBagging(ii, 5)  = calculateLogLikelihood(predictions, MSE_cv5s(1:ii, 3), cross_validation_5, solveFor, data);
  
  MSEBagging(ii, 6) = calculateBaggingMSE(predictions, MSE_trainings(1:ii, 3), trainingPosition, solveFor, data);
  LLBagging(ii, 6)  = calculateLogLikelihood(predictions, MSE_trainings(1:ii, 3), trainingPosition, solveFor, data);
  
  MSEBagging(ii, 7) = mean( MSEBagging(ii, 1:5) );
  LLBagging(ii, 7) = mean( LLBagging(ii, 1:5) );
  
  disp([solveFor ": Naive search with bagging for combining " num2str(ii) ...
  " sets: cross-validation MSE: " num2str( MSEBagging(ii, 7), "%2.4f" ) 
  "; log likelihood: " num2str( LLBagging(ii, 7), "%2.4f" )]);
end
[minMSEBagging minMSEBaggingPosition] = min(MSEBagging(:, 7));
[maxLLBagging maxLLBaggingPosition] = max(LLBagging(:, 7));

timeNaiveBagging = time - timei2;
disp([solveFor ": Time: " num2str(timeNaiveBagging, "%2.2f") ...
      " s; Naive search with bagging that minimized cross-validation MSE (" ...
      num2str(minMSEBagging, "%2.4f") ") had cross-validation log likelihood of " ...
      num2str(LLBagging(minMSEBaggingPosition, 7), "%2.4f") " and combined " num2str(minMSEBaggingPosition) " sets. " ...
      "The MSE training when combining the same number of sets was " ...
      num2str( MSEBagging(minMSEBaggingPosition, 6), "%2.4f" ) " and the log likelihood was" ...
      num2str( LLBagging(minMSEBaggingPosition, 6), "%2.4f" )]);


% saving
save([folder_ensembleData "/naiveBagging" solveFor ".mat"], "-V7", "MSEBagging", "minMSEBagging", ...
      "minMSEBaggingPosition", "LLBagging", "maxLLBagging", "maxLLBaggingPosition", "timeNaiveBagging");
% load([folder_ensembleData "/naiveBagging" solveFor ".mat"])
      
% Random search and bagging
%% Defining CDFs
timei3 = time;
scaling = 1/sum(1./MSE_trainings( 1:numberOfNonNanSets, 1) );
CDF_MSE_trainings = NA*zeros(1, size(MSE_trainings, 1) );
PDF_MSE_trainings = NA*zeros(1, size(MSE_trainings, 1) );
PDF_MSE_trainings(1) = scaling/MSE_trainings(1,1);
CDF_MSE_trainings(1) = PDF_MSE_trainings(1);
for ii = 2:numberOfNonNanSets
  PDF_MSE_trainings(ii) = scaling/MSE_trainings(ii,1);
  CDF_MSE_trainings(ii) = PDF_MSE_trainings(ii) + CDF_MSE_trainings(ii - 1);
end

numberOfGroupsToSample = maxSetsToCombine*5;
randomGroupsMSE = cell(numberOfGroupsToSample, 2);
MSERandomBagging = NA*zeros(numberOfGroupsToSample, 7);
LLRandomBagging = NA*zeros(numberOfGroupsToSample, 7); % log likelihood
% will have (1) MSE, (2) numberOfSets, (3) sets positions
for ii = 1:numberOfGroupsToSample
  numberOfSets = unidrnd(maxSetsToCombine - 1) + 1;
  randomGroupsMSE(ii, 1) = numberOfSets;
  
  % getting the sets id
  randomGroupsMSE(ii, 2) = zeros(1, numberOfSets);
  jj = 1;
  while(jj <= numberOfSets)
    % sampling a set
    set_j_position = find(CDF_MSE_trainings >= unifrnd(0, 1), 1);
    if ( isempty( find(randomGroupsMSE{ii, 2} == set_j_position) ) )
      % if my new sample is not equal to any other sample I already had, I add it
      % and move to the next one
      randomGroupsMSE{ii, 2}(jj) = set_j_position;
      jj += 1;
    end
  end
  
  % now I have selected the positions to use and can calculate MSE
  MSERandomBagging(ii, 1) = calculateBaggingMSE(predictions, MSE_cv1s(randomGroupsMSE{ii, 2}, 3), cross_validation_1c, solveFor, data);
  LLRandomBagging(ii, 1) = calculateLogLikelihood(predictions, MSE_cv1s(randomGroupsMSE{ii, 2}, 3), cross_validation_1c, solveFor, data);
  
  MSERandomBagging(ii, 2) = calculateBaggingMSE(predictions, MSE_cv2s(randomGroupsMSE{ii, 2}, 3), cross_validation_2c, solveFor, data);
  LLRandomBagging(ii, 2) = calculateLogLikelihood(predictions, MSE_cv2s(randomGroupsMSE{ii, 2}, 3), cross_validation_2c, solveFor, data);
  
  MSERandomBagging(ii, 3) = calculateBaggingMSE(predictions, MSE_cv3s(randomGroupsMSE{ii, 2}, 3), cross_validation_3c, solveFor, data);
  LLRandomBagging(ii, 3) = calculateLogLikelihood(predictions, MSE_cv3s(randomGroupsMSE{ii, 2}, 3), cross_validation_3c, solveFor, data);
  
  MSERandomBagging(ii, 4) = calculateBaggingMSE(predictions, MSE_cv4s(randomGroupsMSE{ii, 2}, 3), cross_validation_4c, solveFor, data);
  LLRandomBagging(ii, 4) = calculateLogLikelihood(predictions, MSE_cv4s(randomGroupsMSE{ii, 2}, 3), cross_validation_4c, solveFor, data);
  
  MSERandomBagging(ii, 5) = calculateBaggingMSE(predictions, MSE_cv5s(randomGroupsMSE{ii, 2}, 3), cross_validation_5c, solveFor, data);
  LLRandomBagging(ii, 5) = calculateLogLikelihood(predictions, MSE_cv5s(randomGroupsMSE{ii, 2}, 3), cross_validation_5c, solveFor, data);
  
  MSERandomBagging(ii, 6) = calculateBaggingMSE(predictions, MSE_trainings(randomGroupsMSE{ii, 2}, 3), trainingPosition, solveFor, data);
  LLRandomBagging(ii, 6) = calculateLogLikelihood(predictions, MSE_trainings(randomGroupsMSE{ii, 2}, 3), trainingPosition, solveFor, data);
  
  MSERandomBagging(ii, 7) = mean( MSERandomBagging(ii, 1:5) );
  LLRandomBagging(ii, 7) = mean( LLRandomBagging(ii, 1:5) );
  
  disp([solveFor ": Iteration " num2str(ii) ": Random search with bagging for combining " num2str(numberOfSets) ...
        " sets: cross-validation MSE: " num2str( MSERandomBagging(ii, 7), "%2.4f" ) 
        "; log likelihood: " num2str( LLRandomBagging(ii, 7), "%2.4f" )]);
end
[minMSERandomBagging minMSERandomBaggingPosition] = min(MSERandomBagging(:, 7));
[maxLLRandomBagging maxLLRandomBaggingPosition] = max(LLRandomBagging(:, 7));

timeRandomBagging = time - timei3;
disp([solveFor ": Time: " num2str(timeRandomBagging, "%2.2f") ...
      " s; Random search with bagging that minimized cross-validation MSE (" ...
      num2str(minMSERandomBagging, "%2.4f") ") had cross-validation log likelihood of " ...
      num2str(LLRandomBagging(minMSERandomBaggingPosition, 7), "%2.4f") " for the iteration " ...
      num2str(minMSERandomBaggingPosition) ...
      " which combined " num2str( randomGroupsMSE{minMSERandomBaggingPosition, 1} ) " sets. " ...
      "The MSE training for the same iteration was " ...
      num2str( MSERandomBagging(minMSERandomBaggingPosition, 6), "%2.4f" ) " and the log likelihood was" ...
      num2str( LLRandomBagging(minMSERandomBaggingPosition, 6), "%2.4f" )]);
      
% saving
save([folder_ensembleData "/randomBagging" solveFor ".mat"], "-V7", "MSERandomBagging", "randomGroupsMSE", ...
      "minMSERandomBagging", "minMSERandomBaggingPosition", "LLRandomBagging", "maxLLRandomBagging", ...
      "maxLLRandomBaggingPosition", "CDF_MSE_trainings", ...
      "PDF_MSE_trainings", "timeRandomBagging", "numberOfGroupsToSample");
% load([folder_ensembleData "/randomBagging" solveFor ".mat"])
      
      
      
% Greedy search and bagging
timei4 = time;
MSEGreedyBagging = NA*zeros(maxSetsToCombine, 13, length(lambdaVector));
setsGreedyBagging = NA*zeros(maxSetsToCombine, 6, length(lambdaVector));
LLGreedyBagging = NA*zeros(maxSetsToCombine, 7, length(lambdaVector));
ObjGreedyBagging = NA*zeros(maxSetsToCombine, 6, length(lambdaVector));

%% preparing to run parallel
% call for parallel function
greedyBaggingParallel = @(training, testing, setsToSearch, msg, lambda) greedyBagging(predictions, ... 
                          maxSetsToCombine, training, testing, setsToSearch, solveFor, data, msg, lambda);
                          
trainingCell = cell(6,1); testingCell = cell(6,1); msg = cell(6,1); 
setsToSearchCell = cell(6,1); lambdaVectorCell = cell(6,1);

trainingCell{1} = cross_validation_1; trainingCell{2} = cross_validation_2;
trainingCell{3} = cross_validation_3; trainingCell{4} = cross_validation_4;
trainingCell{5} = cross_validation_5; trainingCell{6} = trainingPosition;

testingCell{1} = cross_validation_1c; testingCell{2} = cross_validation_2c;
testingCell{3} = cross_validation_3c; testingCell{4} = cross_validation_4c;
testingCell{5} = cross_validation_5c; testingCell{6} = trainingPosition;

setsToSearchCell{1} = MSE_cv1s(1:numberOfNonNanSets, 3); setsToSearchCell{2} = MSE_cv2s(1:numberOfNonNanSets, 3);
setsToSearchCell{3} = MSE_cv3s(1:numberOfNonNanSets, 3); setsToSearchCell{4} = MSE_cv4s(1:numberOfNonNanSets, 3);
setsToSearchCell{5} = MSE_cv5s(1:numberOfNonNanSets, 3); setsToSearchCell{6} = MSE_trainings(1:numberOfNonNanSets, 3);

msg{1} = [solveFor ": CV1: "]; msg{2} = [solveFor ": CV2: "]; msg{3} = [solveFor ": CV3: "];
msg{4} = [solveFor ": CV4: "]; msg{5} = [solveFor ": CV5: "]; msg{6} = [solveFor ": Training: "];


for ii = 1:length(lambdaVector)
  lambdaVectorCell{1} = lambdaVector(ii); lambdaVectorCell{2} = lambdaVector(ii);
  lambdaVectorCell{3} = lambdaVector(ii); lambdaVectorCell{4} = lambdaVector(ii);
  lambdaVectorCell{5} = lambdaVector(ii); lambdaVectorCell{6} = lambdaVector(ii);
  [MSEPGBout SPGBout LLPGBout OFPGBout] = parcellfun(6, greedyBaggingParallel, trainingCell, ...
                                           testingCell, setsToSearchCell, msg, lambdaVectorCell);
  
  MSEGreedyBagging(:, 1:2, ii)   = MSEPGBout( 1:maxSetsToCombine, :);
  MSEGreedyBagging(:, 3:4, ii)   = MSEPGBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  MSEGreedyBagging(:, 5:6, ii)   = MSEPGBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  MSEGreedyBagging(:, 7:8, ii)   = MSEPGBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  MSEGreedyBagging(:, 9:10, ii)  = MSEPGBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  MSEGreedyBagging(:, 11:12, ii) = MSEPGBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
  MSEGreedyBagging(:, 13, ii)    = mean( MSEGreedyBagging(:, [2:2:10], ii), 2);
  
  % the algorithm outputs set_ID. -1 converts set_ID to #set
  setsGreedyBagging(:, 1, ii) = SPGBout( 1:maxSetsToCombine ) - 1;
  setsGreedyBagging(:, 2, ii) = SPGBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :) - 1;
  setsGreedyBagging(:, 3, ii) = SPGBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :) - 1;
  setsGreedyBagging(:, 4, ii) = SPGBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :) - 1;
  setsGreedyBagging(:, 5, ii) = SPGBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :) - 1;
  setsGreedyBagging(:, 6, ii) = SPGBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :) - 1;
  
  LLGreedyBagging(:, 1, ii) = LLPGBout( 1:maxSetsToCombine );
  LLGreedyBagging(:, 2, ii) = LLPGBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  LLGreedyBagging(:, 3, ii) = LLPGBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  LLGreedyBagging(:, 4, ii) = LLPGBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  LLGreedyBagging(:, 5, ii) = LLPGBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  LLGreedyBagging(:, 6, ii) = LLPGBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
  LLGreedyBagging(:, 7, ii) = mean( LLGreedyBagging(:, 1:5, ii), 2);
  
  ObjGreedyBagging(:, 1, ii) = OFPGBout( 1:maxSetsToCombine );
  ObjGreedyBagging(:, 2, ii) = OFPGBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  ObjGreedyBagging(:, 3, ii) = OFPGBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  ObjGreedyBagging(:, 4, ii) = OFPGBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  ObjGreedyBagging(:, 5, ii) = OFPGBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  ObjGreedyBagging(:, 6, ii) = OFPGBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
end
  [minMSEGreedyBagging minMSEGreedyBaggingPosition] = min(MSEGreedyBagging(:, 13, 1));

timeGreedyBagging = time - timei4;
disp([solveFor ": Time: " num2str(timeGreedyBagging, "%2.2f") " s; Greedy search with bagging had minimum cross-validation MSE (" ...
      num2str(minMSEGreedyBagging, "%2.4f") ") when combining " num2str(minMSEGreedyBaggingPosition) " sets. " ...
      "The MSE training when combining the same number of sets was " ...
      num2str( MSEGreedyBagging(minMSEGreedyBaggingPosition, 11, 1), "%2.4f" )]);
% saving
save([folder_ensembleData "/greedyBagging" solveFor ".mat"], "-V7", "MSEGreedyBagging", "setsGreedyBagging", ...
     "LLGreedyBagging", "ObjGreedyBagging", "minMSEGreedyBagging", "minMSEGreedyBaggingPosition", ...
     "lambdaVector", "timeGreedyBagging", "trainingCell", "testingCell", "msg", ...
     "setsToSearchCell", "lambdaVectorCell");
% load([folder_ensembleData "/greedyBagging" solveFor ".mat"])

     
     
     
% Naive search and optimum weights
timei5 = time;
MSEOptimumWeights = cell(7, length(lambdaVector));
setsOptimumWeights = cell(6, length(lambdaVector));
LLOptimumWeights = cell(7, length(lambdaVector));
ObjOptimumWeights = cell(7, length(lambdaVector));
naiveOptimumWeights = cell(6, length(lambdaVector));

setsToSearchNOWCell = cell(6,1);
setsToSearchNOWCell{1} = MSE_cv1s(1:maxSetsToCombine, 3); setsToSearchNOWCell{2} = MSE_cv2s(1:maxSetsToCombine, 3);
setsToSearchNOWCell{3} = MSE_cv3s(1:maxSetsToCombine, 3); setsToSearchNOWCell{4} = MSE_cv4s(1:maxSetsToCombine, 3);
setsToSearchNOWCell{5} = MSE_cv5s(1:maxSetsToCombine, 3); setsToSearchNOWCell{6} = MSE_trainings(1:maxSetsToCombine, 3);

naiveOptimumWeightParallel = @(training, testing, setsToSearch, msg, lambda) naiveOptimumWeight(predictions, ... 
                                               training, testing, setsToSearch, solveFor, data, msg, lambda, ...
                                               numberOfInitialSetsOptimumWeights);


gradientObjectiveFunctionOptimumWeightMSE(0, 0, 0, 0, 0, 1);
HessianObjectiveFunctionOptimumWeightMSE(0, 0, 0, 0, 0, 1);

for ii = 1:length(lambdaVector)
  lambdaVectorCell{1} = lambdaVector(ii); lambdaVectorCell{2} = lambdaVector(ii);
  lambdaVectorCell{3} = lambdaVector(ii); lambdaVectorCell{4} = lambdaVector(ii);
  lambdaVectorCell{5} = lambdaVector(ii); lambdaVectorCell{6} = lambdaVector(ii);

  [MSEPNOWBout SPNOWBout NWPOWBout LLNGOWBout OFNGOWBout] = parcellfun(6, naiveOptimumWeightParallel, trainingCell, ...
                                                        testingCell, setsToSearchNOWCell, msg, lambdaVectorCell);
                                                                              
  MSEOptimumWeights{1, ii} = MSEPNOWBout( 1:maxSetsToCombine, :);
  MSEOptimumWeights{2, ii} = MSEPNOWBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  MSEOptimumWeights{3, ii} = MSEPNOWBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  MSEOptimumWeights{4, ii} = MSEPNOWBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  MSEOptimumWeights{5, ii} = MSEPNOWBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  MSEOptimumWeights{6, ii} = MSEPNOWBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
  MSEOptimumWeights{7, ii} = mean( [ MSEOptimumWeights{1, ii}(:,2) ...
                                     MSEOptimumWeights{2, ii}(:,2) ...
                                     MSEOptimumWeights{3, ii}(:,2) ...
                                     MSEOptimumWeights{4, ii}(:,2) ...
                                     MSEOptimumWeights{5, ii}(:,2) ], 2);
  
  setsOptimumWeights{1, ii} = SPNOWBout(1, :);
  setsOptimumWeights{2, ii} = SPNOWBout(2, :);
  setsOptimumWeights{3, ii} = SPNOWBout(3, :);
  setsOptimumWeights{4, ii} = SPNOWBout(4, :);
  setsOptimumWeights{5, ii} = SPNOWBout(5, :);
  setsOptimumWeights{6, ii} = SPNOWBout(6, :);
  
  naiveOptimumWeights{1, ii} = NWPOWBout(1, :);
  naiveOptimumWeights{2, ii} = NWPOWBout(2, :);
  naiveOptimumWeights{3, ii} = NWPOWBout(3, :);
  naiveOptimumWeights{4, ii} = NWPOWBout(4, :);
  naiveOptimumWeights{5, ii} = NWPOWBout(5, :);
  naiveOptimumWeights{6, ii} = NWPOWBout(6, :);
  
  LLOptimumWeights{1, ii} = LLNGOWBout( 1:maxSetsToCombine );
  LLOptimumWeights{2, ii} = LLNGOWBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  LLOptimumWeights{3, ii} = LLNGOWBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  LLOptimumWeights{4, ii} = LLNGOWBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  LLOptimumWeights{5, ii} = LLNGOWBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  LLOptimumWeights{6, ii} = LLNGOWBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
  LLOptimumWeights{7, ii} = mean( [ LLOptimumWeights{1, ii} ...
                                    LLOptimumWeights{2, ii} ...
                                    LLOptimumWeights{3, ii} ...
                                    LLOptimumWeights{4, ii} ...
                                    LLOptimumWeights{5, ii} ], 2);
  
  ObjOptimumWeights{1, ii} = OFNGOWBout( 1:maxSetsToCombine );
  ObjOptimumWeights{2, ii} = OFNGOWBout( (maxSetsToCombine + 1):(2*maxSetsToCombine), :);
  ObjOptimumWeights{3, ii} = OFNGOWBout( (2*maxSetsToCombine + 1):(3*maxSetsToCombine), :);
  ObjOptimumWeights{4, ii} = OFNGOWBout( (3*maxSetsToCombine + 1):(4*maxSetsToCombine), :);
  ObjOptimumWeights{5, ii} = OFNGOWBout( (4*maxSetsToCombine + 1):(5*maxSetsToCombine), :);
  ObjOptimumWeights{6, ii} = OFNGOWBout( (5*maxSetsToCombine + 1):(6*maxSetsToCombine), :);
  ObjOptimumWeights{7, ii} = mean( [ ObjOptimumWeights{1, ii} ...
                                     ObjOptimumWeights{2, ii} ...
                                     ObjOptimumWeights{3, ii} ...
                                     ObjOptimumWeights{4, ii} ...
                                     ObjOptimumWeights{5, ii} ], 2);
                                           
end
[minMSEOptimumWeights minMSEOptimumWeightsPosition] = min(MSEOptimumWeights{7, 31});

timeNaiveOptimumWeight = time - timei5;
disp([solveFor "Time: " num2str(timeNaiveOptimumWeight, "%2.2f") " s; Naive search with optimum weights had minimum cross-validation MSE (" ...
      num2str(minMSEOptimumWeights, "%2.4f") ") when combining " num2str(minMSEOptimumWeightsPosition) " sets. " ...
      "The MSE training when combining the same number of sets was " ...
      num2str( MSEOptimumWeights(minMSEOptimumWeightsPosition, 11, 1), "%2.4f" ) ...
      " and it had " num2str( sum( optimumWeights{minMSEOptimumWeightsPosition, 6, 1} > 0 ) ) ...
      " non-zero weights."]);

% Adapting naiveOptimumWeights for compatibility with R
for ii = 1:length(lambdaVector)
  for jj = 1:19
    for kk = 1:6
      naiveOptimumWeights{kk, ii}{jj} = NA;
    end
  end
end


% saving
save([folder_ensembleData "/naiveOptimumWeight" solveFor ".mat"], "-V7", "MSEOptimumWeights", "naiveOptimumWeights", ...
     "minMSEOptimumWeights", "minMSEOptimumWeightsPosition", "setsOptimumWeights", "LLOptimumWeights", ...
     "ObjOptimumWeights", "timeNaiveOptimumWeight", "setsToSearchNOWCell");
% load([folder_ensembleData "/naiveOptimumWeight" solveFor ".mat"])    

     
     
     
     
% Random search with optimum weight
% Will selected between min and max sets to add
timei6 = time;
MSERandomOptimumWeights = cell(7, length(lambdaVector));
randomOptimumWeights = cell(6, length(lambdaVector));
LLRandomOptimumWeights = cell(7, length(lambdaVector));
ObjRandomOptimumWeights = cell(7, length(lambdaVector));
randomGroupsMSEOptimumWeight = cell(numberOfGroupsToSample, 2);

% will have (1) MSE, (2) numberOfSets, (3) sets positions

for ii = 1:numberOfGroupsToSample
  numberOfSetsRandom = unidrnd(maxSetsToRandomlyCombineOptimumWeight ...
                         - minSetsToRandomlyCombineOptimumWeight + 1) ...
                         + minSetsToRandomlyCombineOptimumWeight - 1;
                         
  randomGroupsMSEOptimumWeight(ii, 1) = numberOfSetsRandom;
  
  % getting the sets id
  randomGroupsMSEOptimumWeight(ii, 2) = zeros(1, numberOfSetsRandom);
  jj = 1;
  while(jj <= numberOfSetsRandom)
    % sampling a set
    set_j_position = find(CDF_MSE_trainings >= unifrnd(0, 1), 1);
    if ( isempty( find(randomGroupsMSEOptimumWeight{ii, 2} == set_j_position) ) )
      % if my new sample is not equal to any other sample I already had, I add it
      % and move to the next one
      randomGroupsMSEOptimumWeight{ii, 2}(jj) = set_j_position;
      jj += 1;
    end
  end
end

randomOptimumWeightParallel = @(training, testing, setsToSearch, msg, lambda) randomOptimumWeight(predictions, ... 
                          training, testing, randomGroupsMSEOptimumWeight, setsToSearch, solveFor, data, msg, lambda);
                          
for ii = 1:length(lambdaVector)
  lambdaVectorCell{1} = lambdaVector(ii); lambdaVectorCell{2} = lambdaVector(ii);
  lambdaVectorCell{3} = lambdaVector(ii); lambdaVectorCell{4} = lambdaVector(ii);
  lambdaVectorCell{5} = lambdaVector(ii); lambdaVectorCell{6} = lambdaVector(ii);
  [MSEPROWBout RWPOWBout LLPROWBout OFPROWBout] = parcellfun(6, randomOptimumWeightParallel, trainingCell, ...
                                                        testingCell, setsToSearchCell, msg, lambdaVectorCell);
                                                                              
  MSERandomOptimumWeights{1, ii} = MSEPROWBout( 1:numberOfGroupsToSample, :);
  MSERandomOptimumWeights{2, ii} = MSEPROWBout( (numberOfGroupsToSample + 1):(2*numberOfGroupsToSample), :);
  MSERandomOptimumWeights{3, ii} = MSEPROWBout( (2*numberOfGroupsToSample + 1):(3*numberOfGroupsToSample), :);
  MSERandomOptimumWeights{4, ii} = MSEPROWBout( (3*numberOfGroupsToSample + 1):(4*numberOfGroupsToSample), :);
  MSERandomOptimumWeights{5, ii} = MSEPROWBout( (4*numberOfGroupsToSample + 1):(5*numberOfGroupsToSample), :);
  MSERandomOptimumWeights{6, ii} = MSEPROWBout( (5*numberOfGroupsToSample + 1):(6*numberOfGroupsToSample), :);
  MSERandomOptimumWeights{7, ii} = mean( [ MSERandomOptimumWeights{1, ii}(:,2) ...
                                           MSERandomOptimumWeights{2, ii}(:,2) ...
                                           MSERandomOptimumWeights{3, ii}(:,2) ...
                                           MSERandomOptimumWeights{4, ii}(:,2) ...
                                           MSERandomOptimumWeights{5, ii}(:,2) ], 2);
  
  randomOptimumWeights{1, ii} = RWPOWBout(1, :);
  randomOptimumWeights{2, ii} = RWPOWBout(2, :);
  randomOptimumWeights{3, ii} = RWPOWBout(3, :);
  randomOptimumWeights{4, ii} = RWPOWBout(4, :);
  randomOptimumWeights{5, ii} = RWPOWBout(5, :);
  randomOptimumWeights{6, ii} = RWPOWBout(6, :);
  
  LLRandomOptimumWeights{1, ii} = LLPROWBout( 1:numberOfGroupsToSample );
  LLRandomOptimumWeights{2, ii} = LLPROWBout( (numberOfGroupsToSample + 1):(2*numberOfGroupsToSample), :);
  LLRandomOptimumWeights{3, ii} = LLPROWBout( (2*numberOfGroupsToSample + 1):(3*numberOfGroupsToSample), :);
  LLRandomOptimumWeights{4, ii} = LLPROWBout( (3*numberOfGroupsToSample + 1):(4*numberOfGroupsToSample), :);
  LLRandomOptimumWeights{5, ii} = LLPROWBout( (4*numberOfGroupsToSample + 1):(5*numberOfGroupsToSample), :);
  LLRandomOptimumWeights{6, ii} = LLPROWBout( (5*numberOfGroupsToSample + 1):(6*numberOfGroupsToSample), :);
  LLRandomOptimumWeights{7, ii} = mean( [ LLRandomOptimumWeights{1, ii} ...
                                          LLRandomOptimumWeights{2, ii} ...
                                          LLRandomOptimumWeights{3, ii} ...
                                          LLRandomOptimumWeights{4, ii} ...
                                          LLRandomOptimumWeights{5, ii} ], 2);
  
  ObjRandomOptimumWeights{1, ii} = OFPROWBout( 1:numberOfGroupsToSample );
  ObjRandomOptimumWeights{2, ii} = OFPROWBout( (numberOfGroupsToSample + 1):(2*numberOfGroupsToSample), :);
  ObjRandomOptimumWeights{3, ii} = OFPROWBout( (2*numberOfGroupsToSample + 1):(3*numberOfGroupsToSample), :);
  ObjRandomOptimumWeights{4, ii} = OFPROWBout( (3*numberOfGroupsToSample + 1):(4*numberOfGroupsToSample), :);
  ObjRandomOptimumWeights{5, ii} = OFPROWBout( (4*numberOfGroupsToSample + 1):(5*numberOfGroupsToSample), :);
  ObjRandomOptimumWeights{6, ii} = OFPROWBout( (5*numberOfGroupsToSample + 1):(6*numberOfGroupsToSample), :);
  ObjRandomOptimumWeights{7, ii} = mean( [ ObjRandomOptimumWeights{1, ii} ...
                                           ObjRandomOptimumWeights{2, ii} ...
                                           ObjRandomOptimumWeights{3, ii} ...
                                           ObjRandomOptimumWeights{4, ii} ...
                                           ObjRandomOptimumWeights{5, ii} ], 2);

end

timeRandomOptimumWeight = time - timei6;
disp([solveFor ": Time: " num2str(timeGreedyOptimumWeight, "%2.2f") " s; Random selection with optimum weights is finished."]);
% saving
save([folder_ensembleData "/randomOptimumWeight" solveFor ".mat"], "-V7", ...
     "randomGroupsMSEOptimumWeight", "MSERandomOptimumWeights", ...
     "randomOptimumWeights", "LLRandomOptimumWeights", "ObjRandomOptimumWeights", "timeRandomOptimumWeight");
% load([folder_ensembleData "/randomOptimumWeight" solveFor ".mat"])

                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
%% Greedy search with optimum weight
timei7 = time;
MSEGreedyOptimumWeight = cell(7, length(lambdaVector));
setsGreedyOptimumWeights = cell(6, length(lambdaVector));
LLGreedyOptimumWeights = cell(7, length(lambdaVector));
ObjGreedyOptimumWeights = cell(7, length(lambdaVector));
greedyWeights = cell(6, length(lambdaVector));


%% preparing to run parallel
% call for parallel function
greedyOptimumWeightParallel = @(training, testing, setsToSearch, msg, lambda) greedyOptimumWeight(predictions, ... 
                          training, testing, setsToSearch, solveFor, data, msg, lambda, numberOfInitialSetsOptimumWeights);

for ii = 1:length(lambdaVector)
  lambdaVectorCell{1} = lambdaVector(ii); lambdaVectorCell{2} = lambdaVector(ii);
  lambdaVectorCell{3} = lambdaVector(ii); lambdaVectorCell{4} = lambdaVector(ii);
  lambdaVectorCell{5} = lambdaVector(ii); lambdaVectorCell{6} = lambdaVector(ii);
  [MSEPGOWBout SPGOWBout GWPOWBout LLPGOWBout OFPGOWBout] = parcellfun(6, greedyOptimumWeightParallel, trainingCell, ...
                                                        testingCell, setsToSearchCell, msg, lambdaVectorCell);
                                                                              
  MSEGreedyOptimumWeight{1, ii} = MSEPGOWBout( 1:numberOfNonNanSets, :);
  MSEGreedyOptimumWeight{2, ii} = MSEPGOWBout( (numberOfNonNanSets + 1):(2*numberOfNonNanSets), :);
  MSEGreedyOptimumWeight{3, ii} = MSEPGOWBout( (2*numberOfNonNanSets + 1):(3*numberOfNonNanSets), :);
  MSEGreedyOptimumWeight{4, ii} = MSEPGOWBout( (3*numberOfNonNanSets + 1):(4*numberOfNonNanSets), :);
  MSEGreedyOptimumWeight{5, ii} = MSEPGOWBout( (4*numberOfNonNanSets + 1):(5*numberOfNonNanSets), :);
  MSEGreedyOptimumWeight{6, ii} = MSEPGOWBout( (5*numberOfNonNanSets + 1):(6*numberOfNonNanSets), :);
  MSEGreedyOptimumWeight{7, ii} = mean( [ MSEGreedyOptimumWeight{1, ii}(:,2) ...
                                          MSEGreedyOptimumWeight{2, ii}(:,2) ...
                                          MSEGreedyOptimumWeight{3, ii}(:,2) ...
                                          MSEGreedyOptimumWeight{4, ii}(:,2) ...
                                          MSEGreedyOptimumWeight{5, ii}(:,2) ], 2);
  
  setsGreedyOptimumWeights{1, ii} = SPGOWBout(1, :);
  setsGreedyOptimumWeights{2, ii} = SPGOWBout(2, :);
  setsGreedyOptimumWeights{3, ii} = SPGOWBout(3, :);
  setsGreedyOptimumWeights{4, ii} = SPGOWBout(4, :);
  setsGreedyOptimumWeights{5, ii} = SPGOWBout(5, :);
  setsGreedyOptimumWeights{6, ii} = SPGOWBout(6, :);
  
  greedyWeights{1, ii} = GWPOWBout(1, :);
  greedyWeights{2, ii} = GWPOWBout(2, :);
  greedyWeights{3, ii} = GWPOWBout(3, :);
  greedyWeights{4, ii} = GWPOWBout(4, :);
  greedyWeights{5, ii} = GWPOWBout(5, :);
  greedyWeights{6, ii} = GWPOWBout(6, :);
  
  LLGreedyOptimumWeights{1, ii} = LLPGOWBout( 1:numberOfNonNanSets );
  LLGreedyOptimumWeights{2, ii} = LLPGOWBout( (numberOfNonNanSets + 1):(2*numberOfNonNanSets), :);
  LLGreedyOptimumWeights{3, ii} = LLPGOWBout( (2*numberOfNonNanSets + 1):(3*numberOfNonNanSets), :);
  LLGreedyOptimumWeights{4, ii} = LLPGOWBout( (3*numberOfNonNanSets + 1):(4*numberOfNonNanSets), :);
  LLGreedyOptimumWeights{5, ii} = LLPGOWBout( (4*numberOfNonNanSets + 1):(5*numberOfNonNanSets), :);
  LLGreedyOptimumWeights{6, ii} = LLPGOWBout( (5*numberOfNonNanSets + 1):(6*numberOfNonNanSets), :);
  LLGreedyOptimumWeights{7, ii} = mean( [ LLGreedyOptimumWeights{1, ii} ...
                                          LLGreedyOptimumWeights{2, ii} ...
                                          LLGreedyOptimumWeights{3, ii} ...
                                          LLGreedyOptimumWeights{4, ii} ...
                                          LLGreedyOptimumWeights{5, ii} ], 2);
  
  ObjGreedyOptimumWeights{1, ii} = OFPGOWBout( 1:numberOfNonNanSets );
  ObjGreedyOptimumWeights{2, ii} = OFPGOWBout( (numberOfNonNanSets + 1):(2*numberOfNonNanSets), :);
  ObjGreedyOptimumWeights{3, ii} = OFPGOWBout( (2*numberOfNonNanSets + 1):(3*numberOfNonNanSets), :);
  ObjGreedyOptimumWeights{4, ii} = OFPGOWBout( (3*numberOfNonNanSets + 1):(4*numberOfNonNanSets), :);
  ObjGreedyOptimumWeights{5, ii} = OFPGOWBout( (4*numberOfNonNanSets + 1):(5*numberOfNonNanSets), :);
  ObjGreedyOptimumWeights{6, ii} = OFPGOWBout( (5*numberOfNonNanSets + 1):(6*numberOfNonNanSets), :);
  ObjGreedyOptimumWeights{7, ii} = mean( [ ObjGreedyOptimumWeights{1, ii} ...
                                           ObjGreedyOptimumWeights{2, ii} ...
                                           ObjGreedyOptimumWeights{3, ii} ...
                                           ObjGreedyOptimumWeights{4, ii} ...
                                           ObjGreedyOptimumWeights{5, ii} ], 2);

end
[minMSEGreedyOptimumWeight minMSEGreedyOptimumWeightPosition] = min(MSEGreedyOptimumWeight{7, 1});

% Adapting naiveOptimumWeights for compatibility with R
for ii = 1:length(lambdaVector)
  for jj = 1:19
    for kk = 1:6
      greedyWeights{kk, ii}{jj} = NA;
    end
  end
end

timeGreedyOptimumWeight = time - timei7;
disp([solveFor ": Time: " num2str(timeGreedyOptimumWeight, "%2.2f") " s; Greedy selection with optimum weights had minimum cross-validation MSE (" ...
      num2str(minMSEGreedyOptimumWeight, "%2.4f") ") in the iteration " num2str(minMSEGreedyOptimumWeightPosition) ...
      ", which combined " num2str( length( greedyWeights{6, 1}{minMSEGreedyOptimumWeightPosition} ) ) " training sets and " ...
      "had MSE training of " ...
      num2str( MSEGreedyOptimumWeight{6, 1}(minMSEGreedyOptimumWeightPosition, 2), "%2.4f" ) ...
      "."]);
% saving
save([folder_ensembleData "/greedyOptimumWeight" solveFor ".mat"], "-V7", "MSEGreedyOptimumWeight", "setsGreedyOptimumWeights", ...
     "greedyWeights", "minMSEGreedyOptimumWeight", "minMSEGreedyOptimumWeightPosition", "timeGreedyOptimumWeight", ...
     "LLGreedyOptimumWeights", "ObjGreedyOptimumWeights");
% load([folder_ensembleData "/greedyOptimumWeight" solveFor ".mat"])


end
