function predictionsEnsemble = ensemblePrediction(sets, weights, training = NA, ...
                                 solveFor = "all", heat_power_in = NA, ...
                                 Ta_pen_in = NA, Ta_brooder_in = NA, ...
                                 Tg_brooder_in = NA, Tr_in = NA, ...
                                 nameEnsemble = "", nparallel = 3, ...
                                 loadEstimations = 1, Ta_brooder_in_increase = 0)
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
% File:   ensemblePrediction.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
% Created on May 19, 2018.
%
%
% Function description:
% Calculates the MSE of the given positions of the specified predictions
%
%
%
% Usage:
% predictionsEnsemble = ensemblePrediction(sets, weights, training, ...
%                                  solveFor, heat_power_in, ...
%                                  Ta_pen_in, Ta_brooder_in, ...
%                                  Tg_brooder_in, Tr_in, nameEnsemble, nparallel, ...
%                                  loadEstimations, Ta_brooder_in_increase)
%
% Input:
% sets: a vector with the #sets.
% weights: a vector with the weigths for every #sets.
% training: Should calculate MSE for training, testing, or both? If 0, calculates 
%           MSE for testing. If 1, calculates MSE for training. If 2, calculates MSE for
%           both, testing and training. Otherwise, will not load any.
% solveFor: Specify if should solve for all models or for a particular model.
%           Valid inputs are:
%            (1) "all": loads solutions from all models.
%            (2) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (3) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
%
% heat_power_in: Power of the supplemental heat source (W). Can be a vector. Used to
%                estimate Ta_pen, Ta_brooder, Tg_brooder, and Tr. If NA, will check
%                if positions is NA. If so, will throw an error. If not, will load the
%                measured dataset and solve the problem for the given positions.
% Ta_pen_in: Air temperature inside the pen (oC). Can be a vector. Used to estimate Ta_brooder. 
%            If NA, will check if Ta_brooder is given. If so, will use the given
%            value of Ta_brooder instead of predicting it. If not, will check if
%            measured dataset should be loaded and solved (tested if heat_power == NA).
%            If not, will throw an error.
% Ta_brooder_in: Air temperature inside the brooder (oC). Can be a vector. Used to
%                estimate Tg_brooder and Tr. If NA, Ta_brooder is predicted given 
%                Ta_pen and heat_power. If not, the given values will be used. In such case, 
%                Ta_brooder will be added by Ta_pen_offset, which is the offset to 
%                compensate for measurement errors.
% Tg_brooder_in: Black globe temperature inside the brooder (oC). Can be a vector.
%                Used to estimate Tr. If NA, Tg_brooder is predicted given Ta_brooder and heat_power.
%                If not, the given values will be used. In such case, 
%                Tg_brooder will be added by Ta_pen_offset, which is the offset to 
%                compensate for measurement errors.
% Tr_in: Rectal temperature (oC). Can be a vector. If NA, Tr is predicted given
%        Ta_brooder, Tg_brooder, and heat_power. If not, the given values will be used.
%        In such case, Tr will be added by Ta_pen_offset, which is the offset to 
%        compensate for measured errors.
% nameEnsemble: name used to salve the ensemble predictions
% nparallel: Number or process to lunch in parallel.
% loadEstimations: when running for training or testing:
%                       0: recalculates the predictions even if they are available
%                       1: load predictions if they are available
% Ta_brooder_in_increase: when Ta_brooder is an input, should it be add by the 
%                         the effect of heat_source_power?
%                           0: Ta_brooder is as the input.
%                           1: Ta_brooder is the input + the added effect of
%                              heat_soure_power in increasing Ta_brooder
%                  
% 
% Output:
% predictionsEnsemble: matrix with the predictions from the model. Rows are the different
%                      inputs and columns are the different variables predicted. Columns are:
%                       1) skin-surface temperature (Ts) weighted by the input "weights"
%                       2) skin-surface heat flux (q_skin) weighted by the input "weights"
%                       3) hair-coat surface temperature (Th) weighted by the input "weights"
%                       4) time it took to solve for this position summed for all calculations
%                       5) value of the objective function. Weighted by the input "weights"
%                       6) maximum value of the constraint function. Weighted by the input "weights"
%                       7) Ta_pen used in the calculation. Weighted by the input "weights"
%                       8) Ta_brooder used in the calculation. Weighted by the input "weights"
%                       9) Tg_brooder used in the calculation. Weighted by the input "weights"
%                      10) Tr used in the calculation. Weighted by the input "weights"
%                      11) skin-surface temperature (Ts) standard deviation. Assumes Ts follows
%                          a normal distribution.
%                      12) skin-surface heat flux (q_skin) standard deviation. Assumes q_skin follows
%                          a normal distribution.
%                      13) hair-coat surface temperature (Th) standard deviation. Assumes Th follows
%                          a normal distribution.


% loading training and testing positions
folder_datasets = locateFolderPath("datasets");
if (training == 1)
  % we will run the model for the training positions
  load([folder_datasets "/trainingPosition.mat"])
  positions = trainingPosition;
  datasetTypeName = "Training";
elseif (training == 0)
  % we will run the model for the testing positions
  load([folder_datasets "/testingPosition.mat"])
  positions = testingPosition;
  datasetTypeName = "Testing";
elseif (training == 2)
% will run for both, testing and training.
  load([folder_datasets "/trainingPosition.mat"])
  load([folder_datasets "/testingPosition.mat"])
  positions = [trainingPosition' testingPosition];
  datasetTypeName = "TestingAndTraining";
else
  positions = NA;
  datasetTypeName = "";
end


% solving for everyone?
if ( strcmp(solveFor, "all") )
  predictionsEnsemble1 = ensemblePrediction(sets, weights, training, ...
                                 "Skin", heat_power_in, ...
                                 Ta_pen_in, Ta_brooder_in, ...
                                 Tg_brooder_in, Tr_in);
  predictionsEnsemble2 = ensemblePrediction(sets, weights, training, ...
                                 "KM", heat_power_in, ...
                                 Ta_pen_in, Ta_brooder_in, ...
                                 Tg_brooder_in, Tr_in);
                                 
  predictionsEnsemble = zeros(size(predictionsEnsemble1, 1), ...
                              size(predictionsEnsemble1, 2), 2);
                              
  predictionsEnsemble(:, :, 1) = predictionsEnsemble1;
  predictionsEnsemble(:, :, 2) = predictionsEnsemble2;
  return
elseif ( !strcmp(solveFor, "Skin") && !strcmp(solveFor, "KM"))
  error(["Wrong input for solveFor (" solveFor "). Accepted inputs are "...
         "\"Skin\" or \"KM\"."]);
end

%% Adding the necessary paths
folder_analyticalModels_KM = locateFolderPath("analyticalModels/KM");
addpath(folder_analyticalModels_KM);
%% Adding the necessary paths
folder_analyticalModels_skin = locateFolderPath("analyticalModels/Skin");
addpath(folder_analyticalModels_skin);

% getting the weights that are non-zero
setsNonZero = NA;
weightsNonZero = NA;
numberOfNonZeros = 0;
for ii = 1:length(weights)
  if (weights(ii) != 0)
    numberOfNonZeros += 1;
    setsNonZero(numberOfNonZeros) = sets(ii);
    weightsNonZero(numberOfNonZeros) = weights(ii);
  end
end

% Running in parallel
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  pkg load parallel
else
  % we are executing in MATLAB
end

if ( strcmp(solveFor, "KM") )
  solveSetParallel = @(sets) solveSetKM(sets, heat_power_in, Ta_pen_in, Ta_brooder_in, ...
                          Tg_brooder_in, Tr_in, positions, datasetTypeName, 1, 0, 0, ...
                          loadEstimations, Ta_brooder_in_increase);
elseif ( strcmp(solveFor, "Skin") )
  solveSetParallel = @(sets) solveSetSkin(sets, heat_power_in, Ta_pen_in, Ta_brooder_in, ...
                          Tg_brooder_in, Tr_in, positions, datasetTypeName, 1, 0, 0, ...
                          loadEstimations, Ta_brooder_in_increase);
end
timei = time;
predictionsEnsembleParallel = pararrayfun(nparallel, solveSetParallel, setsNonZero);

predictionsEnsemble = zeros(size(predictionsEnsembleParallel, 1), 13);
for ii = 1:length(setsNonZero)
  % calculating predictions and variables that are weighted
  predictionsEnsemble(:, [1:3, 5:10] ) += weightsNonZero(ii)*...
                      predictionsEnsembleParallel(:, [1:3, 5:10] + (ii - 1)*10);
  % I will calculated an "average running time" for each point.
  predictionsEnsemble(:, 4) += predictionsEnsembleParallel(:, 4 + (ii - 1)*10); % time is summed
  
  % Variances
  predictionsEnsemble(:, 11:13) += weightsNonZero(ii)*...
                      predictionsEnsembleParallel(:, [1:3] + (ii - 1)*10).^2;
end

predictionsEnsemble(:, 4) *= 1/nparallel;
totalMachineTime = nparallel*sum(predictionsEnsemble(:, 4));

% if we are loading, we will approximate the time. Otherwise, we will use the given time
if (loadEstimations == 1)
  totalTime = totalMachineTime;
else
  totalTime = timei - time;
end
% Finishing computing variances and converting it to standard deviation
predictionsEnsemble(:, 11:13) -= predictionsEnsemble(:, 1:3).^2;
predictionsEnsemble(:, 11:13) = sqrt(predictionsEnsemble(:, 11:13));

% cleaning after myself
rmpath(folder_analyticalModels_KM)
rmpath(folder_analyticalModels_skin)

%% If we are running for training or testing, we will calculate MSE
if (training == 0 || training == 1 || training == 2)
  addpath(folder_datasets)
  data_piglets
  Ts = data(positions,  9); % skin-surface temperature (oC).
  Th = data(positions, 14); % hair-coat temperature (oC)
  
  Th_diff = predictionsEnsemble(:, 3) - Th;
  Ts_diff = predictionsEnsemble(:, 1) - Ts;
  mean_errors = NA*zeros(11, 1);
  %  1: Average time per core
  %  2: MSE for Ts
  %  3: MSE for Th
  %  4: Total MSE
  %  5: Mean standard deviation for Ts
  %  6: Mean standard deviation for Th
  %  7: Mean standard deviation for Ts and Th
  %  8: Mean standard deviation for q_skin
  %  9: Total negative log-likelihood
  % 10: Negative log-likelihood for Ts
  % 11: Negative log-likelihood for Th
  
  mean_errors(1) = totalTime;
  mean_errors(2) = Ts_diff'*Ts_diff/length(positions);
  
  if ( strcmp(solveFor, "KM") )
    mean_errors(3) = Th_diff'*Th_diff/length(positions);
    mean_errors(4) = ( mean_errors(2) + mean_errors(3) )/2;
  else
    mean_errors(4) = mean_errors(2);
  end
  
  mean_errors(5) = mean(predictionsEnsemble(:, 11));
  mean_errors(6) = mean(predictionsEnsemble(:, 13));
  mean_errors(7) = ( mean_errors(5) + mean_errors(6) )/2;
  
  mean_errors(8) = mean(predictionsEnsemble(:, 13));
  
  % calculating log-likelihoods
  folder_ensembeOptimumWeight = locateFolderPath("ensemble/optimumWeight");
  addpath(folder_ensembeOptimumWeight)
  
  half_n_All = size(predictionsEnsembleParallel, 1);
  predictionsAllSets = zeros(2*half_n_All, length(weightsNonZero));
  for ii = 1:length(setsNonZero)
    % calculating predictions and variables that are weighted
    predictionsAllSets(1:half_n_All, ii) = predictionsEnsembleParallel(:, 1 + (ii - 1)*10);
    predictionsAllSets( (half_n_All+1):end, ii) = predictionsEnsembleParallel(:, 3 + (ii - 1)*10);
  end

  % total negative log-likelihoods
  mean_errors(9) = objectiveFunctionOptimumWeightMSE(weightsNonZero', [Ts; Th], ...
                        predictionsAllSets, ...
                        0, length(Ts)*2);
  % negative log-likelihood for Ts
  mean_errors(10) = objectiveFunctionOptimumWeightMSE(weightsNonZero', Ts, ...
                        predictionsAllSets(1:half_n_All, :), ...
                        0, length(Ts));
  % negative log-likelihood for Th
  mean_errors(11) = objectiveFunctionOptimumWeightMSE(weightsNonZero', Th, ...
                        predictionsAllSets( (half_n_All+1):end, :), ...
                        0, length(Th));
                        
  
  disp(["Ensemble " datasetTypeName " Total machine time: " num2str(mean_errors(1), "%2.4f") ...
              "s; MSE (std; NLL): Total = " num2str(mean_errors(4), "%2.4f") ...
               " oC^2 (" num2str(mean_errors(7), "%2.4f") ...
               " oC; " num2str(mean_errors(9), "%2.4f") ...
               "); Ts = " num2str(mean_errors(2), "%2.4f") ...
               " oC^2 (" num2str(mean_errors(5), "%2.4f") ...
               " oC; " num2str(mean_errors(10), "%2.4f") ...
               "); Th = " num2str(mean_errors(3), "%2.4f") ...
               " oC^2 (" num2str(mean_errors(6), "%2.4f") " oC; " ...
               num2str(mean_errors(11), "%2.4f") ")"]);
               
  rmpath(folder_ensembeOptimumWeight)
end

%% saving
% we save in case we are running for training, testing or a name was given to the variable datasetTypeName
if (training == 0 || training == 1 || training == 2 || !isempty(nameEnsemble))
  foldername = [locateFolderPath(["ensembleData"]) "/" solveFor];
  
  filename = [foldername "/" nameEnsemble datasetTypeName ".mat"];

  if (exist("mean_errors", "var") == 1)
    save(filename, "-V7", "mean_errors", "totalTime", "predictionsEnsemble");
  else
    save(filename, "-V7", "totalTime", "predictionsEnsemble");
  end
  
end
