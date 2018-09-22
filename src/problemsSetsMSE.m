function problemsSetsMSE(set_i, loadSets = 0, training = 1, solveFor = "all", ...
                         heat_power_in = NA, Ta_pen_in = NA, Ta_brooder_in = NA, ...
                         Tg_brooder_in = NA, Tr_in = NA)
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
% File:   problemsSetsMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 26, 2018.
%
%
% Function description:
% Given the #set, creates a random set, and calls the functions that solve the 
% problems for training/testing dataset,
% and outputs the mean squared error for that #set for Th and Ts and the 
% time it took to run it
%
%
% Usage:
% 1) Minimum usage. Uses default values
% problemsSetsMSE(set_i);
%
% 2) Loads previously #set if they were already sampled
% problemsSetsMSE(set_i, 1);
%
% 3) Loads previously #set if they were already sampled and calculates the MSE
% for the testing dataset
% problemsSetsMSE(set_i, 1, 0);
%
% 4) Loads previously #set if they were already sampled and calculates the MSE
% for the testing dataset using the "Skin" model
% problemsSetsMSE(set_i, 1, 0, "Skin");
% 
%
% Input:
% set_i: #set
% loadSets: Should load previously sampled sets or sample new sets? If 0, sample
%           new sets. If 1, uses previously sampled sets whenever available.
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
% heat_power: Power of the supplemental heat source (W). Can be a vector. Used to
%             estimate Ta_pen, Ta_brooder, Tg_brooder, and Tr. If NA, will check
%             if positions is NA. If so, will throw an error. If not, will load the
%             measured dataset and solve the problem for the given positions.
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
%
% Output:
% None
%

% timing
time_i = time;

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  more off
else
  % we are executing in MATLAB
end

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



% obtaining the sets
folder_sets = locateFolderPath("sets");
filename = [folder_sets "/S" num2str(set_i) ".mat"];
if (loadSets == 1 && exist(filename) == 2)
  load(filename);
  
else
  if (set_i == 0 || set_i == 1)
    % if set_i == 0 or set_i == 1, we are getting deterministic values
    % They don't have offsets but set_i = 0 uses all the measured data while
    % set_i == 1 uses the best estimations
    [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
    d, epsilong, dg, ...
    z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
    To, qo, ua, Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getInputs();
    
    
  else 
    % otherwise, we are getting random sets
    [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
    d, epsilong, dg, ...
    z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
    To, qo, ua, Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getRandomInputs();
  end
  
  % saving the set and the calculated properties
  save(filename, "-V7", "set_i", "k", "kh", "wb", "rhob", "cb", "Tb_m", ...
  "qtotal", "L", "N", "D", "HL", "epsilon", "d", "epsilong", "dg", "z", ...
  "Lt", "TMR_m", "h_m", "h_m_skin", "omega", "phi", ...
  "NMuscleLayers", "NFatLayers", "NSkinLayers", "NHairCoatLayers", ...
  "To", "qo", "ua", "Ta_pen_stderr", "Ta_brooder_stderr", ...
  "Tg_brooder_stderr", "Tr_stderr");
  
end


%% Adding the necessary paths
folder_analyticalModels_KM = locateFolderPath("analyticalModels/KM");
addpath(folder_analyticalModels_KM);
%% Adding the necessary paths
folder_analyticalModels_skin = locateFolderPath("analyticalModels/Skin");
addpath(folder_analyticalModels_skin);

% function that solves KM
if ( strcmp(solveFor, "all") || strcmp(solveFor, "KM") )
  solveSetKM(set_i, heat_power_in, Ta_pen_in, Ta_brooder_in, ...
             Tg_brooder_in, Tr_in, positions, datasetTypeName, loadSets, 0);
end
% function that solves Skin
if (strcmp(solveFor, "all") || strcmp(solveFor, "Skin") )
  solveSetSkin(set_i, heat_power_in, Ta_pen_in, Ta_brooder_in, ...
               Tg_brooder_in, Tr_in, positions, datasetTypeName, loadSets, 0);
end

% cleaning after myself
rmpath(folder_analyticalModels_KM)
rmpath(folder_analyticalModels_skin)


folder_setsData = locateFolderPath("setsData");

totalTime = time - time_i; % time taken to run a set for all the problems and datapoints
filename = [folder_setsData "/S" num2str(set_i) "data" datasetTypeName solveFor "TotalTime.mat"];
save(filename, "-V7", "totalTime");

end
