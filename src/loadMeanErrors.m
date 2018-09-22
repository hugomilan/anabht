function mean_errors_out = loadMeanErrors(sets, loadModel = "all", training = 1)
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
% File:   loadMeanErrors.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 13, 2018.
%
%
% Function description:
% Contains the input parameters and its calculations for implementing the 
% analytical model. Outputs the variables required for the computations. 
% Modify this function or create another one if you want to solve for another problem.
%
%
%
% Usage:
% mean_errors = loadMeanErrors(sets, loadModel, training)
%
% Input:
% sets: a vector with the #sets to load.
% loadModel: name of the model to load. Valid inputs are:
%            (1) "all": loads solutions from all models.
%            (2) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (3) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
% training: Should calculate MSE for training, testing, or both? If 0, calculates 
%           MSE for testing. If 1, calculates MSE for training. If 2, calculates MSE for
%           both, testing and training.
%
% Output:
% mean_errors: matrix of mean errors. Lines are different sets and columns are the
% errors/data of the set

if (strcmp(loadModel, "all"))
  mean_errors_out.KM = loadMeanErrors(sets, "KM", training);
  mean_errors_out.Skin = loadMeanErrors(sets, "Skin", training);
  return
end  

% loading training and testing positions
if (training == 1)
  % we will run the model for the training positions
%  load("datasets/trainingPosition.mat")
%  positions = trainingPosition;
  datasetTypeName = "Training";
elseif (training == 0)
  % we will run the model for the testing positions
%  load("datasets/testingPosition.mat")
%  positions = testingPosition;
  datasetTypeName = "Testing";

else 
% will run for both, testing and training.
%  load("datasets/trainingPosition.mat");
%  load("datasets/testingPosition.mat")
%  positions = [trainingPosition' testingPosition];
  datasetTypeName = "TestingAndTraining";
end

mean_errors_out = zeros(length(sets), 4);
folder_setsData = locateFolderPath("setsData");
for ii = 1:length(sets)
  set_i = sets(ii);
  filename = [folder_setsData "/" loadModel "/S" num2str(set_i) "data" datasetTypeName ".mat"];
  if (exist(filename) != 2)
    disp(["Could not find the file " filename]);
    mean_errors_out(ii, :) = NA;
    continue
  end
  load(filename);

  mean_errors_out(ii, :) = mean_errors;
end

end
