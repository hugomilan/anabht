function predictions = loadPredictions(sets, loadModel = "all", training = 1)
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
% File:   loadPredictions.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Load the predictions of the specified sets
%
%
%
% Usage:
% predictions = loadPredictions(sets, loadModel, training)
%
% Input:
% sets: a vector with the #sets to load.
% training: flag that indicates if should load predictions using the training set (1)
%            or the testing set (0).
% loadModel: name of the model to load. Valid inputs are:
%            (1) "all": loads solutions from all models.
%            (2) "Skin": loads the solutions from the problem that
%                        considers skin only
%            (3) "KM": loads the solutions from the problem that solves 
%                        the Kowalski-Mitchell formulation for radiation inside
%                        the hair-coat layer
%
% Output:
% predictions: matrix that contains the predictions. The first dimension is the
%              the sample and has size 200. The second dimension is the prediction
%              variable and has size 3 for (1) Th, (2) Ts, (3) qs. The third dimension
%              is the set_ID of the ith input of the variable "sets".


if (strcmp(loadModel, "all"))
  predictions.KM = loadMeanErrors(sets, "KM", training);
  predictions.Skin = loadMeanErrors(sets, "Skin", training);
  return
end 

% running for training and testing positions
if (training == 1)
  % we will run for the training positions
  datasetTypeName = "Training";
elseif (training == 0)
  % we will run for the testing positions
  datasetTypeName = "Testing";
else 
% will run for both, testing and training.
  datasetTypeName = "TestingAndTraining";
end

predictions = zeros(200, 3, length(sets));
for ii = 1:length(sets)
  set_i = sets(ii);
  filename = ["setsData/" loadModel "/S" num2str(set_i) "data" datasetTypeName ".mat"];
  if (exist(filename) != 2)
    disp(["Could not find the file " filename]);
    predictions(:, :, ii) = NA;
    continue
  end
  load(filename);

  predictions(:, :, ii) = data_out(:,[3 4 9]);
end

end
