function testHairCoatDivision(minN, maxN, training = 1)
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
% File:   testHairCoatDivision.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on May 5, 2018.
%
%
% Function description:
% Calculates Ts and Th for chaning the division of the hair coat between minN and maxN
%
%
% Usage:
% 1) Calculates Ts and Th for chaning the division of the hair coat between minN and maxN
% testHairCoatDivision(minN, maxN);
% 
%
% Input:
% minN: starting number of divisions of the hair coat
% maxN: final number of divisions of the hair coat
% training: Should calculate for training, testing, or both? If 0, calculates 
%           for testing. If 1, calculates for training. If 2, calculates for
%           both, testing and training.
% 
% Output:
% None
%

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
else 
% will run for both, testing and training.
  load([folder_datasets "/trainingPosition.mat"])
  load([folder_datasets "/testingPosition.mat"])
  positions = [trainingPosition' testingPosition];
  datasetTypeName = "TestingAndTraining";
end



%% Adding the necessary paths
set_i = 0;
% obtaining the sets
folder_sets = locateFolderPath("sets");
filename = [folder_sets "/S" num2str(set_i) ".mat"];
predictionsSet = cell(maxN - minN, 1);
dataDivisionOut = zeros(maxN - minN, 11);

for ii = minN:maxN
  disp(["Calculating for " num2str(ii) " divisions of the hair coat."])
  [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
    d, epsilong, dg, ...
    z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
    To, qo, ua, Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getInputs(ii);
    
  % saving the set and the calculated properties
  save(filename, "-V7", "set_i", "k", "kh", "wb", "rhob", "cb", "Tb_m", ...
  "qtotal", "L", "N", "D", "HL", "epsilon", "d", "epsilong", "dg", "z", ...
  "Lt", "TMR_m", "h_m", "h_m_skin", "omega", "phi", ...
  "NMuscleLayers", "NFatLayers", "NSkinLayers", "NHairCoatLayers", ...
  "To", "qo", "ua", "Ta_pen_stderr", "Ta_brooder_stderr", ...
  "Tg_brooder_stderr", "Tr_stderr");
  
  % solving for the specified number of layers of the hair coat
  timei = time;
  predictionsSet{ ii - minN + 1 } = solveSetKM(0, NA, NA, NA, NA, NA, positions, datasetTypeName, 0, 0);
  dataDivisionOut( ii - minN + 1, 11) = time - timei;
  dataDivisionOut( ii - minN + 1, 1) = ii;
  Ts = predictionsSet{ ii - minN + 1 }(find( !isnan( predictionsSet{ ii - minN + 1 }(:,1)  ) ) , 1);
  qs = predictionsSet{ ii - minN + 1 }(find( !isnan( predictionsSet{ ii - minN + 1 }(:,2)  ) ) , 2);
  Th = predictionsSet{ ii - minN + 1 }(find( !isnan( predictionsSet{ ii - minN + 1 }(:,3)  ) ) , 3);
  dataDivisionOut( ii - minN + 1,  2) = min(Ts);
  dataDivisionOut( ii - minN + 1,  3) = mean(Ts);
  dataDivisionOut( ii - minN + 1,  4) = max(Ts);
  dataDivisionOut( ii - minN + 1,  5) = min(qs);
  dataDivisionOut( ii - minN + 1,  6) = mean(qs);
  dataDivisionOut( ii - minN + 1,  7) = max(qs);
  dataDivisionOut( ii - minN + 1,  8) = min(Th);
  dataDivisionOut( ii - minN + 1,  9) = mean(Th);
  dataDivisionOut( ii - minN + 1, 10) = max(Th);
  disp( num2str(dataDivisionOut( ii - minN + 1, :)) );
  disp(["Done."])
  
end


folder_dataOutput = locateFolderPath("dataOutput");

filename = [folder_dataOutput "/HairCoatDivision_" num2str(minN) "_" num2str(maxN) "data" datasetTypeName "KM.mat"];
save(filename, "dataDivisionOut");

end
