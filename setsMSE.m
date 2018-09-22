function mean_errors = setsMSE(set_i, loadSets = 0, training = 1)
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
% File:   setsMSE.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 2, 2018.
%
%
% Function description:
% Given the #set, creates a random set, solves it for training/testing dataset,
% and outputs the mean squared error for that #set for Th and Ts and the 
% time it took to run it
%
%
% Usage:
% 1) Minimum usage. Uses default values
% mean_errors = setsMSE(set_i);
%
% 2) Loads previously #set if they were already sampled
% mean_errors = setsMSE(set_i, 1);
%
% 2) Loads previously #set if they were already sampled and calculates the MSE
% for the testing dataset
% mean_errors = setsMSE(set_i, 1, 0);
% 
%
% Input:
% set_i: #set
% loadSets: Should load previously sampled sets or sample new sets? If 0, sample
% new sets. If 1, uses previously sampled sets whenever available.
% training: Should calculate MSE for training, testing, or both? If 0, calculates 
% MSE for testing. If 1, calculates MSE for training. If 2, calculates MSE for
% both, testing and training.
% 
% Output:
% mean_errors: mean squared errors for (1) both, Th and Ts, (2) Ts, (3) and Th,
% (4) and the time it took to run this function.
%

% Including the folder paths in the search directories
if exist("datasets") != 0
  addpath("datasets");
end
if exist("src") != 0
  addpath("src");
end
% adding the location of the java class
javaaddpath("java/");
% adding the dependency
javaaddpath("java/h2o-genmodel.jar");
% instantiating the prediction model
predictTaTgTr = javaObject("PredictTaTgTr");

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  more off
else
  % we are executing in MATLAB
end

% Load input data
data_piglets % this file generates data matrix
% Columns are:
%  1: Day
%  2: Time_in_degress
%  3: Hour
%  4: Brooder
%  5: Heat_source_power (W)
%  6: Sequence
%  7: Piglet
%  8: Tr (oC)
%  9: Ts (oC)
% 10: Ta_brooder (oC)
% 11: RH_brooder (%)
% 12: Tg_brooder (oC)
% 13: Ta_pen (oC)
% 14: Th_IR (oC)
% 15: Irradiance (W/m2)

% loading training and testing positions
if (training == 1)
  % we will run the model for the training positions
  load("datasets/trainingPosition.mat")
  positions = trainingPosition;
  datasetTypeName = "Training";
elseif (training == 0)
  % we will run the model for the testing positions
  load("datasets/testingPosition.mat")
  positions = testingPosition;
  datasetTypeName = "Testing";

else 
% will run for both, testing and training.
  load("datasets/trainingPosition.mat");
  load("datasets/testingPosition.mat")
  positions = [trainingPosition' testingPosition];
  datasetTypeName = "TestingAndTraining";
end


% This will contain the number of the sample, the predicted skin surface and
% hair coat surface temperatures, and their errors
data_out = zeros(200, 12);
%  1: #sample
%  2: #set
%  3: Th_ana
%  4: Ts_ana
%  5: Th_ana_diff
%  6: Ts_ana_diff
%  7: Th_ana_percentage
%  8: Ts_ana_percentage
%  9: q_skin_ana
% 10: Time taken to run for the #sample
% 11: Value of the objective function
% 12: Maximum value of the constraint function
% 13: Ta_pen used in the calculation
% 14: Ta_brooder used in the calculation
% 15: Tg_brooder used in the calculation
% 16: Tr used in the calculation
data_out(:, :) = NA;
data_out(:, 1) = 1:200;
data_out(:, 2) = set_i;


mean_errors = zeros(1, 5);
%  1: set number
%  2: Mean squared error for skin surface temperature
%  3: Mean squared error for hair-coat surface temperature
%  4: Mean squared error for skin surface and hair-coat surface temperature
%  5: Total time it took to run
mean_errors(:,:) = NA;

time_i = time;

% obtaining the sets
filename = ["sets/S" num2str(set_i) ".mat"];
if (loadSets == 1 && exist(filename) == 2)
  load(filename);
  
else
  if (set_i == 0 || set_i == 1)
    % if set_i == 0 or set_i == 1, we are getting deterministic values
    % They don't have offsets but set_i = 0 uses all the measured data while
    % set_i == 1 uses the best estimations
    [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
    rho, alpha, rhoh, alphah, d, epsilong, dg, ...
    z, Lt, TMR_m, h_m, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
    To, qo, ua, Ta_pen_offset, Ta_brooder_offset, Tg_brooder_offset, Tr_offset] = getInputs();
    
    % Offsets for the predicted temperatures
    % offset for the air temperature measured in the pen
    Ta_pen_offset = set_i*normrnd(0, 0.35);
    % offset for the predicted air temperature inside the brooder
    Ta_brooder_offset = set_i*normrnd(0, 1.64);
    % offset for the predicted black globe temperature inside the brooder
    Tg_brooder_offset = set_i*normrnd(0, 0.69);
    % offset for the predicted rectal temperature
    Tr_offset = set_i*normrnd(0, 0.55);
  else 
    % otherwise, we are getting random sets
    [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
    rho, alpha, rhoh, alphah, d, epsilong, dg, ...
    z, Lt, TMR_m, h_m, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
    To, qo, ua, Ta_pen_offset, Ta_brooder_offset, Tg_brooder_offset, Tr_offset] = getRandomInputs();
  end

  % calculating parameters that will be used during calculations
  NTissueLayers = NMuscleLayers + NFatLayers + NSkinLayers;
  Lh = sum(L((end - NHairCoatLayers + 1):end));
  alpha_eff_KM = 2/(3*pi)*(N.*HL/Lh).*D; %effective absorption coefficient of the hair coat for the KM model
  % N_effective accordingly to Kowalski (1978)
  
  [A, b, A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin, ...
  matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
  getProblemIndependentMatrices(k, L, qtotal, wb, rhob, cb, ...
  To, qo, NTissueLayers, alpha_eff_KM);
  
  % obtaining a first estimation of beta
  T_initial = 37;
  temp_offset_initial = updateTempParametersProblemDependent(...
       temp_offset, wb(1:NTissueLayers), T_initial*ones(1, length(temp_offset)) );
  beta_initial = matrix_temp\(T_initial - temp_offset_initial);
  
  % saving the set and the calculated properties
  save(filename, "-V7", "set_i", "k", "kh", "wb", "rhob", "cb", "Tb_m", ...
  "qtotal", "L", "N", "D", "HL", "epsilon", "rho", "alpha", "rhoh", ...
  "alphah", "d", "epsilong", "dg", "z", "Lt", "TMR_m", "h_m", ...
  "NMuscleLayers", "NFatLayers", "NSkinLayers", "NHairCoatLayers", ...
  "NTissueLayers", "Lh", "To", "qo", "ua", "alpha_eff_KM", "A", "b", ...
  "A_FD", "E", "b_FD1", "A_q_skin", "E_q_skin", "b_q_skin", ...
  "matrix_temp", "temp_offset", "matrix_heat", "heat_offset", "T_initial", ...
  "beta_initial", "temp_offset_initial", "Ta_pen_offset", "Ta_brooder_offset", ...
  "Tg_brooder_offset", "Tr_offset");
  
end



Tb = zeros(1, NTissueLayers + NHairCoatLayers);
sigma = 5.670373e-8; % Stefan-Boltzmann constant
% tolerances in case the objective does not converge
objTol = 10;
eqTol = 10;
for position_i = 1:length(positions)
  time_set_i = time;
  sample_i = positions(position_i);
  
  % loading the environmental variables
  Ta_pen = data(sample_i, 13) + Ta_pen_offset; % Air temperature in the pen (oC)
  heat_power = data(sample_i, 5); % heat source power (W)
  
  % obtaining the predicted temperatures inside the brooder and the rectal temperature
  if (set_i == 0)
    Ta_brooder = data(sample_i, 10); % Air temperature inside the brooder(oC)
    Tg_brooder = data(sample_i, 12); % Black globe temperature inside the brooder(oC)
    Tr = data(sample_i,  8); % Rectal Temperature (oC)
  else
    predictTaTgTr.predict(Ta_pen, heat_power);
    Ta_brooder = predictTaTgTr.Ta_brooder + Ta_brooder_offset;
    Tg_brooder = predictTaTgTr.Tg_brooder + Tg_brooder_offset;
    Tr = predictTaTgTr.Tr + Tr_offset;
  end
  % values for validation. Only used to calculate the model error
  Ts = data(sample_i,  9); % skin-surface temperature (oC).
  Th = data(sample_i, 14); % hair-coat temperature (oC)
  
  % saving input variables
  data_out(sample_i, 13) = Ta_pen;
  data_out(sample_i, 14) = Ta_brooder;
  data_out(sample_i, 15) = Tg_brooder;
  data_out(sample_i, 16) = Tr;
  
  % Blood temperature
  for ii = 1:NTissueLayers
    Tb(ii) = Tr*Tb_m(ii);
  end
  % mean radiant temperature
  TMR = getTMR(Tg_brooder, Ta_brooder, ua, dg, epsilong, z, Lt);
  TMR *= TMR_m;
  
  % Effective hair coat thermal conductivity
  for ii = 1:NHairCoatLayers
    k(NTissueLayers + ii) = getEffectiveHairCoatThermalConductivity(Ta_brooder, N(ii), D(ii), HL(ii), Lh, kh(ii), z, Lt);
  end
  
  % Updating the matrices
  [b_updated, A_FD_updated, b_FD_updated, A_q_skin_updated, ...
  b_q_skin_updated, temp_offset_updated] = ...
  updateProblemDependentMatrices(b, A_FD, b_FD1, A_q_skin, ...
  b_q_skin, temp_offset, k, wb, Tb, Tr, TMR, NTissueLayers);
  
  
  % Kowalski and Mitchel (1979) for long-wave radiation:        (1)
  [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
      solveAnalyticalModelKM(Ta_brooder, k, L, qo, d, z, Lt, h_m, ua, NTissueLayers, ...
   A, b_updated, A_FD_updated, E, b_FD_updated, A_q_skin_updated, E_q_skin, b_q_skin_updated, ...
   matrix_temp, temp_offset_updated, matrix_heat, heat_offset, ...
   beta_initial, T_initial);
  data_out(sample_i, 11) = obj;
  data_out(sample_i, 12) = max(abs(eqConst));
  
  %% if info != 101, we got problems.
  if (info == 101)
    % success. Do nothing
  elseif (info == 104)
    % the stepsize has become too small. Check if the obj and eqConst are OK
    if (obj > objTol || max(abs(eqConst)) > eqTol)
      % something was not correct. Continue and it will save NA
      disp(["#set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
      continue
    end
  elseif (info == 102)
    % The BFGS update falied
      % something was not correct. Continue and it will save NA
      disp(["#set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
    continue
  elseif (info == 103)
    % the maximum number of iterations was reached. Check if the obj and eqConst are OK
    if (obj > objTol || max(abs(eqConst)) > eqTol)
      % something was not correct. Continue and it will save NA
      disp(["#set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
      continue
    end
  end
  
  % calculate errors
  Th_ana = T_int(end);
  Ts_ana = T_int(NTissueLayers + 1);
  data_out(sample_i, 3) = Th_ana;
  data_out(sample_i, 4) = Ts_ana;

  Th_diff = Th_ana - Th;
  Ts_diff = Ts_ana - Ts;
  data_out(sample_i, 5) = Th_diff;
  data_out(sample_i, 6) = Ts_diff;

  Th_error = Th_diff/Th*100;
  Ts_error = Ts_diff/Ts*100;
  data_out(sample_i, 7) = Th_error;
  data_out(sample_i, 8) = Ts_error;
  
  % Calculating heat flux at the skin
  data_out(sample_i, 9) = q_int(NTissueLayers + 1); % q_skin_ana;
  % Time it took to run for this sample
  data_out(sample_i, 10) = time - time_set_i;
  
%  disp(["#Set: " num2str(set_i) "; #Sample: " num2str(sample_i) ...
%  "; q_skin = " num2str(q_int(NTissueLayers + 1), "%2.4f") " W/m2" ...
%  "; Errors: Ts = " num2str(Ts_diff, "%2.4f") ...
%  " oC; Th = " num2str(Th_diff, "%2.4f") " oC; Time: " ...
%  num2str(data_out(sample_i, 10), "%2.2f") " s"]);
%
%  [T_plot q_plot x_plot] = getPlotValues(beta, 100, k, L, ...
%   qtotal, wb, rhob, cb, Tb, T_int, q_int, NTissueLayers);
% 
%  plot(x_plot, T_plot, '-', "linewidth", 4)
%  hold on
%  Lm = sum(L(1:(NMuscleLayers)));
%  Lf = sum(L((NMuscleLayers + 1):(NMuscleLayers + NFatLayers)));
%  Ls = sum(L((NMuscleLayers + NFatLayers + 1):(NMuscleLayers + NFatLayers + NSkinLayers)));
%  plot([Lm, Lm], [20, 50], '--', "linewidth", 3)
%  plot([Lf + Lm, Lf + Lm], [20, 50], '--', "linewidth", 3)
%  plot([Ls + Lf + Lm, Ls + Lf + Lm], [20, 50], '--', "linewidth", 3)
%  ylim([20,45])
%  xlim([0, Lm + Lf + Ls + Lh])
%  plot([Lm + Lf + Ls, Lm + Lf + Ls + Lh], [Ts, Th], '*', "linewidth", 5, "markersize", 10)
%  pause(0.01)
%  hold off

end
totalTime = time - time_i; % time taken to run all models

mean_errors(1) = set_i;
mean_errors(2) = data_out(positions, 6)'*data_out(positions, 6)/length(positions);
mean_errors(3) = data_out(positions, 5)'*data_out(positions, 5)/length(positions);
mean_errors(4) = (mean_errors(2) + mean_errors(3))/2;
mean_errors(5) = totalTime;

filename = ["setsData/S" num2str(set_i) "data" datasetTypeName ".mat"];
save(filename, "-V7", "data_out", "mean_errors");

disp([" #Set: " num2str(set_i) "; Time: " num2str(mean_errors(5), "%2.4f") ...
      " s; MSE: Total = " num2str(mean_errors(4), "%2.4f") ...
      " oC^2; Ts = " num2str(mean_errors(2), "%2.4f") ...
      " oC^2; Th = " num2str(mean_errors(3), "%2.4f") " oC^2"]);


end
