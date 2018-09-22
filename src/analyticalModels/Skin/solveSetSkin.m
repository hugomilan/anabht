function predictionsSet = solveSetSkin(set_i, heat_power_in = NA, Ta_pen_in = NA, ...
                        Ta_brooder_in = NA, Tg_brooder_in = NA, Tr_in = NA, ...
                        positions = NA, datasetTypeName = "", loadSets = 1, ...
                        dispPartialResults = 1, plotPartialResults = 0, ...
                        loadEstimations = 0, Ta_brooder_in_increase = 0)
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
% File:   solveSetSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 27, 2018.
%
%
% Function description:
% Given the #set, loads the parameters of that set and calculates Ts and q_skin.
% If running for the measured dataset, will calculate the mean squared error.
%
%
% Usage:
% predictionsSet = solveSetSkin(set_i, heat_power_in, Ta_pen_in, Ta_brooder_in, ...
%                      Tg_brooder_in, Tr_in, positions, datasetTypeName, loadSets, ...
%                      dispPartialResults, plotPartialResults, loadEstimations, ...
%                      Ta_brooder_in_increase)
% 
%
% Input:
% set_i: #set
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
% positions: Vector with positions of the measured dataset to use. Used only if 
%            heat_power = NA
% datasetTypeName: String that append to the name of the file saved if computed for
%            the measured dataset.
% loadSets: Should load previously data previously calculated or load the base data
%           and do the calculations? If 0, load the base data and do the calculations.
%           If 1, uses the previously calculated dataset.
% dispPartialResults: should display what are the partial computations? 1: yes, 2: no.
% plotPartialResults: should plot the profile of the partial computations? 0: no,
%                     anything different than zero means yes and indicates how long
%                     should the pause be for every plot.
% loadEstimations: when running for training or testing:
%                       0: recalculates the predictions even if they are available
%                       1: load predictions if they are available
% Ta_brooder_in_increase: when Ta_brooder is an input, should it be add by the 
%                         the effect of heat_source_power?
%                           0: Ta_brooder is as the input.
%                           1: Ta_brooder is the input + the added effect of
%                              heat_soure_power in increasing Ta_brooder
% 
% Output:
% predictionsSet: matrix with the predictions from the model. Rows are the different
%                 inputs and columns are the different variables predicted. Columns are:
%                   1) skin-surface temperature (Ts)
%                   2) skin-surface heat flux (q_skin)
%                   3) NA (this would be hair-coat surface temperature, Th, but this model does
%                      not have it; it assumes the hair-coat is negligible)
%                   4) time it took to solve for this position
%                   5) value of the objective function
%                   6) maximum value of the constraint function
%                   7) Ta_pen used in the calculation
%                   8) Ta_brooder used in the calculation
%                   9) Tg_brooder used in the calculation
%                  10) Tr used in the calculation.
%

%% Adding the necessary paths
folder_analyticalModels = locateFolderPath("analyticalModels");
addpath(folder_analyticalModels);

% testing inputs
[loadMeasuredData, data, data_out, mean_errors, predictionsSet, positions] = ...
                      testInputs(set_i, heat_power_in, Ta_pen_in, ...
                        Ta_brooder_in, Tg_brooder_in, Tr_in, positions);


% testing if we are loading previous estimations
if (loadEstimations == 1)
  folder_setsData = locateFolderPath("setsData");
  foldername = [folder_setsData "/KM"];
  filename = [foldername "/S" num2str(set_i) "data" datasetTypeName ".mat"];
  % if the file exists, we load it and return
  if (exist(filename) == 2)
    load(filename)
    % Converting data_out to predictionsSet
    predictionsSet = data_out(positions, [4, 9, 3, 10:16]);
    
%  1: skin-surface temperature (Ts)            %  4: Ts_ana
%  2: skin-surface heat flux (q_skin)          %  9: q_skin_ana
%  3: hair-coat surface temperature (Th)       %  3: Th_ana
%  4: time it took to solve for this position  % 10: Time taken to run for the #sample
%  5: value of the objective function          % 11: Value of the objective function
%  6: maximum value of the constraint function % 12: Maximum value of the constraint function
%  7: Ta_pen used in the calculation           % 13: Ta_pen used in the calculation
%  8: Ta_brooder used in the calculation       % 14: Ta_brooder used in the calculation
%  9: Tg_brooder used in the calculation       % 15: Tg_brooder used in the calculation
% 10: Tr used in the calculation.              % 16: Tr used in the calculation
    return
  end
end


%%%%% Not using machine learning anymore
% adding the location of the java class. Used to predict Ta_brooder, Tg_brooder and Tr
%folder_java = locateFolderPath("java");
%javaaddpath([folder_java "/"]);
% adding the dependency
%javaaddpath([folder_java "/h2o-genmodel.jar"]);
% instantiating the prediction model
%predictTaTgTr = javaObject("PredictTaTgTr");
folder_LM_matrix = locateFolderPath("R/ML_training/Best ML models/LM_matrix");
predictTaTgTr.LM_Ta_brooder = load([folder_LM_matrix "/Ta_brooder.mat"]);
predictTaTgTr.LM_Ta_brooder_inc = load([folder_LM_matrix "/Ta_brooder_inc.mat"]);
predictTaTgTr.LM_Tg_brooder = load([folder_LM_matrix "/Tg_brooder.mat"]);
predictTaTgTr.LM_Tg_brooder_Ta_pen = load([folder_LM_matrix "/Tg_brooder_Ta_pen.mat"]);
predictTaTgTr.LM_Tr = load([folder_LM_matrix "/Tr.mat"]);

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  more off
else
  % we are executing in MATLAB
end

time_i = time;

% obtaining the sets
folder_sets = locateFolderPath("sets");
filename = [folder_sets "/Skin/S" num2str(set_i) ".mat"];
if (loadSets && exist(filename) == 2)
    % if the specific load file of the Skin model does not exist, we try to load
    % the general file. If it does not exist, we throw an error.
  load(filename);
  
else
  filename = [folder_sets "/S" num2str(set_i) ".mat"];
  if ( exist(filename) != 2 )
    error(["Skin: The file with the data for set " num2str(set_i) " was not found"]);
  end

  load(filename);
  % calculating parameters that will be used during calculations
  NTissueLayers = NMuscleLayers + NFatLayers + NSkinLayers;
  % resizing the parameters
  k = k(1:NTissueLayers);
  L = L(1:NTissueLayers);
  qtotal = qtotal(1:NTissueLayers);
  wb = wb(1:NTissueLayers);
  rhob = rhob(1:NTissueLayers);
  cb = cb(1:NTissueLayers);
  To = To(1:(NTissueLayers + 1)  );
  qo = qo(1:(NTissueLayers + 1)  );
  
  [A, b, matrix_temp, temp_offset, matrix_heat, heat_offset] = ...
  getProblemIndependentMatricesSkin(k, L, qtotal, wb, rhob, ...
  cb, To, qo);
  
  % obtaining a first estimation of beta
  T_initial = 37;
  temp_offset_initial = updateTempParametersProblemDependent(...
       temp_offset, wb, T_initial*ones(1, length(temp_offset)) );
  beta_initial = matrix_temp\(T_initial - temp_offset_initial);
  
  % saving the set and the calculated properties
  foldername = [folder_sets "/Skin"];
  if ( exist(foldername) != 7)
    mkdir(foldername);
  end
  filename = [foldername "/S" num2str(set_i) ".mat"];
  save(filename, "-V7", "set_i", "k", "wb", "rhob", "cb", "Tb_m", ...
  "qtotal", "L", "N", "D", "HL", "epsilon", ...
  "d", "epsilong", "dg", "z", "Lt", "TMR_m", "h_m", ...
  "NMuscleLayers", "NFatLayers", "NSkinLayers", ...
  "NTissueLayers", "To", "qo", "ua", "A", "b", ...
  "matrix_temp", "temp_offset", "matrix_heat", "heat_offset", "T_initial", ...
  "beta_initial", "temp_offset_initial", "Ta_pen_stderr", "Ta_brooder_stderr", ...
  "Tg_brooder_stderr", "Tr_stderr");
  
end

% adding the environment folder into the path
folder_environment = locateFolderPath("environment");
addpath(folder_environment);

Tb = zeros(1, NTissueLayers);
sigma = 5.670373e-8; % Stefan-Boltzmann constant
% tolerances in case the objective does not converge
objTolPerc = 2;
eqTol = 1e-3;
for position_i = 1:length(positions)
  time_position_i = time;
  sample_i = positions(position_i);
  
  % getting predictions or loading measured data
  [heat_power, Ta_pen, Ta_brooder, Tg_brooder, Tr] = ...
            getPredictions(set_i, sample_i, position_i, loadMeasuredData, ...
              data, heat_power_in, Ta_pen_in, Ta_brooder_in, Tg_brooder_in, ...
              Tr_in, predictTaTgTr, Ta_pen_stderr, Ta_brooder_stderr, ...
              Tg_brooder_stderr, Tr_stderr, Ta_brooder_in_increase);
              
              
  % if we are loading measured data
  if (loadMeasuredData)    
    % values for validation. Only used to calculate the model error
    Ts = data(sample_i,  9); % skin-surface temperature (oC).
    Th = data(sample_i, 14); % hair-coat temperature (oC)
    
    % saving input variables
    data_out(sample_i, 13) = Ta_pen;
    data_out(sample_i, 14) = Ta_brooder;
    data_out(sample_i, 15) = Tg_brooder;
    data_out(sample_i, 16) = Tr;
  end
  
  % saving data into the output variable
  predictionsSet(position_i,  7) = Ta_pen;
  predictionsSet(position_i,  8) = Ta_brooder;
  predictionsSet(position_i,  9) = Tg_brooder;
  predictionsSet(position_i, 10) = Tr;
  
  % Blood temperature
  for ii = 1:NTissueLayers
    Tb(ii) = Tr*Tb_m(ii);
  end
  % mean radiant temperature
  TMR = getTMR(Tg_brooder, Ta_brooder, ua, dg, epsilong, z, Lt);
  TMR *= TMR_m;
  
  % Updating the matrices
  [b_updated, temp_offset_updated] = ...
  updateProblemDependentMatricesSkin(b, temp_offset, wb, Tb, Tr);
  
  
  % Solve for when the hair-coat layer is neglected
  [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
      solveAnalyticalModelSkin(Ta_brooder, d, z, Lt, h_m, ua, A, b_updated, ...
   TMR, epsilon, qo, matrix_temp, temp_offset_updated, matrix_heat, heat_offset, beta_initial);
  
  objTol = (objTolPerc*q_int(end)/100)^2;
  predictionsSet(position_i, 5) = obj;
  predictionsSet(position_i, 6) = max(abs(eqConst));
  if (loadMeasuredData)
    data_out(sample_i, 11) = obj;
    data_out(sample_i, 12) = predictionsSet(position_i, 6);
  end
  
  %% if info != 101, we got problems.
  if (info == 101)
    % success. Do nothing
  elseif (info == 104)
    % the stepsize has become too small. Check if the obj and eqConst are OK
    if (obj > objTol || max(abs(eqConst)) > eqTol)
      % something was not correct. Continue and it will save NA
      disp(["Skin: #set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
      continue
    end
  elseif (info == 102)
    % The BFGS update falied
      % something was not correct. Continue and it will save NA
      disp(["Skin: #set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
    continue
  elseif (info == 103)
    % the maximum number of iterations was reached. Check if the obj and eqConst are OK
    if (obj > objTol || max(abs(eqConst)) > eqTol)
      % something was not correct. Continue and it will save NA
      disp(["Skin: #set: " num2str(set_i) " #sample: " num2str(sample_i) ...
      " Did not converge. Info = " num2str(info) ...
      "; Objective value = " num2str(obj, "%2.4f") ...
      "; Maximum constraint error = " num2str(max(abs(eqConst)), "%2.4f")]);
      continue
    end
  end
  
  Ts_ana = T_int(end);
  predictionsSet(position_i, 1) = Ts_ana;
  predictionsSet(position_i, 2) = q_int(end); % q_skin_ana;
  predictionsSet(position_i, 3) = NA;
  predictionsSet(position_i, 4) = time - time_position_i;
  
  if (loadMeasuredData)
    % calculate errors
    Th_ana = NA;
    data_out(sample_i, 3) = Th_ana;
    data_out(sample_i, 4) = Ts_ana;

    Th_diff = NA;
    Ts_diff = Ts_ana - Ts;
    data_out(sample_i, 5) = Th_diff;
    data_out(sample_i, 6) = Ts_diff;

    Th_error = NA;
    Ts_error = Ts_diff/Ts*100;
    data_out(sample_i, 7) = Th_error;
    data_out(sample_i, 8) = Ts_error;
  
    % Calculating heat flux at the skin
    data_out(sample_i, 9) = q_int(end); % q_skin_ana;
    % Time it took to run for this sample
    data_out(sample_i, 10) = predictionsSet(position_i, 4);
  end
  if (dispPartialResults)
    if(loadMeasuredData)
      msg = ["KM: #Set: " num2str(set_i) "; #Sample: " num2str(sample_i) ...
          "; q_skin = " num2str(q_int(NTissueLayers + 1), "%2.4f") " W/m2" ...
          "; Error: Ts = " num2str(Ts_diff, "%2.4f") ...
          " oC; Time: " num2str(data_out(sample_i, 10), "%2.2f") " s"];
    else
      msg = ["KM: #Set: " num2str(set_i) "; Position " num2str(position_i) "/" ...
          num2str(length(positions)) "; Ts = " num2str(Ts_ana, "%2.4f") ...
          " oC; q_skin = " num2str(q_int(NTissueLayers + 1), "%2.4f") " W/m2" ...
          "; Time: " num2str(predictionsSet(position_i, 4), "%2.2f") " s"];
    end
    disp(msg);
  end
  
  if (plotPartialResults)
    [T_plot q_plot x_plot] = getPlotValues(beta, 100, k, L, ...
     qtotal, wb, rhob, cb, Tb, T_int, q_int, NTissueLayers, 1);
 
    plot(x_plot, T_plot, '-', "linewidth", 4)
    hold on
    Lm = sum(L(1:(NMuscleLayers)));
    Lf = sum(L((NMuscleLayers + 1):(NMuscleLayers + NFatLayers)));
    Ls = sum(L((NMuscleLayers + NFatLayers + 1):(NMuscleLayers + NFatLayers + NSkinLayers)));
    plot([Lm, Lm], [20, 50], '--', "linewidth", 3)
    plot([Lf + Lm, Lf + Lm], [20, 50], '--', "linewidth", 3)
    plot([Ls + Lf + Lm, Ls + Lf + Lm], [20, 50], '--', "linewidth", 3)
    ylim([20,45])
    xlim([0, Lm + Lf + Ls + Lh])
    if (loadMeasuredData)
      plot([Lm + Lf + Ls], Ts, '*', "linewidth", 5, "markersize", 10)
    end
    pause(plotPartialResults)
    hold off
  end
end
totalTime = time - time_i; % time taken to run all models

if (loadMeasuredData)
  mean_errors(1) = set_i;
  mean_errors(2) = data_out(positions, 6)'*data_out(positions, 6)/length(positions);
  mean_errors(3) = NA;
  mean_errors(4) = mean_errors(2);
  
  folder_setsData = locateFolderPath("setsData");
  foldername = [folder_setsData "/Skin"];
  if ( exist(foldername) != 7)
    mkdir(foldername);
  end
  filename = [foldername "/S" num2str(set_i) "data" datasetTypeName ".mat"];
  save(filename, "-V7", "data_out", "mean_errors", "totalTime");

  disp([" Skin: #Set: " num2str(set_i) "; Time: " num2str(totalTime, "%2.4f") ...
        " s; MSE: Total = " num2str(mean_errors(4), "%2.4f") ...
        " oC^2; Ts = " num2str(mean_errors(2), "%2.4f") ...
        " oC^2; Th = " num2str(mean_errors(3), "%2.4f") " oC^2"]);
end

% cleaning after myself
rmpath(folder_analyticalModels);
rmpath(folder_environment);

end
