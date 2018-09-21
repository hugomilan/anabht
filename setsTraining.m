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
% File:   setsTraining.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Calls the functions that calculate the meas squared error on the training dataset.

% select the sets to train
% 0: deterministic set. Other numbers: random set
sets = 0:10001;

% If loadSets = 1, and if the #set was already sampled, the algorithm will use the
% previously sampled set
loadSets = 0; 

% solve for a specific model or for all
solveFor = "KM"; %"all"



%%%%%%%%%%%%%% NO NEED TO MODIFY BELOW THIS LINE %%%%%%%%%%%%%%

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  pkg load parallel
  more off
else
  % we are executing in MATLAB
end

% loading the path
if exist("src") != 0
  addpath("src");
end

time_i = time;
fun = @(x) problemsSetsMSE(x, loadSets, 1, solveFor);
error = @(x) disp(x);

% running
pararrayfun(3, fun, sets, 'ErrorHandler', error);

run_time = time - time_i
mean_run_time = run_time/length(sets);

% loading the MSE
MSEs = loadMeanErrors(sets, solveFor);
if (strcmp(solveFor, "all"))
  minMSEError = inf;
  minMSEErrorPosition = 0;
  bestModelName = NA;
  for [MSEModel_i, NameModel_i] = MSEs
    [minErrors minErrorsPosition] = min( MSEModel_i(:,4) );
    if (minErrors < minMSEError)
      minMSEError = minErrors;
      minMSEErrorPosition = minErrorsPosition;
      bestModelName = NameModel_i;
    end
  end
  
  
else
  [minErrors minErrorsPosition] = min( MSEs(:,4) );
  bestModelName = solveFor;
end

disp(["Time: mean = " num2str(mean_run_time, "%2.4f") " s; " ...
      "total = " num2str(run_time, "%2.4f") " s; " ...
      "Minimum MSE: " num2str(minErrors, "%2.4f") " oC^2 " ...
      "for #set " num2str(sets(minErrorsPosition)) " and " ...
      bestModelName " model"])

% Saving data
filename = ["setsData/TimeS" num2str(sets(1)) "-" num2str(sets(end)) ".mat"];
save(filename, "-V7", "run_time", "mean_run_time");
