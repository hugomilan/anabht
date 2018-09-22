function [] = energyBalance(name, Ta_pen_surface_initial, heat_power_surface_initial, ...
                               weights = 1:0.1:20, CI_interval = 0.95, ...
                               N_random_samples = 1e4, printText = false)
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
% File:   energyBalance.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
% Created on June 23, 2018.
%
%
% Function description:
% Calculates the energy (in)balance between heat flowing to the environment and
% and generated heat. Calls the parallel version energyBalanceParallel.m to speed up
% computations. Confidence intervals are bootstrapped.
%
%
%

% Usage:
% energyBalance(name, Ta_pen_surface_initial, heat_power_surface_initial, ...
%                          weights, CI_interval, N_random_samples, printText)
%
% Input:
% name: name of the file with the ensemble prediction without ".mat"
% Ta_pen_surface_initial: Input temperature of the ensemble prediction.
% heat_power_surface: Input supplemental heat of the ensemble prediction.
% weights: vector with weights of the animal to do the energy balance.
% CI_interval: Percentage of the confidence interval. E.g., 0.95.
% N_random_samples: Number of random samples to do the bootstrap confidence interval calculation
% printText: print text for partial calculations
%                  
% 
% Output:
% None.


folder_ensembleData = locateFolderPath("ensembleData");
load([folder_ensembleData, "/KM/", name, ".mat"])
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

% disping stuff as they are computed
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  pkg load parallel
  more off
else
  % we are executing in MATLAB
end

% These values need to be updated accordingly to what was used in ensemblesTestingMSE
% If heat_power_in had been saved, that wouldn't been necessary
Ta_pen_surface =             repmat(Ta_pen_surface_initial, 1, length(heat_power_surface_initial));
heat_power_surface = reshape(repmat(heat_power_surface_initial,  length(Ta_pen_surface_initial), 1), 1, []);



% matrix of differences. Mean, lower x% CI, and upper x% CI
mean_mat_diff = zeros(length(weights), size(predictionsEnsemble, 1));
median_mat_diff = zeros(length(weights), size(predictionsEnsemble, 1));
lower_mat_diff = zeros(length(weights), size(predictionsEnsemble, 1));
upper_mat_diff = zeros(length(weights), size(predictionsEnsemble, 1));

% uncertainties
A_random        = normrnd(1, 0.05, 1, N_random_samples); % multiplies the calculated area
Ta_ideal_random = normrnd(0, 1,    1, N_random_samples); % adds to the calculated Ta_ideal
THP_random      = normrnd(0, 0.15, 1, N_random_samples); % adds to the calculated THP
LHP_random      = unifrnd(0, 15,   1, N_random_samples); % LHP values (W/m2). See energyBalanceParrallel for a discussion.
q_skin_random   = normrnd(0, 1,    1, N_random_samples); % unitary std. Gets multiplied by the q_skin_std


time_i = time;
fun = @(weights) energyBalanceParallel(weights, predictionsEnsemble, ...
                 Ta_pen_surface, heat_power_surface, CI_interval, ...
                 N_random_samples, A_random, Ta_ideal_random, THP_random, LHP_random, ...
                 q_skin_random);
                 

% running
[lower_mat_diff, mean_mat_diff, median_mat_diff, upper_mat_diff] = pararrayfun(8, fun, weights);
lower_mat_diff  = reshape( reshape(lower_mat_diff, length(Ta_pen_surface), [])', ...
                         length(weights), length(Ta_pen_surface_initial), []);
mean_mat_diff   = reshape( reshape(mean_mat_diff,  length(Ta_pen_surface), [])', ...
                         length(weights), length(Ta_pen_surface_initial), []);
median_mat_diff = reshape( reshape(median_mat_diff,  length(Ta_pen_surface), [])', ...
                         length(weights), length(Ta_pen_surface_initial), []);
upper_mat_diff  = reshape( reshape(upper_mat_diff, length(Ta_pen_surface), [])', ...
                         length(weights), length(Ta_pen_surface_initial), []);

lower_mat_diff_temp = abs(lower_mat_diff);
[max_val_lower, max_pos_lower] = min(lower_mat_diff_temp, [], 3);

mean_mat_diff_temp = abs(mean_mat_diff);
[max_val_mean, max_pos_mean] = min(mean_mat_diff_temp, [], 3);

median_mat_diff_temp = abs(median_mat_diff);
[max_val_median, max_pos_median] = min(median_mat_diff_temp, [], 3);

upper_mat_diff_temp = abs(upper_mat_diff);
[max_val_upper, max_pos_upper] = min(upper_mat_diff_temp, [], 3);


min_val_lower_real    = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));
min_power_lower_real  = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));

min_val_mean_real     = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));
min_power_mean_real   = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));

min_val_median_real   = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));
min_power_median_real = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));

min_val_upper_real    = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));
min_power_upper_real  = zeros(size(lower_mat_diff, 1), size(lower_mat_diff, 2));
for i1 = 1:size(lower_mat_diff, 1)
  for i2 = 1:size(lower_mat_diff, 2)
    min_val_lower_real(i1, i2)    = lower_mat_diff(i1, i2, max_pos_lower(i1, i2));
    min_power_lower_real(i1, i2)  = heat_power_surface_initial(max_pos_lower(i1, i2));
    
    min_val_mean_real(i1, i2)     = mean_mat_diff(i1, i2, max_pos_mean(i1, i2));
    min_power_mean_real(i1, i2)   = heat_power_surface_initial(max_pos_mean(i1, i2));
    
    min_val_median_real(i1, i2)   = median_mat_diff(i1, i2, max_pos_median(i1, i2));
    min_power_median_real(i1, i2) = heat_power_surface_initial(max_pos_median(i1, i2));
    
    min_val_upper_real(i1, i2)    = upper_mat_diff(i1, i2, max_pos_upper(i1, i2));
    min_power_upper_real(i1, i2)  = heat_power_surface_initial(max_pos_upper(i1, i2));
  end
end

time_run = time - time_i;
save([folder_ensembleData, "/KM/", name, "_EnergyFluxBalance.mat"], "-V7", "mean_mat_diff", ...
     "median_mat_diff", "lower_mat_diff", "upper_mat_diff", "time_run", "weights", "Ta_pen_surface_initial", ...
     "heat_power_surface_initial", "max_val_lower", "max_pos_lower", "max_val_mean", "max_pos_mean", ...
     "max_val_median", "max_pos_median", "max_val_upper", "max_pos_upper", "min_val_lower_real", ...
     "min_power_lower_real", "min_val_mean_real", "min_power_mean_real", "min_val_median_real", ...
     "min_power_median_real", "min_val_upper_real", "min_power_upper_real");
end