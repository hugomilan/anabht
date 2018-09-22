function [lower_diff, mean_diff, median_diff, upper_diff] = energyBalanceParallel(weight, ...
           predictionsEnsemble, Ta_pen_surface, heat_power_surface, CI_interval, N_random_samples, ...
           A_random, Ta_ideal_random, THP_random, LHP_random, q_skin_random, printText = false)
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
% File:   energyBalanceParallel.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
% Created on June 23, 2018.
%
%
% Function description:
% Calculates the energy (in)balance between heat flowing to the environment and
% and generated heat. Confidence interval is bootstrapped.
%
%
%
%
% Usage:
% [lower_diff, mean_diff, median_diff, upper_diff] = energyBalanceParallel(weight, ...
%           predictionsEnsemble, Ta_pen_surface, heat_power_surface, CI_interval, N_random_samples, ...
%           A_random, Ta_ideal_random, THP_random, LHP_random, q_skin_random, printText)
%
% Input:
% weight: the weight of the animal to do the energy balance
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
%
% Ta_pen_surface: Input temperature of the ensemble prediction.
% heat_power_surface: Input supplemental heat of the ensemble prediction.
% CI_interval: Percentage of the confidence interval. E.g., 0.95.
% N_random_samples: Number of random samples to do the bootstrap confidence interval calculation
% A_random: random values to multiply the area for each sample
% Ta_ideal_random: random values of the Ta_ideal to add for each sample
% THP_random: random values of the THP to add for each sample
% LHP_random: random values of the LHP to add for each sample
% q_skin_random: random values of the std of q_skin. To multiply q_skin std and add
%                for each sample.
% printText: print text for partial calculations
%                  
% 
% Output:
% lower_diff: lower value for the confidence interval of the difference between
%             energy flowing out vs. generated energy.
% mean_diff: mean value for the confidence interval of the difference between
%            energy flowing out vs. generated energy.
% median_diff: median value for the confidence interval of the difference between
%              energy flowing out vs. generated energy.
% upper_diff: upper value for the confidence interval of the difference between
%             energy flowing out vs. generated energy.

% matrix of differences. Mean, lower x% CI, and upper x% CI
CI_lower_position = round( (1 - CI_interval)/2*N_random_samples );
mean_diff = zeros(1, size(predictionsEnsemble, 1));
lower_diff = zeros(1, size(predictionsEnsemble, 1));
upper_diff = zeros(1, size(predictionsEnsemble, 1));
vec_diff = zeros(1, N_random_samples); % vector with temporary data

A = 0.087*weight^(2/3); % surface area
A_i2 = A_random*A; % including uncertainty in the surface area

for i1 = 1:size(predictionsEnsemble, 1)
%  Ta         = predictionsEnsemble(i1, 8); % Ta_brooder
  q_skin     = predictionsEnsemble(i1, 2); % q_skin
  q_skin_std = predictionsEnsemble(i1, 12); % q_skin standard deviation
  q_skin_i2  = q_skin_std*q_skin_random + q_skin; % including uncertainty in q_skin
  H          = heat_power_surface(i1);
  
  % obtaining the equivalent temperature. The one that would provide the same heat flux
  % when supplemental heat is zero.
  [min_val, min_pos] = min( abs( q_skin - predictionsEnsemble(1:length(Ta_pen_surface), 2) ) );
  Ta = predictionsEnsemble(min_pos, 8);
  
  THP = 10.^(0.715 - 0.0025*Ta + 0.0211*log10(weight));
  THP_i2    = THP_random + THP; % including uncertainty in THP

%  LHP    = -2.26 + 0.194*Ta + 0.0679*weight - 0.0034*Ta*weight;
%  LHP_i2 = LHP_random + LHP; % including uncertainty in LHP
% this model for
% latent heat production was not used because 1) it considers the total latent heat
% production in the pen, including from urine, feces, water dripping, etc., so it 
% is an upper limit to the total latent heat production of piglets and is likely
% much higher than what actual values would be. 2) no evaluation metrics of this
% equation were shown in the paper (i.e., no R^2 or anything similar).
%
% Hence, latent heat production was assumed to be a value between 0 and 15 W/m2,
% based on values published in the literature for other animal species.
  LHP_i2 = LHP_random/weight.*A_i2;
  
  vec_diff  = q_skin_i2 - (THP_i2 - LHP_i2)*weight./A_i2;
    
  vec_diff = sort(vec_diff);
  mean_diff(i1)  = mean(vec_diff);
  median_diff(i1)  = median(vec_diff);
  lower_diff(i1) = vec_diff(CI_lower_position);
  upper_diff(i1) = vec_diff(N_random_samples - CI_lower_position);
    
  if (printText)
    disp(["Weight = " num2str(weight, "%2.1f") " kg; Ta = " num2str(Ta_pen_surface(i1), "%2.1f") " oC; H = " num2str(H, "%3.0f") " W; " ...
          num2str(CI_interval*100, "%2.2f") "%CI: lower " num2str(lower_diff(i1), "%2.2f") ...
          " W; mean = " num2str(mean_diff(i1), "%2.2f") ...
          " W; median = " num2str(median_diff(i1), "%2.2f") " W; upper = " ...
          num2str(upper_diff(i1), "%2.2f") " W."])
  end
end


     
end
