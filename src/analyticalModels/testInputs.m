function [loadMeasuredData, data, data_out, mean_errors, predictionsSet, positions] = ...
                      testInputs(set_i, heat_power_in, Ta_pen_in, ...
                        Ta_brooder_in, Tg_brooder_in, Tr_in, positions)
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
% File:   testInputs.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 27, 2018.
%
%
% Function description:
% Will test if the inputs are valid and, if not, will throw an error.
%
%
% Usage:
% [loadMeasuredData, data, data_out, mean_errors, predictionsSet] = ...
%                       testInputs(set_i, heat_power_in, Ta_pen_in, ...
%                         Ta_brooder_in, Tg_brooder_in, Tr_in, positions)
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
% 
% Output:
% loadMeasuredData: if 1, it had load the measured data. If zero, we are running
%                   for the given inputs.
% data: the measured data, if loaded, or NA.
% data_out: matrix that will store calculation data for the measured data, if 
%           loaded, or NA.
% mean_errors: matrix that will store calculation data for the measured data, if 
%              loaded, or NA.
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
% positions: Vector with positions of the measured dataset to use OR a list with
%            the positions to run the calculations
%


%% Checking if the input conditions are met
loadMeasuredData = 0; % 0: don't load measured data. 1: load measured data
% checking conditions for heat_power
if ( isnan(heat_power_in) )
  if ( isnan(positions) )
    error(["One or all positions of heat_power and positions have NA.\n" ...
           "heat_power or positions should not have NA."]);
  else
    loadMeasuredData = 1;
  end
else
  % we will use the variable positions to get the input data. That is, if we are
  % loading the measured data, positions contains the position of the measured data
  % to load. Otherwise, it is just a loop variable.
  positions = 1:length(heat_power_in);
end

% checking conditions for Ta_pen
if ( isnan(Ta_pen_in) & isnan(Ta_brooder_in) )
  % will be true when a position of both, Ta_pen_in and Ta_brooder_in is NA
  if (loadMeasuredData == 0)
    error(["Ta_pen and Ta_brooder have a position where both are NA " ...
           "and the measured data should not be load (positions is different to NA)"]);
  end
end


if (loadMeasuredData)
  % Including the folder paths in the search directories
  
  folder_datasets = locateFolderPath("datasets");
  addpath(folder_datasets);
    
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


  mean_errors = zeros(1, 4);
  %  1: set number
  %  2: Mean squared error for skin surface temperature
  %  3: Mean squared error for hair-coat surface temperature
  %  4: Mean squared error for skin surface and hair-coat surface temperature
  mean_errors(:) = NA;
  
  % cleaning after myself
  rmpath(folder_datasets);
  
else
  data = NA;
  data_out = NA;
  mean_errors = NA;
  
end

predictionsSet = zeros(length(positions), 10);
%  1: skin-surface temperature (Ts)
%  2: skin-surface heat flux (q_skin)
%  3: hair-coat surface temperature (Th)
%  4: time it took to solve for this position
%  5: value of the objective function
%  6: maximum value of the constraint function
%  7: Ta_pen used in the calculation
%  8: Ta_brooder used in the calculation
%  9: Tg_brooder used in the calculation
% 10: Tr used in the calculation.

end
