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
% File:   example.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on September 22, 2018.
%
%
% Function description:
% Example on how to run the code

% before you run this code, you will need the following octave packages:
% struct; install it with "pkg install -forge struct"
% parallel; install it with "pkg install -forge parallel"

% adding paths of the algorithms
if exist("src") != 0
  addpath("src");
  addpath("src/ensemble");
end

Ta = 25:30; % temperatures from 25-30 oC
H = zeros(1,6); % 0 W supplemental heat

% if this does not work, try commenting out line 452 of parcellfun.m
% you can use the following command to open it up: open parcellfun.m
% see discussion on https://savannah.gnu.org/bugs/index.php?51448
% 
predictions = bestEnsemblePrediction(H, NA, Ta);
% predictions: matrix with the predictions from the model. Rows are the different
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

plot(Ta, predictions(:,1))
hold on
plot(Ta, predictions(:,1) + predictions(:,11))
plot(Ta, predictions(:,1) - predictions(:,11))
hold off