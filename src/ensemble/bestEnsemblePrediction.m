function predictionsEnsemble = bestEnsemblePrediction(heat_power_in = NA, ...
                                 Ta_pen_in = NA, Ta_brooder_in = NA, ...
                                 Tg_brooder_in = NA, Tr_in = NA, training = NA, ...
                                 nameEnsemble = "", nparallel = 3, loadEstimations = 1, ...
                                 Ta_brooder_in_increase = 0)
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
% File:   bestEnsemblePrediction.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on Jun 3, 2018.
%
%
% Function description:
% Computes predictions using the best ensemble as described in the paper

%
%
% Usage:
% predictionsEnsemble = bestEnsemblePrediction(heat_power_in, Ta_pen_in, ...
%                                Ta_brooder_in, Tg_brooder_in, Tr_in, training, ...
%                                nameEnsemble, nparallel, loadEstimations, ...
%                                Ta_brooder_in_increase)
%
% Input:
% heat_power_in: Power of the supplemental heat source (W). Can be a vector. Used to
%                estimate Ta_brooder, Tg_brooder, and Tr. If NA, will check
%                if positions is NA. If so, will throw an error. If not, will load the
%                measured dataset and solve the problem for the given positions.
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
% nameEnsemble: name used to salve the ensemble predictions
% training: Should calculate MSE for training, testing, or both? If 0, calculates 
%           MSE for testing. If 1, calculates MSE for training. If 2, calculates MSE for
%           both, testing and training. Otherwise, will not load any.
% nparallel: Number or process to lunch in parallel.
% loadEstimations: when running for training or testing:
%                       0: recalculates the predictions even if they are available
%                       1: load predictions if they are available
% Ta_brooder_in_increase: when Ta_brooder is an input, should it be add by the 
%                         the effect of heat_source_power?
%                           0: Ta_brooder is as the input.
%                           1: Ta_brooder is the input + the added effect of
%                              heat_soure_power in increasing Ta_brooder
%                  
% 
% Output:
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

%%%% NO NEED TO MODIFY THIS FILE


sets = [7651, 2180, 3410];
predictionsEnsemble = ensemblePrediction( sets, ...
                                 ones(1, length(sets))/length(sets), training, ...
                                 "KM", heat_power_in, Ta_pen_in, ...
                                 Ta_brooder_in, Tg_brooder_in, Tr_in, ...
                                 nameEnsemble, nparallel, loadEstimations, ...
                                 Ta_brooder_in_increase);

