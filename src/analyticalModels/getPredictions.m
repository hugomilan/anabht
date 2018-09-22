function [heat_power, Ta_pen, Ta_brooder, Tg_brooder, Tr] = ...
            getPredictions(set_i, sample_i, position_i, loadMeasuredData, ...
              data, heat_power_in, Ta_pen_in, Ta_brooder_in, Tg_brooder_in, ...
              Tr_in, predictTaTgTr, Ta_pen_stderr, Ta_brooder_stderr, ...
              Tg_brooder_stderr, Tr_stderr, Ta_brooder_in_increase = 0, ...
              T_measured_uncertainty = 0.35)
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
% File:   getPredictions.m
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
% [heat_power, Ta_pen, Ta_brooder, Tg_brooder, Tr] = ...
%             getPredictions(set_i, sample_i, position_i, loadMeasuredData, ...
%               data, heat_power_in, Ta_pen_in, Ta_brooder_in, Tg_brooder_in, ...
%               Tr_in, predictTaTgTr, Ta_pen_stderr, Ta_brooder_stderr, ...
%               Tg_brooder_stderr, Tr_stderr, Ta_brooder_in_increase, ...
%               T_measured_uncertainty)
% 
%
% Input:
% set_i: #set
% sample_i: #sample
% position_i: #position
% loadMeasuredData: flag that indicates if the measured data was loaded (1).
% data: measured data (if loadMeasuredData == 1) or NA.
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
% predictTaTgTr: variable with the structures to make predictions
% Ta_pen_stderr: standard error multiply for including uncertainty into Ta_pen.
% Ta_brooder_stderr: standard error multiply for including uncertainty into Ta_brooder.
% Tg_brooder_stderr: standard error multiply for including uncertainty into Tg_brooder.
% Tr_stderr: standard error multiply for including uncertainty into Tr.
% Ta_brooder_in_increase: when Ta_brooder is an input, should it be add by the 
%                         the effect of heat_source_power?
%                           0: Ta_brooder is as the input.
%                           1: Ta_brooder is the input + the added effect of
%                              heat_soure_power in increasing Ta_brooder
% T_measured_uncertainty: value of the uncertainty in the measured temperatures.
%                         If a single value, the same uncertainty is used for all
%                         measured temperatures. Otherwise, if a vector, will use
%                         positions [1] for Ta_pen, [2] for Ta_brooder, 
%                         [3] for Tg_brooder, and [4] for Tr
%
%
% 
% Output:
% heat_power: intensity of the supplemental heat source (W). This can be the input
%             intensity or the measured data.
% Ta_pen: predicted air temperature inside the pen (oC). This can be the input 
%         temperature plus an offset or the measured data.
% Ta_brooder: predicted air temperature inside the brooder (oC).This can be the input 
%         temperature plus an offset or the measured data.
% Tg_brooder: predicted black globe temperature (oC).This can be the input 
%         temperature plus an offset or the measured data.
% Tr: predicted rectal temperature (oC).This can be the input 
%         temperature plus an offset or the measured data.
%

% Testing if the vector of uncertainty is at the right size
if (length(T_measured_uncertainty) == 1)
  T_measured_uncertainty *= ones(1,4);
elseif (length(T_measured_uncertainty) != 4)
  error(["The vector with uncertainty in the measured temperature values "...
         "(T_measured_uncertainty) must be size 1 or 4. The given vector has size "
         num2str(length(T_measured_uncertainty)) "."])
end

% if we are loading measured data
if (loadMeasuredData)
  % loading the environmental variables
  Ta_pen = data(sample_i, 13) + Ta_pen_stderr*T_measured_uncertainty(1); % Air temperature in the pen (oC)
  heat_power = data(sample_i, 5); % heat source power (W)
  
  % obtaining the predicted temperatures inside the brooder and the rectal temperature
  if (set_i == 0)
    Ta_brooder = data(sample_i, 10); % Air temperature inside the brooder(oC)
    Tg_brooder = data(sample_i, 12); % Black globe temperature inside the brooder(oC)
    Tr = data(sample_i,  8); % Rectal Temperature (oC)
  else
    % assume that we are only given Ta_pen
    
    % predicting Ta_brooder
    Ta_brooder_x = [1 Ta_pen heat_power]';
    Ta_brooder_sigma = predictTaTgTr.LM_Ta_brooder.sigma*sqrt( 1 + ...
                       Ta_brooder_x'*predictTaTgTr.LM_Ta_brooder.invXX*Ta_brooder_x );
    Ta_brooder = predictTaTgTr.LM_Ta_brooder.coeff*Ta_brooder_x + Ta_brooder_stderr*Ta_brooder_sigma;
    
    % predicting Tg_brooder
    Tg_brooder_x = [1 Ta_pen heat_power]';
    Tg_brooder_sigma = predictTaTgTr.LM_Tg_brooder_Ta_pen.sigma*sqrt( 1 + ...
                       Tg_brooder_x'*predictTaTgTr.LM_Tg_brooder_Ta_pen.invXX*Tg_brooder_x );
    Tg_brooder = predictTaTgTr.LM_Tg_brooder_Ta_pen.coeff*Tg_brooder_x + Tg_brooder_stderr*Tg_brooder_sigma;
    
    % predicting Tr
    Tr = predictTaTgTr.LM_Tr.coeff + Tr_stderr*predictTaTgTr.LM_Tr.sigma;
    
  end   
    
% we are running for given input values
else
  % assuming heat_power_in was given for all positions (no NAs)
  heat_power = heat_power_in(position_i);
  Ta_pen = NA; % initial value of Ta_pen, which is changed if needed
    
  % testing if we need to estimate Ta_brooder
  estimate_Ta_brooder = 1;
  % checking if Ta_brooder_in was given or if we have to estimate it
  if ( position_i <= length(Ta_brooder_in) )
    % we have a vector
    if ( !isnan(Ta_brooder_in(position_i)) )
      % It is not NA, we use the given value
      
      % do we have to add the effect of heat power on Ta_brooder? Is the input 
      % heat power different than zero?
      if (Ta_brooder_in_increase && heat_power != 0)
        % yes, we do; yes, the input is different than zero.
        % sigma from the estimation
        Ta_brooder_sigma = predictTaTgTr.LM_Ta_brooder.sigma*sqrt( 1 + ...
                             heat_power'*predictTaTgTr.LM_Ta_brooder_inc.invXX*heat_power );
        % adding sigma from uncertainty
        Ta_brooder_sigma = sqrt(Ta_brooder_sigma^2 + T_measured_uncertainty(2)^2);
        
        Ta_brooder = Ta_brooder_in(position_i) + ...
                     predictTaTgTr.LM_Ta_brooder_inc.coeff*heat_power + ...
                     Ta_brooder_stderr*Ta_brooder_sigma;
      
    else
      % no, we use the given value plus the uncertainty
        Ta_brooder = Ta_brooder_in(position_i) + Ta_brooder_stderr*T_measured_uncertainty(2);
        
      end
      estimate_Ta_brooder = 0;
    end
  end
    
  if (estimate_Ta_brooder)
    % we have to estimate Ta_brooder. For this, we need Ta_pen
    % checking if Ta_pen_in was given as a vector
    if ( position_i <= length(Ta_pen_in) )
      % we checked before hand. If Ta_brooder(i) is NA, Ta_pen(i) is not
      Ta_pen = Ta_pen_in(position_i) + Ta_pen_stderr*T_measured_uncertainty(1); % Air temperature in the pen (oC)
      % predicting Ta_brooder
      Ta_brooder_x = [1 Ta_pen heat_power]';
      Ta_brooder_sigma = predictTaTgTr.LM_Ta_brooder.sigma*sqrt( 1 + ...
                         Ta_brooder_x'*predictTaTgTr.LM_Ta_brooder.invXX*Ta_brooder_x );
      Ta_brooder = predictTaTgTr.LM_Ta_brooder.coeff*Ta_brooder_x + Ta_brooder_stderr*Ta_brooder_sigma;
    else
      if ( position_i > length(Ta_brooder_in) )
        error(["Solving for position " num2str(position_i) " but Ta_brooder_in " ...
               "has length " num2str(length(Ta_brooder_in)) " and the length of " ...
               "Ta_pen_in is " num2str(length(Ta_pen_in)) ". One of them needs " ...
               "to be specified at position " num2str(position_i) "."])
      else
        error(["Solving for position " num2str(position_i) " but Ta_brooder_in " ...
               "at this position is NA and Ta_pen_in has length " num2str(length(Ta_pen_in)) ...
               ". One of them needs to be specified at position " num2str(position_i) "."])
      end
    end
  end
      
  % testing if we need to estimate Tg_brooder
  estimate_Tg_brooder = 1;
  % checking if Tg_brooder_in was given or if we have to estimate it
  if ( position_i <= length(Tg_brooder_in) )
    % we have a vector
    if ( !isnan(Tg_brooder_in(position_i)) )
      % It is not NA, we use the given value
      Tg_brooder = Tg_brooder_in(position_i) + Tg_brooder_stderr*T_measured_uncertainty(3);
      estimate_Tg_brooder = 0;
    end
  end
    
  if (estimate_Tg_brooder)
    % we have to estimate Tg_brooder
    if (estimate_Ta_brooder)
      % Ta_brooder was estimated as well. In that case, we use Ta_pen to estimate
      % Tg_brooder
      Tg_brooder_x = [1 Ta_pen heat_power]';
      Tg_brooder_sigma = predictTaTgTr.LM_Tg_brooder_Ta_pen.sigma*sqrt( 1 + ...
                       Tg_brooder_x'*predictTaTgTr.LM_Tg_brooder_Ta_pen.invXX*Tg_brooder_x );
      Tg_brooder = predictTaTgTr.LM_Tg_brooder_Ta_pen.coeff*Tg_brooder_x + Tg_brooder_stderr*Tg_brooder_sigma;
    else
      % Ta_brooder was given. We use it to estimate Tg_brooder
      Tg_brooder_x = [1 Ta_brooder heat_power]';
      Tg_brooder_sigma = predictTaTgTr.LM_Tg_brooder.sigma*sqrt( 1 + ...
                       Tg_brooder_x'*predictTaTgTr.LM_Tg_brooder.invXX*Tg_brooder_x );
      Tg_brooder = predictTaTgTr.LM_Tg_brooder.coeff*Tg_brooder_x + Tg_brooder_stderr*Tg_brooder_sigma;
    end
  end
      
  % testing if we need to estimate Tr
  estimate_Tr = 1;
  % checking if Tr_in was given or if we have to estimate it
  if ( position_i <= length(Tr_in) )
    % we have a vector
    if ( !isnan(Tr_in(position_i)) )
      % It is not NA, we use the given value
      Tr = Tr_in(position_i) + Tr_stderr*T_measured_uncertainty(4);
    end
  end
    
  if (estimate_Tr)
    % we have to estimate Tr
    Tr = predictTaTgTr.LM_Tr.coeff + Tr_stderr*predictTaTgTr.LM_Tr.sigma;
  end
    
end
  

end
