function temp_offset = updateTempParametersProblemDependent(temp_offset, wb, Tb)
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
% File:   updateTempParametersProblemDependent.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com

% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrix and vector offset for calculating interface
% temperatures given betas for a specific problem.
%
%
% Usage:
% temp_offset = updateTempParametersProblemDependent(temp_offset, wb, Tb)
%
% Input:
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% Tb: vector column with blood temperature values for layers (oC)
%
%
% Output:
% temp_offset: vector containing offset for the relationship between betas and the temperatures
%

% Calculations:

%temp_offset = zeros(size(k,2) + 1, 1);
%matrix_temp = zeros(size(k,2) + 1, 2*size(k,2));
% First element
a2 = 1;
if wb(a2) == 0
  % matrix_temp(a2, 2*a2 - 1) = 0; % beta_1i
%  matrix_temp(a2, 2*a2 + 0) = 1; % beta_2i
  % temp_offset(a2) = 0;
  
else % wb(a2) != 0
  % matrix_temp(a2, 2*a2 - 1) = 0; % beta_1i
%  matrix_temp(a2, 2*a2 + 0) = 1; % beta_2i
  temp_offset(a2) += Tb(a2); % + qtotal(a2)/( wb(a2)*rhob(a2)*cb(a2) );
end

for a2 = 1:(size(wb,2) - 1)
  % Calculating at the intersection i using layer i; x_i = L_i
  %
  % divide them by two because the temperatures should be the same no matter
  % what function I'm using to calculate it
  if wb(a2) == 0 % calculating T as calculated from this layer
%    matrix_temp(a2 + 1, 2*a2 - 1) = L(a2)/2; % beta_1i
%    matrix_temp(a2 + 1, 2*a2 + 0) = 1/2; % beta_2i
    % matrix_temp(a2 + 1, 2*a2 + 1) = 0; % beta_1(i+1)
    % matrix_temp(a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
%    temp_offset(a2 + 1) = - qtotal(a2)/(2*k(a2))*L(a2)^2/2;
    
  else % wb(a2) != 0 % calculating T as calculated from this layer
%    matrix_temp(a2 + 1, 2*a2 - 1) = sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) )/2; % beta_1i
%    matrix_temp(a2 + 1, 2*a2 + 0) = cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) )/2; % beta_2i
    % matrix_temp(a2 + 1, 2*a2 + 1) = 0; % beta_1(i+1)
    % matrix_temp(a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    temp_offset(a2 + 1) += Tb(a2)/2; % + ( qtotal(a2)/( wb(a2)*rhob(a2)*cb(a2) ) )/2;
    
  end
  
  % Calculating at the intersection i using layer i + 1; x_(i + 1) = 0
  if wb(a2 + 1) == 0 % calculating T as calculated from the next layer
    % matrix_temp(a2 + 1, 2*a2 - 1) += 0; % beta_1i
    % matrix_temp(a2 + 1, 2*a2 + 0) += 0; % beta_2i
    % matrix_temp(a2 + 1, 2*a2 + 1) += 0; % beta_1(i+1)
%     matrix_temp(a2 + 1, 2*a2 + 2) += 1/2; % beta_2(i+1)
%    temp_offset(a2 + 1) += To(a2 + 1)/2;
    
  else  % wb(a2 + 1) != 0 % calculating T as calculated from the next layer
    % matrix_temp(a2 + 1, 2*a2 - 1) += 0; % beta_1i
    % matrix_temp(a2 + 1, 2*a2 + 0) += 0; % beta_2i
    % matrix_temp(a2 + 1, 2*a2 + 1) += 0; % beta_1(i+1)
%     matrix_temp(a2 + 1, 2*a2 + 2) += 1/2; % beta_2(i+1)
     temp_offset(a2 + 1) += Tb(a2 + 1)/2; % ( qtotal(a2+1)/( wb(a2+1)*rhob(a2+1)*cb(a2+1) ) )/2 + To(a2 + 1)/2;
    
  end
end

% last temperature of last element
a2 = size(wb,2);
if wb(a2) == 0
%  matrix_temp(a2 + 1, 2*a2 - 1) = L(a2); % beta_1i
%  matrix_temp(a2 + 1, 2*a2 + 0) = 1; % beta_2i
%  temp_offset(a2 + 1) = -qtotal(a2)/( 2*k(a2) )*L(a2)^2 + To(a2 + 1);
  
else % wb(a2) != 0
%  matrix_temp(a2 + 1, 2*a2 - 1) = sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
%  matrix_temp(a2 + 1, 2*a2 + 0) = cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
  temp_offset(a2 + 1) += Tb(a2); % + qtotal(a2)/( wb(a2)*rhob(a2)*cb(a2) ) + To(a2 + 1);
  
end
