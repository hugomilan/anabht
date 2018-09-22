function [matrix_heat heat_offset] = ...
        getHeatFluxParametersProblemIndependent(k, L, qtotal, wb, rhob, cb, qo)
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
% File:   getHeatFluxParametersProblemIndependent.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com

% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrix and vector offset for calculating interface
% heat flux given betas. This function calculates the components that are
% problem independent.
%
%
% Usage:
% [matrix_heat heat_offset] = ...
%        getHeatFluxParametersProblemIndependent(k, L, qtotal, wb, rhob, cb, qo)
%
% Input:
% k: vector column with thermal conductivity values for layers (W/(m oC))
% L: vector column with legth of layers (m)
% qtotal: vector column with total heat source values for layers (W/m3)
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% rhob: vector column with blood density values for layers (kg/m2)
% cb: vector column with blood specific heat values for layers (J/(kg oC))
% qo: vector column with heat fluxes offsets between layers (W/m2)
%
%
% Output:
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
%

% Calculations:

heat_offset = zeros(size(k,2) + 1, 1);
matrix_heat = zeros(size(k,2) + 1, 2*size(k,2));
% First interface
a2 = 1;
if wb(a2) == 0
  matrix_heat(a2, 2*a2 - 1) = -k(a2); % beta_1i
  % matrix_heat(a2, 2*a2 + 0) = 0; % beta_2i
  % heat_offset(a2) = 0;
  
else % wb(a2) != 0
  matrix_heat(a2, 2*a2 - 1) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) ); % beta_1i
  % matrix_heat(a2, 2*a2) = 0; % beta_2i
  % heat_offset(a2) = 0;
end

for a2 = 1:(size(k,2) - 1)
  % Calculating at the intersection i using layer i; x_i = L_i
  %
  % divide them by two because the heat fluxes should be the same no matter
  % what function (left or right) I'm using to calculate it
  if wb(a2) == 0 % calculating q as calculated from this layer
    matrix_heat(a2 + 1, 2*a2 - 1) = -k(a2)/2; % beta_1i
    % matrix_heat(a2 + 1, 2*a2 + 0) = 0; % beta_2i
    % matrix_heat(a2 + 1, 2*a2 + 1) = 0; % beta_1(i+1)
    % matrix_heat(a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    heat_offset(a2 + 1) = qtotal(a2)*L(a2)/2;
    
  else % wb(a2) != 0 % calculating q as calculated from this layer
    matrix_heat(a2 + 1, 2*a2 - 1) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) )/2; % beta_1i
    matrix_heat(a2 + 1, 2*a2 + 0) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) )/2; % beta_2i
    % matrix_heat(a2 + 1, 2*a2 + 1) = 0; % beta_1(i+1)
    % matrix_heat(a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    % heat_offset(a2 + 1) = 0;
  end
  
  % Calculating at the intersection i using layer i + 1; x_(i + 1) = 0
  if wb(a2 + 1) == 0 % calculating q as calculated from the next layer
    % matrix_heat(a2 + 1, 2*a2 - 1) += 0; % beta_1i
    % matrix_heat(a2 + 1, 2*a2 + 0) += 0; % beta_2i
    matrix_heat(a2 + 1, 2*a2 + 1) += -k(a2 + 1)/2; % beta_1(i+1)
    % matrix_heat(a2 + 1, 2*a2 + 2) += 0; % beta_2(i+1)
    heat_offset(a2 + 1) += qo(a2 + 1)/2;
    
  else % wb(a2 + 1) != 0 % calculating q as calculated from the next layer
    % matrix_heat(a2 + 1, 2*a2 - 1) += 0; % beta_1i
    % matrix_heat(a2 + 1, 2*a2 + 0) += 0; % beta_2i
    matrix_heat(a2 + 1, 2*a2 + 1) += -sqrt( k(a2+1)*wb(a2+1)*rhob(a2+1)*cb(a2+1) )/2; % beta_1(i+1)
    % matrix_heat(a2 + 1, 2*a2 + 2) += 0; % beta_2(i+1)
    heat_offset(a2 + 1) += qo(a2 + 1)/2;
    
  end
end

% last heat flux of last element
a2 = size(k,2);
if wb(a2) == 0
  matrix_heat(a2 + 1, 2*a2 - 1) = -k(a2); % beta_1i
  % matrix_heat(a2 + 1, 2*a2 + 0) = 0; % beta_2i
  heat_offset(a2 + 1) = qtotal(a2)*L(a2) + qo(a2 + 1);
  
else % wb(a2) != 0
  matrix_heat(a2 + 1, 2*a2 - 1) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
  matrix_heat(a2 + 1, 2*a2 + 0) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
  heat_offset(a2 + 1) = qo(a2 + 1);
  
end
