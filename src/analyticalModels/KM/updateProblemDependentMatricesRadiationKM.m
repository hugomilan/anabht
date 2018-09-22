function [A_FD, b_FD, A_q_skin, b_q_skin] = ...
            updateProblemDependentMatricesRadiationKM(TMR, k, A_FD, b_FD1, A_q_skin, b_q_skin)
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
% File:   updateProblemDependentMatricesRadiationKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrix and vector offset for the KM model given thermal parameters that
% are problem independent
%
%
% Usage:
% [A_FD, b_FD, A_q_skin, b_q_skin] = ...
%            updateProblemDependentMatricesRadiationKM(TMR, k, A_FD, b_FD1, A_q_skin, b_q_skin)
%
% Input:
% TMR: Mean radiant temperature (oC)
% k: vector column with thermal conductivity values for layers (W/(m oC))
% A_FD: matrix that multiplies T for the KM model
% b_FD1: matrix that contains two vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
%
%
% Output:
% A_FD: matrix that multiplies T for the KM model
% b_FD: matrix that contains the vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% 
    
NHairCoatLayers = size(A_FD, 1) - 1;
TMRK = TMR + 273.15;
%  sigma = 5.670373e-8; % Stefan-Boltzmann constant

% A_FD is matrix that multiplies the temperatures
%  A_FD = zeros( NHairCoatLayers + 1 );
% E is matrix that multiplies the temperatures to the 4th degree
%  E = zeros( NHairCoatLayers + 1 );
% b is vector to add
b_FD = b_FD1(:, 2);
b_FD += b_FD1(:, 1)*TMRK^4;

% vectors to obtain q_skin
%  A_q_skin = zeros( NHairCoatLayers + 1, 1);
%  E_q_skin = zeros( NHairCoatLayers + 1, 1);
%  b_q_skin = zeros( 1 );

% changing parameters so that they conform with the mathematical derivation
%  alpha_eff_KM(2:(end + 1)) = alpha_eff_KM;
k(2:(end + 1)) = k;
%  qtotal(2:(end + 1)) = qtotal;
for ii = 2:NHairCoatLayers
  % element i and Ts
%    E(ii, 1) = 2*sigma*( + alpha_eff_KM( 2)*expintn(2, ( ii - 3/2 )*deltax*alpha_eff_KM(2) ) ...
%                         + alpha_eff_KM(ii)*expintn(2, ( ii - 1   )*deltax*alpha_eff_KM(ii) ) ...
%                         - alpha_eff_KM( 2)*expintn(2, ( ii - 1   )*deltax*alpha_eff_KM(2) ) );
                           
  % element i and j, with j different from Ts and Th
%    for jj = [ 2:(ii - 1), (ii + 1):NHairCoatLayers ]
%      E(ii, jj) = 2*sigma*alpha_eff_KM(jj)*( + expintn(2, ( abs(ii - jj) - 1/2 )*deltax*alpha_eff_KM(jj) ) ...
%                                             - expintn(2, ( abs(ii - jj)       )*deltax*alpha_eff_KM(jj) ) ...
%                                             + expintn(2, ( abs(ii - jj)       )*deltax*alpha_eff_KM(jj + 1) )
%                                             - expintn(2, ( abs(ii - jj) + 1/2 )*deltax*alpha_eff_KM(jj + 1) ) );
%    end
  % element i and itself
%    E(ii, ii) = -2*sigma*( + alpha_eff_KM(ii    )*expintn(2, 1/2*deltax*alpha_eff_KM(ii)) ...
%                           + alpha_eff_KM(ii + 1)*expintn(2, 1/2*deltax*alpha_eff_KM(ii + 1)) );

  % element i and Th  
%    E(ii, end) = 2*sigma*alpha_eff_KM(end)*( expintn(2, ( NHairCoatLayers - ii + 1/2 )*deltax*alpha_eff_KM(end) ) - ...
%                                             expintn(2, ( NHairCoatLayers - ii + 1   )*deltax*alpha_eff_KM(end) ) );
  
  % b_FD(ii) += b_FD1(ii, 1)*TMRK^4; % calculated above
  % b_FD1(ii, 1) = 2*sigma*alpha_eff_KM(ii)*expintn(2, (NHairCoatLayers - ii + 1)*deltax*alpha_eff_KM(ii) ); % *TMRK^4
  % b_FD1(ii, 2) = (qtotal(ii) + qtotal(ii + 1))/2; % As I'm solving at a point, I average the heat sources from both mediums

  A_FD(ii, (ii - 1):(ii + 1) ) *= k(ii); %1/deltax^2*[1 -2 1];
  
end

% setting up first boundary conditions
% A_FD(1,1) = 1;
% b_FD(1) = -TsK; % TsK is dynamically calculated

% setting up last boundary conditions
% A_FD(end, end) = -k(end)/deltax - h; % h is dynamically calculated--depends on Th
A_FD(end, end - 1) *= k(end); %/deltax
% b_FD(end) = h*TaK; % h is dynamically calculated--depends on Th

% Obtaining the vectors for calculating q_skin
A_q_skin(1:2) *= k(1); % 1/deltax*[1 -1];
b_q_skin *= TMRK^4; % *(-2*sigma*expintn(3, alpha_eff_KM(2)*NHairCoatLayers*deltax))
% E_q_skin(1) = 2*sigma*expintn(3, alpha_eff_KM(2)*deltax/2);
%  for ii = 2:NHairCoatLayers
%    E_q_skin(ii) = - 2*sigma*( + expintn(3, (ii - 1/2)*deltax*alpha_eff_KM(ii) ) ...
%                               - expintn(3, (ii      )*deltax*alpha_eff_KM(ii) ) ...
%                               + expintn(3, (ii      )*deltax*alpha_eff_KM(ii + 1) ) ...
%                               - expintn(3, (ii + 1/2)*deltax*alpha_eff_KM(ii + 1) ) );
%  end
% for Th
%  E_q_skin(end) = - 2*sigma*( + expintn(3, (NHairCoatLayers - 1/2)*deltax*alpha_eff_KM(end) ) ...
%                              - expintn(3, (NHairCoatLayers      )*deltax*alpha_eff_KM(end) ) );
end
