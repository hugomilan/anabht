function [A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin] = ...
            getProblemIndependentMatricesRadiationKM(alpha_eff_KM, deltax, qtotal)
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
% File:   getProblemIndependentMatricesRadiationKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrices componentes for the KM model that are independent of the problem.
%
%
% Usage:
% [A_FD, E, b_FD1, A_q_skin, E_q_skin, b_q_skin] = ...
%            getProblemIndependentMatricesRadiationKM(alpha_eff_KM, deltax, qtotal)
%
% Input:
% alpha_eff_KM: effective absorption coefficient of the hair coat for the KM model
% deltax: space between points (m)
% qtotal: vector column with total heat source values for layers (W/m3)
%
%
% Output:
% A_FD: matrix that multiplies T for the KM model
% E: matrix that multiplies T^4 for the KM model
% b_FD1: matrix that contains two vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% E_q_skin: vector that multiplies T^4 to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% 
  
NHairCoatLayers = length(alpha_eff_KM);
sigma = 5.670373e-8; % Stefan-Boltzmann constant

% A_FD is matrix that multiplies the temperatures
A_FD = zeros( NHairCoatLayers + 1 );
% E is matrix that multiplies the temperatures to the 4th degree
E = zeros( NHairCoatLayers + 1 );
% b is vector to add
b_FD1 = zeros( NHairCoatLayers + 1, 2);

% vectors to obtain q_skin
A_q_skin = zeros( NHairCoatLayers + 1, 1);
E_q_skin = zeros( NHairCoatLayers + 1, 1);
b_q_skin = zeros( 1 );

% changing parameters so that they conform with the mathematical derivation
alpha_eff_KM(2:(end + 1)) = alpha_eff_KM;
% Note that alpha_eff_KM(i) = alpha_eff_KM(j) for all j and i. You will find
% some terms that could be cancelled out in this code but they were maintained
% for compatibility with previous versions of this code.

% k(2:(end + 1)) = k;
qtotal(2:(end + 1)) = qtotal;
for ii = 2:NHairCoatLayers
  % element i and Ts
  E(ii, 1) = 2*sigma*( + alpha_eff_KM( 2)*expintn(2, ( ii - 3/2 )*deltax*alpha_eff_KM(2) ) ...
                       + alpha_eff_KM(ii)*expintn(2, ( ii - 1   )*deltax*alpha_eff_KM(ii) ) ...
                       - alpha_eff_KM( 2)*expintn(2, ( ii - 1   )*deltax*alpha_eff_KM(2) ) );
                           
  % element i and j, with j different from Ts and Th
  for jj = [ 2:(ii - 1), (ii + 1):NHairCoatLayers ]
    E(ii, jj) = 2*sigma*alpha_eff_KM(jj)*( + expintn(2, ( abs(ii - jj) - 1/2 )*deltax*alpha_eff_KM(jj) ) ...
                                           - expintn(2, ( abs(ii - jj)       )*deltax*alpha_eff_KM(jj) ) ...
                                           + expintn(2, ( abs(ii - jj)       )*deltax*alpha_eff_KM(jj + 1) )
                                           - expintn(2, ( abs(ii - jj) + 1/2 )*deltax*alpha_eff_KM(jj + 1) ) );
  end
  % element i and itself
  E(ii, ii) = -2*sigma*( + alpha_eff_KM(ii    )*expintn(2, 1/2*deltax*alpha_eff_KM(ii)) ...
                         + alpha_eff_KM(ii + 1)*expintn(2, 1/2*deltax*alpha_eff_KM(ii + 1)) );

  % element i and Th  
  E(ii, end) = 2*sigma*alpha_eff_KM(end)*( expintn(2, ( NHairCoatLayers - ii + 1/2 )*deltax*alpha_eff_KM(end) ) - ...
                                           expintn(2, ( NHairCoatLayers - ii + 1   )*deltax*alpha_eff_KM(end) ) );

  b_FD1(ii, 1) = 2*sigma*alpha_eff_KM(ii)*expintn(2, (NHairCoatLayers - ii + 1)*deltax*alpha_eff_KM(ii) ); % *TMRK^4
  b_FD1(ii, 2) = (qtotal(ii) + qtotal(ii + 1))/2; % As I'm solving at a point, I average the heat sources from both mediums
  % should do b_FD(ii) = b_FD1(ii, 1)*TMRK^4 + b_FD1(ii, 2);

  A_FD(ii, (ii - 1):(ii + 1) ) = 1/deltax^2*[1 -2 1]; % *k(ii)
  
end

% setting up first boundary conditions
A_FD(1,1) = 1;
% b_FD(1) = -TsK; % TsK is dynamically calculated

% setting up last boundary conditions
% A_FD(end, end) = -k(end)/deltax - h; % h is dynamically calculated--depends on Th
A_FD(end, end - 1) = 1/deltax; % *k(end)
% b_FD(end) = h*TaK; % h is dynamically calculated--depends on Th

% Obtaining the vectors for calculating q_skin
A_q_skin(1:2) = 1/deltax*[1 -1]; % *k(1)
b_q_skin = - 2*sigma*expintn(3, alpha_eff_KM(2)*NHairCoatLayers*deltax); % *TMRK^4
E_q_skin(1) = 2*sigma*expintn(3, alpha_eff_KM(2)*deltax/2);
for ii = 2:NHairCoatLayers
  E_q_skin(ii) = - 2*sigma*( + expintn(3, (ii - 1/2)*deltax*alpha_eff_KM(ii) ) ...
                             - expintn(3, (ii      )*deltax*alpha_eff_KM(ii) ) ...
                             + expintn(3, (ii      )*deltax*alpha_eff_KM(ii + 1) ) ...
                             - expintn(3, (ii + 1/2)*deltax*alpha_eff_KM(ii + 1) ) );
end
% for Th
E_q_skin(end) = - 2*sigma*( + expintn(3, (NHairCoatLayers - 1/2)*deltax*alpha_eff_KM(end) ) ...
                            - expintn(3, (NHairCoatLayers      )*deltax*alpha_eff_KM(end) ) );
end
