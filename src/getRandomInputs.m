function [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
  d, epsilong, dg, z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, ...
  NSkinLayers, NHairCoatLayers, To, qo, ua, ...
  Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getRandomInputs()
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
% File:   getRandomInputs.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on March 31, 2018.
%
%
% Function description:
% Contains the input parameters and its calculations for implementing the 
% analytical model. Outputs the variables required for the computations. 
% Modify this function or create another one if you want to solve for another problem.
%
%
%
% Usage:
% [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
%   d, epsilong, dg, z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, ...
%   NSkinLayers, NHairCoatLayers, To, qo, ua, ...
%   Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getRandomInputs()
%
% Input:
% None.
%
% Output:
% k: vector column with thermal conductivity values for layers (W/(m oC))
% kh: vector column with hair-coat thermal conductivity values for layers (W/(m oC)).
%     kh is used to calculate keff (effective hair-coat thermal conductivity), which
%     depends on air temperature.
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% rhob: vector column with blood density values for layers (kg/m2)
% cb: vector column with blood specific heat values for layers (J/(kg oC))
% Tb_m: vector column with the multipliers for blood temperature values for layers (oC)
% qtotal: vector column with total heat source values for layers (W/m3)
% L: vector column with length of layers (m)
% N: #hairs/m2
% D: diameter of hairs (m)
% HL: lenght of hairs (m)
% epsilon: emissivity of the animal surface
% d: diameter of the animal (m)
% epsilong: black globe emissivity
% dg: emissivity of the black globe
% z: elevation of the location (m).
% Lt: latitute of the location (decimal degrees).
% TMR_m: multiplicative factor for TMR (mean radiant temperature) to consider 
%        for possible measurement errors. TMR is not calculated here because it
%        depends on Ta, Tg, and ua
% h_m: multiplicative factor for h to consider for possible measurement errors
% h_m_skin: multiplicative factor for h to consider for possible measurement errors
%           for convection heat flux at the skin.
% omega: proportion of convection heat transfer at the hair-coat.
% phi: additional proportion of convection heat transfer at the skin.
% NMuscleLayers: number of sublayers of the muscle layer
% NFatLayers: number of sublayers of the muscle layer
% NSkinLayers: number of sublayers of the muscle layer
% NHairCoatLayers: number of sublayers of the muscle layer
% To: vector of offsets temperatures (oC)
% qo: vector of offsets heat fluxes (W/m2)
% ua: air velocity (m/s)
% Ta_pen_stderr: standard error multiply for including uncertainty into Ta_pen.
% Ta_brooder_stderr: standard error multiply for including uncertainty into Ta_brooder.
% Tg_brooder_stderr: standard error multiply for including uncertainty into Tg_brooder.
% Tr_stderr: standard error multiply for including uncertainty into Tr.

% Calculations:

% data of the location (Jaboticabal, Brazil)
z = 595; % elevation (m)
Lt = 21.2583; % latitude (degrees)

% number of sublayers of the layers
NMuscleLayers = round(chi2rnd(10) + 1); %number of sublayers
NFatLayers = round(chi2rnd(10) + 1); %number of sublayers
NSkinLayers = round(chi2rnd(10) + 1); %number of sublayers
NHairCoatLayers = round(chi2rnd(20) + 1); %number of sublayers

NTissueLayers = NMuscleLayers + NFatLayers + NSkinLayers;
NAllLayers = NTissueLayers + NHairCoatLayers;

% length of layers
while( (Lm = normrnd(2, 1)*1e-3) <= 0) end; % length of muscle layer (m). Only values greater than zero are acceptable
while( (Lf = normrnd(5, 2.5)*1e-3) <= 0) end; % length of fat layer (m). Only values greater than zero are acceptable
while( (Ls = normrnd(1, 0.7)*1e-3) <= 0) end; % length of skin layer (m). Only values greater than zero are acceptable
while( (Lh = normrnd(2, 1)*1e-3) <= 0) end; % length of hair coat layer (m). Only values greater than zero are acceptable

% vector of length of layers
L = [Lm/NMuscleLayers*ones(1,NMuscleLayers) ...
     Lf/NFatLayers*ones(1,NFatLayers) ...
     Ls/NSkinLayers*ones(1,NSkinLayers) ...
     Lh/NHairCoatLayers*ones(1,NHairCoatLayers)]; % length of layers (m)
     
% Obtaining physical properties of the hair coat
while( (N = normrnd(150, 50)*1e4) <= 0) end; % #hairs/m2. Only values greater than zero are acceptable
while( (D = normrnd(95, 20)*1e-6) <= 0) end; % hair diameter (m). Only values greater than zero are acceptable
while( (HL = normrnd(12, 2.5)*1e-3) <= 0) end; % hair length (m). Only values greater than zero are acceptableZz
N = N*ones(1, NHairCoatLayers);
D = D*ones(1, NHairCoatLayers);
HL = HL*ones(1, NHairCoatLayers);


% Thermal conductivity
km = zeros(1, NMuscleLayers);
for i1 = 1:NMuscleLayers
  while( (km(i1) = normrnd(0.53, 0.1)) < 0) end;
  % W/mK. Only values greater than zero are acceptable
end
kf = zeros(1, NFatLayers);
for i1 = 1:NFatLayers
  while( (kf(i1) = normrnd(0.29, 0.12)) < 0) end;
  % W/mK. Only values greater than zero are acceptable
end
ks = zeros(1, NSkinLayers);
for i1 = 1:NSkinLayers
  while( (ks(i1) = normrnd(0.21, 0.05)) < 0) end;
  % W/mK. Only values greater than zero are acceptable
end
kh = zeros(1, NHairCoatLayers);
for i1 = 1:NHairCoatLayers
  while( (kh(i1) = normrnd(0.6, 0.35)) < 0) end;
  % W/mK. Only values greater than zero are acceptable
end
% since keff depends on environmental conditions, we don't calculate it here.
k = [km kf ks zeros(1, NHairCoatLayers)]; % thermal conductivity (W/m2)hD Dissertation. Paper Davis and Birkerbak (1972)

% Blood perfusion
wbm = zeros(1, NMuscleLayers);
for i1 = 1:NMuscleLayers
  while( (wbm(i1) = normrnd(2.5, 0.9)* 1e-3) <= 0) end;
end
wbf = zeros(1, NFatLayers);
for i1 = 1:NFatLayers
  while( (wbf(i1) = normrnd(2.5, 1.34)* 1e-3) <= 0) end;
end
wbs = zeros(1, NSkinLayers);
for i1 = 1:NSkinLayers
  while( (wbs(i1) = normrnd(1.55, 0.5)* 1e-3) <= 0) end;
end
wb = [wbm wbf wbs zeros(1,NHairCoatLayers)]; % blood perfusion (m3/(m3 s))

% Blood specific heat and density
% cb is being assigned the value of the volumetric heat capacity, which is cb*rhob. Hence rhob = 1
cb = zeros(1, NAllLayers); % blood specific heat ((J/(kg oC))
rhob = [ones(1,NTissueLayers) zeros(1,NHairCoatLayers)]; % blood density (kg/m3)

% Only values greater than zero are acceptable
for i1 = 1:(NTissueLayers)
  while( (cb(i1) = 1e6*normrnd(4.1, 0.26)) <= 0) end
end

% -) Blood temperature multiplicative factor; accumulative
Tb_m = zeros(1, NAllLayers);
Tb_m(1) = unifrnd(0.98, 1.01);
for i1 = 2:(NMuscleLayers + NFatLayers + NSkinLayers)
  Tb_m(i1) = Tb_m(i1 - 1)*unifrnd(0.98, 1.01);
end

% metabolic heat generation
qm = zeros(1, NMuscleLayers);
for i1 = 1:NMuscleLayers
  while( (qm(i1) = normrnd(684, 137)) <= 0) end; % W/m3. Only values greater than zero are acceptable
end
qf = zeros(1, NFatLayers);
for i1 = 1:NFatLayers
  while( (qf(i1) = normrnd(368, 74)) <= 0) end; % W/m3. Only values greater than zero are acceptable
end
qs = zeros(1, NSkinLayers);
for i1 = 1:NSkinLayers
  while( (qs(i1) = normrnd(368, 74)) <= 0) end; % W/m3. Only values greater than zero are acceptable
end

qtotal = [qm qf qs zeros(1, NHairCoatLayers)]; % Volumetric heat source (W/m3)

% emissivity of the animal surface
epsilon = unifrnd(0.95, 1);

% animal diameter (m)
while( (d = normrnd(10.5, 3.5)*1e-2) < 0) end; % Only values greater than zero are acceptable

% emissivity of the black globe
epsilong = unifrnd(0.92, 0.99);

% black globe diameter (m)
dg = 0.15;

% To consider for possible experimental measurement errors
while( (TMR_m = normrnd(1, 0.1)) <= 0) end; % Only values greater than zero are acceptable
while( (h_m = normrnd(1, 0.1)) <= 0) end; % Only values greater than zero are acceptable
while( (h_m_skin = normrnd(1, 0.1)) <= 0) end; % Only values greater than zero are acceptable
% to consider for possible model deviations
omega = unifrnd(0, 1);
phi = unifrnd(0, omega);

% vector of offsets temperatures and heat fluxes
To = zeros(1, NMuscleLayers + NFatLayers + NSkinLayers + NHairCoatLayers + 1);
qo = zeros(1, NMuscleLayers + NFatLayers + NSkinLayers + NHairCoatLayers + 1);

% Air velocity
while( (ua = normrnd(0.5, 0.25)) <= 0) end; % Only values greater than zero are acceptable

% Standard error multiply for the predicted temperatures. This values multiply
% the standard error of the predictions or the uncertainty of the measurement
% Standard error multiply for the air temperature measured in the pen
Ta_pen_stderr = normrnd(0, 1);
% Standard error multiply for the predicted air temperature inside the brooder
Ta_brooder_stderr = normrnd(0, 1);
% Standard error multiply for the predicted black globe temperature inside the brooder
Tg_brooder_stderr = normrnd(0, 1);
% Standard error multiply for the predicted rectal temperature
Tr_stderr = normrnd(0, 1);

end
