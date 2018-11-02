function [k, kh, wb, rhob, cb, Tb_m, qtotal, L, N, D, HL, epsilon, ...
  d, epsilong, dg, ...
  z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
  To, qo, ua, Ta_pen_stderr, Ta_brooder_stderr, Tg_brooder_stderr, Tr_stderr] = getInputs(NHairCoatLayersInput = 20)
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
% File:   getInputs.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 5, 2018.
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
%  rho, alpha, rhoh, alphah, d, epsilong, dg, ...
%  z, Lt, TMR_m, h_m, h_m_skin, omega, phi, NMuscleLayers, NFatLayers, NSkinLayers, NHairCoatLayers, ...
%  To, qo, ua] = getInputs()
%
% Input:
% NHairCoatLayersInput: number of divisions of the hair coat layer.
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
NMuscleLayers = 1; % number of sublayers
NFatLayers = 1; % number of sublayers
NSkinLayers = 1; % number of sublayers
NHairCoatLayers = NHairCoatLayersInput; % number of sublayers

NTissueLayers = NMuscleLayers + NFatLayers + NSkinLayers;
NAllLayers = NTissueLayers + NHairCoatLayers;

% length of layers
Lm = 2e-3; % length of muscle layer (m)
Lf = 5e-3; % length of fat layer (m)
Ls = 1e-3; % length of skin layer (m)
Lh = 2e-3; % length of hair coat layer (m)

% vector of length of layers
L = [Lm/NMuscleLayers*ones(1,NMuscleLayers) ...
     Lf/NFatLayers*ones(1,NFatLayers) ...
     Ls/NSkinLayers*ones(1,NSkinLayers) ...
     Lh/NHairCoatLayers*ones(1,NHairCoatLayers)]; % length of layers (m)
     
% Obtaining physical properties of the hair coat
HL = 12e-3*ones(1, NHairCoatLayers); % hair length (m)
N = 150e4*ones(1, NHairCoatLayers); % #hairs/m2
D = 95e-6*ones(1, NHairCoatLayers); % hair diameter (m)


% Thermal conductivity
km = 0.53*ones(1,NMuscleLayers); % W/mK
kf = 0.29*ones(1,NFatLayers); % W/mK
ks = 0.21*ones(1,NSkinLayers); % W/mK
kh = 0.30*ones(1, NHairCoatLayers); % W/mK
% since keff depends on environmental conditions, we don't calculate it here.
k = [km kf ks zeros(1, NHairCoatLayers)]; % thermal conductivity (W/m2)

% Blood perfusion
wbm = 2.5e-3*ones(1,NMuscleLayers); % m3/(m3*s)
wbf = 2.5e-3*ones(1,NFatLayers); % m3/(m3*s)
wbs = 1.55e-3*ones(1,NSkinLayers); % m3/(m3*s)
wb = [wbm wbf wbs zeros(1,NHairCoatLayers)]; % blood perfusion (m3/(m3 s))

% Blood specific heat and density
% cb is being assigned the value of the volumetric heat capacity, which is cb*rhob. Hence rhob = 1
cb = 4.1e6*[ones(1,NTissueLayers) zeros(1,NHairCoatLayers)]; % blood specific heat ((J/(kg oC))
rhob = [ones(1,NTissueLayers) zeros(1,NHairCoatLayers)]; % Setting volumetric heat capacity on cb. Hence, rhob = 1

% -) Blood temperature multiplicative factor;
Tb_m = 1*ones(1,NAllLayers);

% metabolic heat generation
qm = 684*ones(1,NMuscleLayers); % W/m3
qf = 368*ones(1,NFatLayers); % W/m3
qs = 368*ones(1,NSkinLayers); % W/m3
qtotal = [qm qf qs zeros(1, NHairCoatLayers)]; % Volumetric heat source (W/m3)

% emissivity of the animal surface
epsilon = 0.98;
% animal diameter (m)
d = 10.5e-2; % from data

% emissivity of the black globe
epsilong = 0.95;
% black globe diameter (m)
dg = 0.15;

% To consider for possible experimental measurement errors
TMR_m = 1;
h_m = 1;
h_m_skin = 1;
% to consider for possible model deviations
omega = 1;
phi = 0;

% vector of offsets temperatures and heat fluxes
To = zeros(1, NAllLayers + 1);
qo = zeros(1, NAllLayers + 1);

% Air velocity
ua = 0.5;

% Standard error multiply for the predicted temperatures. This values multiply
% the standard error of the predictions or the uncertainty of the measurement
% Standard error multiply for the air temperature measured in the pen
Ta_pen_stderr = 0;
% Standard error multiply for the predicted air temperature inside the brooder
Ta_brooder_stderr = 0;
% Standard error multiply for the predicted black globe temperature inside the brooder
Tg_brooder_stderr = 0;
% Standard error multiply for the predicted rectal temperature
Tr_stderr = 0;
end
