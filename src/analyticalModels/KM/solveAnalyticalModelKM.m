function [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
   solveAnalyticalModelKM(Ta, k, L, qo, d, z, Lt, h_m, h_m_skin, omega, phi, ua, NTissueLayers, ...
   A, b, A_FD, E, b_FD, A_q_skin, E_q_skin, b_q_skin, ...
   matrix_temp, temp_offset, matrix_heat, heat_offset, ...
   beta_initial, T_initial = 37, tol = sqrt(eps), maxiter = 1e3)
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
% File:   solveAnalyticalModelKM.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 3, 2018.
%
%
% Function description:
% Given the inputs that describe the environment and the problem, this function
% solves the analytical problem.
%
%
% Usage:
% [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
%   solveAnalyticalModelKM(Ta, k, L, qo, d, z, Lt, h_m, h_m_skin, omega, phi, ua, NTissueLayers, ...
%   A, b, A_FD, E, b_FD, A_q_skin, E_q_skin, b_q_skin, ...
%   matrix_temp, temp_offset, matrix_heat, heat_offset, ...
%   beta_initial, T_initial, tol, maxiter)
%
%
% Input:
% Ta: air temperature (oC)
% k: vector column with thermal conductivity values for layers (W/(m oC))
% L: vector column with length of layers (m)
% qo: vector of offsets heat fluxes (W/m2)
% d: diameter of the animal (m)
% z: elevation of the location (m).
% Lt: latitute of the location (decimal degrees).
% h_m: multiplicative factor for h to consider for possible measurement errors
% h_m_skin: multiplicative factor for h to consider for possible measurement errors
%           for convection heat flux at the skin.
% omega: proportion of convection heat transfer at the hair-coat.
% phi: additional proportion of convection heat transfer at the skin.
% ua: air velocity (m/s).
% NTissueLayers: number of tissue layers.
% A: matrix containing the relationship between betas
% b: vector containing offset for the relationship between betas
% A_FD: matrix that multiplies T for the KM model
% E: matrix that multiplies T^4 for the KM model
% b_FD: matrix that contains the vectors that will add to the KM model equation
% A_q_skin: vector that multiplies T to calculate heat flux at the surface of the skin from the KM model
% E_q_skin: vector that multiplies T^4 to calculate heat flux at the surface of the skin from the KM model
% b_q_skin: vector that adds to calculate heat flux at the surface of the skin from the KM model
% matrix_temp: matrix containing the relationship between betas and the temperatures
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
% beta_initial: First guess of betas
% T_initial: First guess of temperatures
% tol: minimum tolerance value to assume convergence of the optimization functions.
% maxiter: number of maximum iterations to run the solver
%
%
% Output:
% T_int: temperature points at the interfaces (oC)
% q_int: heat flux points at the interfaces (W/m2)
% beta: vector of beta
% obj: value of the objective function at convergence
% info: code of the solver at convergence
% iter: number of iterations to converge
% nf: output from the solver
% lambda: output from the solver
% eqConst: residues from the equality constraints
% 
% Assumptions:
% 1) one-dimentional
% 2) thermal parameters are static.
% 3) Pennes' bioheat equation to represent heat transfer through blood flow
%

%
% For Kowalski and Mitchell's model, this problem is divided in two. The first
% is the problem of the temperature inside the tissue. This has as boundary conditions
% rectal temperature at x = 0 and heat flux q_skin at x = L (where L is at the 
% skin surface). Then, give T_skin, we compute the temperature inside the hair coat
% using finite differences.
        
% Obtaining initial conditions
beta = beta_initial;
T = 273.15 + T_initial*ones(length(b_FD), 1);
Tsnew = matrix_temp(end,:)*beta + temp_offset(end);

for iteration = 1:maxiter
  % setting up the new skin-surface temperature
  Ts = Tsnew;
  T(1) = Ts + 273.15;
  
  % calculating the convection heat flux at the hair-coat surface
  Th = T(end) - 273.15;
  h = getConvHeatTransCoeff(Ta, Th, ua, d, z, Lt);
  h *= h_m*omega;
  
  A_FD(end, end) = -A_FD(end, end - 1) - h; % -k(end)/deltax - h; % A_FD(end, end - 1) = k(end)/deltax;
  b_FD(1) = -T(1);
  b_FD(end) = h*(Ta + 273.15);

  % setting up objective functions, etc.
  objFunction = @(T)objectiveFunctionKM(T, A_FD, E, b_FD);
  gradObjFunction = @(T)gradientObjectiveFunctionKM(T, A_FD, E, b_FD);
  HessObjFunction = @(T)HessianObjectiveFunctionKM(T, A_FD, E, b_FD);

  eqConstFunction = @(T)equalityConstraintFunctionKM(T, A_FD, E, b_FD);
  gradEqConstFunction = @(T)gradientEqualityConstraintFunctionKM(T, A_FD, E, b_FD);

  % Minimizing the proble to solve for T
  [T, obj, info, iter, nf, lambda] = sqp(T, ...
  {objFunction, gradObjFunction, HessObjFunction}, {eqConstFunction, gradEqConstFunction}, [], ... 
  -realmax, realmax, maxiter, tol);

  % Now we have a new estimation of skin surface heat flux due to radiation.
  % Set this estimation as a boundary condition and calculate the new skin-surface
  % temperature
  
  % calculating the convection heat flux at the skin-surface
  h = getConvHeatTransCoeff(Ta, Ts, ua, d, z, Lt);
  h *= h_m_skin*(1 - omega + phi);
  qconv = h*(Ts - Ta);
  
  % heat flux boundary condition
  b(end) = A_q_skin'*T + E_q_skin'*T.^4 + b_q_skin + qo(NTissueLayers + 1) + qconv;
  % Finding new betas
  beta = A\b;
  
  % Checking if we converged
  Tsnew = matrix_temp(end,:)*beta + temp_offset(end);
  if ( (Tsnew - Ts)^2 < tol)
    break;
  end
end

% calculating outputs
eqConst = eqConstFunction(T);
T_int = [matrix_temp*beta + temp_offset; T(2:end) - 273.15];
q_int = [matrix_heat*beta + heat_offset; 
          k((NTissueLayers + 1):end)'.*( T(1:(end-1)) - T(2:end) )/L(end)];


end
