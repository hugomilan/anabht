function [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
   solveAnalyticalModelSkin(Ta, d, z, Lt, h_m, ua, A, b, TMR, epsilon, qo, ...
   matrix_temp, temp_offset, matrix_heat, heat_offset, ...
   beta_initial, tol = sqrt(eps), maxiter = 1e3)
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
% File:   solveAnalyticalModelSkin.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 28, 2018.
%
%
% Function description:
% Given the inputs that describe the environment and the problem, this function
% solves the analytical problem.
%
%
% Usage:
% [T_int, q_int, beta, obj, info, iter, nf, lambda, eqConst] = ...
%    solveAnalyticalModelSkin(Ta, d, z, Lt, h_m, ua, A, b, TMR, epsilon, qo, ...
%    matrix_temp, temp_offset, matrix_heat, heat_offset, ...
%    beta_initial, tol = sqrt(eps), maxiter = 1e3)
%
%
% Input:
% Ta: air temperature (oC)
% d: diameter of the animal (m)
% z: elevation of the location (m).
% Lt: latitute of the location (decimal degrees).
% h_m: multiplicative factor for h to consider for possible measurement errors
% ua: air velocity (m/s).
% A: matrix containing the relationship between betas
% b: vector containing offset for the relationship between betas
% TMR: Mean radiant temperature (oC)
% epsilon: emissivity of the animal surface
% qo: vector of offsets heat fluxes (W/m2)
% matrix_temp: matrix containing the relationship between betas and the temperatures
% temp_offset: vector containing offset for the relationship between betas and the temperatures
% matrix_heat: matrix containing the relationship between betas and heat fluxes
% heat_offset: vector containing offset for the relationship between betas and heat fluxes
% beta_initial: First guess of betas
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
 

objFunction = @(beta)objectiveFunctionSkin(beta, matrix_temp, temp_offset, ...
      matrix_heat, heat_offset, Ta, TMR, ua, d, epsilon, qo(end), z, Lt, h_m);
gradObjFunction = @(beta)gradientObjectiveFunctionSkin(beta, matrix_temp, temp_offset, ...
      matrix_heat, heat_offset, Ta, TMR, ua, d, epsilon, qo(end), z, Lt, h_m);
HessObjFunction = @(beta)HessianObjectiveFunctionSkin(beta, matrix_temp, temp_offset, ...
      matrix_heat, heat_offset, Ta, TMR, ua, d, epsilon, qo(end), z, Lt, h_m);
      
eqConstFunction = @(beta)equalityConstraintFunctionSkin(beta, A, b);
gradEqConstFunction = @(beta)gradientEqualityConstraintFunctionSkin(beta, A, b);

[beta, obj, info, iter, nf, lambda] = sqp(beta_initial, ...
{objFunction, gradObjFunction, HessObjFunction}, {eqConstFunction, gradEqConstFunction}, [], ...
-realmax, realmax, maxiter, tol);

eqConst = eqConstFunction(beta);

T_int = matrix_temp*beta + temp_offset;
q_int = matrix_heat*beta + heat_offset;


end
