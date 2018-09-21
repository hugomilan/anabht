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
% File:   ensembleTraining.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Calls the functions that train the ensembles

% select the sets to train
% 0: deterministic set. Other numbers: random set
sets = 0:10001;
% solve for a specific model or for all
solveFor = "KM"; %"all"

%%%%%%%%%%%%%% NO NEED TO MODIFY BELOW THIS LINE %%%%%%%%%%%%%%

if exist("src") != 0
  addpath("src");
end

ensembleSetsMSE(sets, solveFor);