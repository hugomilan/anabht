function MSE_selected = naiveSetsSelection(MSE, numberOfSets)
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
% File:   naiveSetsSelection.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 19, 2018.
%
%
% Function description:
% Selects best sets. The naive set selection will simple return MSE(1:numberOfSets, :)
%
%
%
% Usage:
% MSE_selected = naiveSetsSelection(MSE, numberOfSets)
%
% Input:
% MSE: matrix containing the MSE, where the first column is the MSE and the 
%      second column is the number of the set.
% numberOfSets = number of sets to return.
%
% Output:
% MSE_selected: sets selected

MSE_selected = MSE(1:numberOfSets, :);

end
