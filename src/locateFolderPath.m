function folderPath = locateFolderPath(folderName)
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
% File:   locateFolderPath.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 26, 2018.
%
%
% Function description:
% Given the name of the folder, it locates where the folder is and returns the
% location of the folder.
%
%
%
% Usage:
% folderPath = locateFolderPath(folderName)
%
% Input:
% folderName: Name of the folder
%
% Output:
% folderPath: string containing the location of the folder

% adding portability between MATLAB and Octave
if (exist('OCTAVE_VERSION', 'builtin') ~= 0)
  % we are executing in Octave
  testFunction = @(name)!isempty( stat(name) );
else
  % we are executing in MATLAB
  testFunction = @(name) (exist(name) == 7);
end

folderPath = folderName;
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

% may be we need to add "src"
folderPath = ["src/" folderName];
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

% maybe we need to go some directories down
folderPath = ["../" folderName];
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

folderPath = ["../" folderPath];
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

folderPath = ["../" folderPath];
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

folderPath = ["../" folderPath];
% if the name corresponds to a folder, we have found it
if testFunction(folderPath)
   return 
end

error(["Did not find the location of the folder \"" folderName "\""]);

end
