function val = expintn(n, x)
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
% File:   expintn.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 26, 2018.
%
%
% Function description:
% Calculates the exponential integral of order n given x
%
%
% Usage:
% val = expintn(n, x)
%
% Input:
% n: order of the exponential integral
% x: value of x to calculate the integral
%
% Output:
% val: the value of the exponential integral of order n given x
% 

  if (n == 0)
    val = exp(-x)./x;
  else
    val = expint(x);
    for ii = 1:(n - 1)
      val = 1/ii*( exp(-x) - x.*val);
    end
  end
   
end