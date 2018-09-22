function [T_plot q_plot x_plot] = getPlotValues(beta, nPointsWithinLayer, k, L, ...
   qtotal, wb, rhob, cb, Tb, T_int, q_int, NTissueLayers = length(wb), FD_Calc = 0);
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
% File:   getPlotValues.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com

% Created on February 11, 2018.
%
%
% Function description:
% Calculates values to plot for temperature and heat flux
%
%
% Usage:
% [T_plot q_plot x_plot] = getPlotValues(beta, nPointsWithinLayer, k, L, ...
%   qtotal, wb, rhob, cb, Tb)
%
% Input:
% beta: value of betas for each layer
% nPointsWithinLayer: number of points to sample from a layer
% k: vector column with thermal conductivity values for layers (W/(m oC))
% L: vector column with legth of layers (m)
% qtotal: vector column with total heat source values for layers (W/m3)
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% rhob: vector column with blood density values for layers (kg/m2)
% cb: vector column with blood specific heat values for layers (J/(kg oC))
% Tb: vector column with blood temperature values for layers (oC)
% T_int: temperature points at the interfaces (oC)
% q_int: heat flux points at the interfaces (W/m2)
% NTissueLayers: number of tissue layers.
% FD_Calc: variable that tells if heat fluxes and temperatures after NTissueLayer should be
%          computed using linear approximations (FD_Calc = 1; when the problem
%          inside the hair coat is solved with finite differences) or if it should
%          be computed using the analytical equations.
%
%
% Output:
% T_plot: temperature points extracted from layers (oC)
% q_plot: heat flux points extracted from layers (W/m2)
% x_plot: position in space of the sampled temperature points (m). First layer is at x = 0
% 

% Calculations:

T_plot = zeros(size(k,2)*nPointsWithinLayer + 1,1);
q_plot = zeros(size(k,2)*nPointsWithinLayer + 1,1);
x_plot = zeros(size(k,2)*nPointsWithinLayer + 1,1);

x_medium = zeros(nPointsWithinLayer,1);
for a1 = 1:size(k,2)
  x_medium(1) = 0;
  if (a1 != 1)
    x_plot(1 + (a1-1)*nPointsWithinLayer) = L(a1-1)/nPointsWithinLayer + ...
                                             x_plot((a1-1)*nPointsWithinLayer);
  end
  
  for a2 = 2:nPointsWithinLayer
    x_medium(a2) = x_medium(a2-1) + L(a1)/nPointsWithinLayer;
    x_plot(a2 + (a1-1)*nPointsWithinLayer) = L(a1)/nPointsWithinLayer + ...
                                             x_plot(a2 - 1 + (a1-1)*nPointsWithinLayer);
  end

  if (a1 > NTissueLayers && FD_Calc)
    T_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
          T_int(a1) + (T_int(a1 + 1) - T_int(a1))/x_medium(end)*x_medium;
          
    q_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
          q_int(a1) + (q_int(a1 + 1) - q_int(a1))/x_medium(end)*x_medium;
  elseif (wb(a1) == 0)
    T_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
        -qtotal(a1)/(2*k(a1))*(x_medium).^2 + beta(2*a1 - 1).*x_medium + beta(2*a1);
        
    q_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
        qtotal(a1)*x_medium + -k(a1)*beta(2*a1 - 1);
        
  else % wb(a1) != 0
    T_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
        Tb(a1) + qtotal(a1)/(wb(a1)*rhob(a1)*cb(a1)) + ...
        beta(2*a1 - 1)*sinh(x_medium*sqrt( wb(a1)*rhob(a1)*cb(a1)/k(a1) ) ) ...
        + beta(2*a1)*cosh(x_medium*sqrt( wb(a1)*rhob(a1)*cb(a1)/k(a1) ) );
        
    q_plot( ((a1-1)*nPointsWithinLayer + 1):((a1)*nPointsWithinLayer) ) = ...
        - beta(2*a1 - 1)*sqrt( k(a1)*wb(a1)*rhob(a1)*cb(a1) )*cosh(x_medium*sqrt( wb(a1)*rhob(a1)*cb(a1)/k(a1) ) ) ...
        - beta(2*a1)*sqrt( k(a1)*wb(a1)*rhob(a1)*cb(a1) )*sinh(x_medium*sqrt( wb(a1)*rhob(a1)*cb(a1)/k(a1) ) );
  end
end

% solving for the last point
% a1 = size(k,2);
x_medium_last = x_medium(end) + L(a1)/nPointsWithinLayer;
x_plot(end) = L(a1)/nPointsWithinLayer + x_plot(end - 1);

T_plot(end) = T_int(end);
q_plot(end) = q_int(end);

end
