function b = updateBetaParametersProblemDependent(b, wb, Tb, Tr)
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
% File:   getBetaParametersProblemIndependent.m
% Author: Hugo Fernando Maia Milan
% Email:  hugofernando@gmail.com
%
% Created on April 11, 2018.
%
%
% Function description:
% Calculates the matrix and vector offset for beta given thermal parameters that
% are problem independent
%
%
% Usage:
% b = updateBetaParametersProblemDependent(b, wb, Tb, Tr)
%
% Input:
% b: vector containing offset for the relationship between betas
% wb: vector column with blood perfusion values for layers (m3/(s m3))
% Tb: vector column with blood temperature values for layers (oC)
% Tr: Rectal temperature (oC)
%
%
% Output:
% b: vector containing offset for the relationship between betas
% 

% Calculations:

% first boundary condition
% applying first boundary condition
%A(1, 2) = 1; % beta_2i
b(1) += Tr - Tb(1); % + To(1) - qtotal(1)/( wb(1)*rhob(1)*cb(1) )

% equations for ensuring the same temperature and heat flux at the interface 
% between the mediums
for a2 = 1:(size(wb,2) - 1)
  if (wb(a2) == 0 && wb(a2 + 1) == 0)
    % forcing continuity of temperature
    % A(2*a2, 2*a2 - 1) = L(a2); % beta_1i
    % A(2*a2, 2*a2 + 0) = 1; % beta_2i
    % A(2*a2, 2*a2 + 1) = 0; % beta_1(i+1)
    % A(2*a2, 2*a2 + 2) = -1; % beta_2(i+1)
    % b(2*a2) = qtotal(a2)/( 2*k(a2) )*L(a2)^2 + To(a2 + 1);
    
    % forcing continuity of heat flux
    % A(2*a2 + 1, 2*a2 - 1) = -k(a2); % beta_1i
    % % A(2*a2 + 1, 2*a2 + 0) = 0; % beta_2i
    % A(2*a2 + 1, 2*a2 + 1) = k(a2+1); % beta_1(i+1)
    % A(2*a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    % b(2*a2 + 1) = -qtotal(a2)*L(a2) + qo(a2 + 1);
    
  elseif (wb(a2) == 0 && wb(a2 + 1) != 0)
    % forcing continuity of temperature
    % A(2*a2, 2*a2 - 1) = L(a2); % beta_1i
    % A(2*a2, 2*a2 + 0) = 1; % beta_2i
    % A(2*a2, 2*a2 + 1) = 0; % beta_1(i+1)
    % A(2*a2, 2*a2 + 2) = -1; % beta_2(i+1)
    b(2*a2) += Tb(a2 + 1); % + qtotal(a2)/( 2*k(a2) )*L(a2)^2 + qtotal(a2+1)/( wb(a2+1)*rhob(a2+1)*cb(a2+1) ) + To(a2 + 1)
    
    % forcing continuity of heat flux
    % A(2*a2 + 1, 2*a2 - 1) = -k(a2); % beta_1i
    % A(2*a2 + 1, 2*a2 + 0) = 0; % beta_2i
    % A(2*a2 + 1, 2*a2 + 1) = sqrt( k(a2+1)*wb(a2+1)*rhob(a2+1)*cb(a2+1) ); % beta_1(i+1)
    % A(2*a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    % b(2*a2 + 1) = -qtotal(a2)*L(a2) + qo(a2 + 1);
    
  elseif (wb(a2) != 0 && wb(a2 + 1) == 0)
    % forcing continuity of temperature
    % A(2*a2, 2*a2 - 1) = sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
    % A(2*a2, 2*a2 + 0) = cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
    % A(2*a2, 2*a2 + 1) = 0; % beta_1(i+1)
    % A(2*a2, 2*a2 + 2) = -1; % beta_2(i+1)
    b(2*a2) += - Tb(a2); % - qtotal(a2)/( wb(a2)*rhob(a2)*cb(a2) ) + To(a2 + 1);
    
    % forcing continuity of heat flux
    % A(2*a2 + 1, 2*a2 - 1) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
    % A(2*a2 + 1, 2*a2 + 0) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
    % A(2*a2 + 1, 2*a2 + 1) = k(a2+1); % beta_1(i+1)
    % A(2*a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    % b(2*a2 + 1) = qo(a2 + 1);
    
  else % (wb(a2) != 0 & wb(a2 + 1) != 0)
    % forcing continuity of temperature
    % A(2*a2, 2*a2 - 1) = sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
    % A(2*a2, 2*a2 + 0) = cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
    % A(2*a2, 2*a2 + 1) = 0; % beta_1(i+1)
    % A(2*a2, 2*a2 + 2) = -1;
    b(2*a2) += - Tb(a2) + Tb(a2+1); % - qtotal(a2)/( wb(a2)*rhob(a2)*cb(a2) ) + qtotal(a2+1)/( wb(a2+1)*rhob(a2+1)*cb(a2+1) ) + To(a2 + 1);
    
    % forcing continuity of heat flux
    % A(2*a2 + 1, 2*a2 - 1) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*cosh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_1i
    % A(2*a2 + 1, 2*a2 + 0) = -sqrt( k(a2)*wb(a2)*rhob(a2)*cb(a2) )*sinh( L(a2)*sqrt( wb(a2)*rhob(a2)*cb(a2)/k(a2) ) ); % beta_2i
    % A(2*a2 + 1, 2*a2 + 1) = sqrt( k(a2+1)*wb(a2+1)*rhob(a2+1)*cb(a2+1) ); % beta_1(i+1)
    % A(2*a2 + 1, 2*a2 + 2) = 0; % beta_2(i+1)
    % b(2*a2 + 1) = qo(a2 + 1);
  end    
end

end
