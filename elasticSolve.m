function fem = elasticSolve(fem,tspan)
%ELASTICSOLVE  Solve the elastic equations by FEMLAB.
%
% SYNOPSIS :
%    fem = elasticSolve(fem,tspan)
%    Note: This function takes the 'fem' structure assembled by
%       ELMODELASSEMBLE and calles the FEMLAB PDE solver that solves the 
%       equation based on information stored in 'fem'. The solution is also 
%       stored in 'fem'.
%
% INPUT :
%    tspan : A 1D array of time steps (strictly increasing or decreasing) over 
%       which the time dependent elastic equation is integrated and the
%       solution is output. For static problem, pass [].

%Check if 'tspan' is a strictly increasing (or decreasing) array of numerical
% values:
if ~isempty(tspan)
   if ~isnumeric(tspan) | ndims(tspan) > 1
      error(['''tspan'' should be a 1D array of numerical values.  ' ...
         'See ELASTICSOLVE.']);
   elseif length(tspan) > 1
      if tspan(2)-tspan(1) > 0 % 'tspan' is increasing.
         sgn = 1; % see code 5 lines below for the purpose of 'sgn'.
      elseif tspan(2)-tspan(1) == 0
         error(['''tspan'' should be a strictly increasing or decreasing ' ...
            '1D array.  See ELASTICSOLVE.']);
      else % 'tspan' is decreasing.
         sgn = -1;
      end

      for k = 2:length(tspan)
         if sgn*(tspan(k)-tspan(k-1)) <= 0
            error(['''tspan'' should be a strictly increasing or ' ...
               'decreasing array.  See ELASTICSOLVE.']);
         end
      end
   end
end

% Evaluate initial condition
%fem.sol = femlin(fem,'Context','main');
fem.sol = femlin(fem);

