%[pdf, cdf] = stepModelODE(k, t) calculates the lifetime distribution for a 
% multi-step model of the form:
%
%    k1     k2      kn-1
% S1 --> S2 --> ... --> Sn
%
% INPUTS
%     k : vector of rate constants
%     t : time vector
%
% OUTPUTS
%   pdf : pdf of lifetimes
%   cdf : cdf of lifetimes
%
% Note: ODE solver-based implementation

% Francois Aguet, 2013

function [pdf, varargout] = stepModelODE(k, t)
S0 = [1 zeros(1,numel(k))];
sol = ode45(@(t,y) ksteps(t, y, k(:)), [0 t(end)], S0);
Y = deval(sol, t);
pdf = Y(end-1,:);
pdf = pdf/sum(pdf);%/(x(2)-x(1));
if nargout>1
    varargout{1} = Y(end,:);
end

function dy = ksteps(~, y, k)
S = -diag([k; 0]) + diag(k,-1);
dy = S*y;