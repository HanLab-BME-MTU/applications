function [dS, dT, dL] = stateTransition(state)
% STATETRANSITION performs a one-step transiton in kinetic Monte Carlo
% simulation.
% 
% SYNOPSIS [dS, dT, dL] = stateTransition(state)
%
% INPUT
%   Mandatory
%       state           : Structure comprising the following fields
%           .name       : String representing the state name
%           .speed      : Mean speed [microns/minute]
%                         positive: growth, negative: shrinkage
%           .nTransitions : Number of possible transitions from the
%                         state
%           .transit    : Array of transitions, each of which containing
%                         the following fields
%               .rate   : Transition frequency [1/minute]
%               .dS     : Integer j denoting a destination state, model(j)
%
% OUTPUT
%       dS              : Integer denoting the destination state after
%                         one-step transition
%       dT              : Time interval [seconds] of the source state, 
%                         assuming exponential distribution. 
%       dL              : Length increment [microns] during dT
%
% Pei-hsin Hsu, December 2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the cumulative rates c(i)       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = state.nTransits;

c(1) = state.transit(1).rate;
if n > 1
    for i = 2:n, c(i) = c(i-1) + state.transit(i).rate; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the selection of destination state dS   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = find(rand - c/c(n) < 0, 1, 'first');
dS = state.transit(j).dS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the time interval of the source state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implementation notes:
%
% (1) Theorem:
%     If X is a continuous random variable, and Y = g(X) is strictly
%     increasing on the range of X, then their respective cumulative
%     distribution function (c.d.f), Fx and Fy, satisfies
%     Fy(y) = Fx((inverse g)(y))
% (2) Corollary:
%     If 0 < g(x) < 1 is strictly increasing for all x, then
%     X = (inverse g)(rand) has c.d.f g(x)
% (3) g(t) = 1 - exp(-r * t) is the c.d.f of an exponential random variable
%     with parameter r. Let y = g(t), then 
%     t = (inverse g) (y) = - (1/r) * log(1 - y)

t = - (1/c(n)) * log(rand);

% Calculate the microtubule length increment
dL = state.speed * t;

% Change the unit of time interval to [second]
dT = t * 60;