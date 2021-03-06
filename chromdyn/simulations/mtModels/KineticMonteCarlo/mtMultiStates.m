function [traj, stateIndex, errFlag] = mtMultiStates(model, initLen, time)
% MTMULTISTATES uses kinetic Monte Carlo simulation to generate an MT
% trajectory, assuming the growth, shrinkage and pause time are exponentially
% distributed.
%
% SYNOPSIS [traj, stateIndex, errFlag] = mtMultiStates(model, initLen, time)
%
% INPUT 
%   Mandatory
%       model           : Array of states, each of which containing the
%                         following fields
%           .name       : String representing the state name
%           .speed      : Mean speed [microns/minute]
%                         positive: growth, negative: shrinkage
%           .nTransitions : Number of possible transitions from the
%                         current state
%           .transit    : Array of transitions, each of which containing
%                         the following fields
%               .rate   : Transition frequency [1/minute]
%               .dS     : Integer j denoting a destination state, model(j)
%       initLen         : Initial MT length [microns]
%       time            : Total simulation time [seconds]
%
% OUTPUT
%       traj            : Array comprising two columns 
%                         1st column: time at transitions [seconds]
%                         2nd column: MT length [microns]
%       stateIndex      : Integer array chronicling the transitions
%
% Pei-hsin Hsu, December 2008

errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input. Initialize output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 3
    disp('--mtKinetMonteCarlo: Wrong number of arguments');
    errFlag = 1;
    return
end

traj(1, :) = [0, initLen];
stateIndex(1) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kinetic Monte Carlo simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;

% End the simulation when cumulative time exceeds the argument time
while (traj(i, 1) < time)
    
    % Perform a one-step state transition
    [dS, dT, dL] = stateTransition(model(stateIndex(i)));
    
    % Update output
    stateIndex(i + 1) = dS;
    traj(i + 1, :) = [traj(i, 1) + dT, traj(i, 2) + dL];
    
    i = i + 1;
end
