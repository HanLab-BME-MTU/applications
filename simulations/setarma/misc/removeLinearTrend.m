function [trajectoryOut, linFit] = removeLinearTrend(trajectory, sigmaTrajectory)
%REMOVELINEARTREND removes a linear trend from a data series 
%   by doing a robust linear fit and subsequent subtraction. Furthermore,
%   it sets the mean to zero
%
% SYNOPSIS  outTrajectory = removeLinearTrend(trajectory, sigmaTrajectory)
%
% INPUT     trajectory     : any 1D data series vector. Can contain NaNs
%           sigmaTrajectory: (optional) The corresponding uncertainties
%
% OUTPUT    trajectoryOut: the above data series without linear trend and
%                          with zero mean
%           linFit       : intercept and slope of linear fit
%
%
% currently, there is no error propagation implemented.
%
%c: jonas 9/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

% get size of original trajectory
trajSizeOri = size(trajectory);

% make into column vector
trajectory = returnRightVector(trajectory);

% and get the new size
trajSize = size(trajectory);

% check for sigma
if nargin > 1 && ~isempty(sigmaTrajectory)
    sigmaTrajectory = returnRightVector(sigmaTrajectory);
else
    sigmaTrajectory = ones(trajSize);
end
    
%============


%=================
% LINEAR FIT
%=================

% we fit A*x=B+E

% prepare matrices
A = [ones(trajSize),[1:trajSize(1)]'];
B = trajectory;

% robust linear fit
linFit = linearLeastMedianSquares(A,B,diag(sigmaTrajectory));

%=================

%===================
% SUBTRACT & RETURN
%===================

% subtract trend
trajectory = trajectory - A*linFit;

% subtract mean
trajectoryOut = trajectory - nanmean(trajectory);
