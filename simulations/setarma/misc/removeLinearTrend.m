function trajectoryOut = removeLinearTrend(trajectory)
%REMOVELINEARTREND removes a linear trend from a data series 
%   by doing a robust linear fit and subsequent subtraction. Furthermore,
%   it sets the mean to zero
%
% SYNOPSIS  outTrajectory = removeLinearTrend(trajectory)
%
% INPUT     trajectory: any 1D data series vector. Can contain NaNs
%
% OUTPUT    trajectoryOut: the above data series without linear trend and
%                          with zero mean
%
%
%c: jonas 9/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

% get size of original trajectory
trajSizeOri = size(trajectory);

% make into column vector
trajectory = trajectory(:);

% and get the new size
trajSize = size(trajectory);

%============


%=================
% LINEAR FIT
%=================

% we fit A*x=B+E

% prepare matrices
A = [ones(trajSize),[1:trajSize(1)]'];
B = trajectory;

% robust linear fit
xLin = linearLeastMedianSquares(A,B);

%=================

%===================
% SUBTRACT & RETURN
%===================

% subtract trend
trajectory = trajectory - A*xLin;

% subtract mean
trajectoryOut = trajectory - nanmean(trajectory);