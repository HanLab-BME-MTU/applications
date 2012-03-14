% Spatial sampling functions
%
% Jeffrey L. Werbin 
% Harvard Medical School
%
% Last update: 6/6/2011
%
% This file is one of several functions for artifically resampling point
% process data. 
%
% All the functions take the point process as a list of points of 2D
% positions and return a list in the same format.
%

%Up sampling, STORM style each true inital point can be revisualized
%several times potentially with slightly different positions.
%The following function randomly resamples the initial list and gives each
%point a small delta in localized position.
function [list]=PointP_upSample(StartList, multiplier, xy_Sig)
%PointP_upSample takes a point process and resamples it into a larger
%number of points. The points are chosen at random from startlist (multiple
%times) and each point is given a small randomly chosen displacement.

%Sets the size of the output list
n = size(StartList,1);
nF = floor(multiplier*n);

%Choosing random points in indexed form
ind = random('unid',n,[nF,1]);

list = StartList(ind,:)+random('norm',0,xy_Sig,[nF,2]);

end


