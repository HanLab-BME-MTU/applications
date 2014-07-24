function [list]=PointP_upSample(StartList, multiplier, xy_Sig)
%PointP_upSample takes a point process and resamples it into a larger
%number of points. The points are chosen at random from startlist (multiple
%times) and each point is given a small randomly chosen displacement.
%
%Up sampling, STORM style each true inital point can be revisualized
%several times potentially with slightly different positions.
%The following function randomly resamples the initial list and gives each
%point a small delta in localized position.
%
% Jeffrey L. Werbin 
% Harvard Medical School
%
% Last update: 6/6/2011

%sets the size of the output list
n = size(StartList,1);
nF = floor(multiplier*n);

%Choosing random points in indexed form
ind = random('unid',n,[nF,1]);

list = StartList(ind,:)+random('norm',0,xy_Sig,[nF,2]);

end


