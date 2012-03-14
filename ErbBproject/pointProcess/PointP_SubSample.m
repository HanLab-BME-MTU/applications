% Spatial sampling functions
%
% Jeffrey L. Werbin 
% Harvard Medical School
%
% Last update: 7/7/2011
%
% This file includes several functions for artifically resampling point
% process data. 
%
% All the functions take the point process as a list of points of 2D
% positions and return a list in the same format.
%


function [list] = PointP_SubSample(StartList, percent)
%PointP_SubSample takes a list of coordinates defining a point process and
%picks and returns only perecent*numpts(startlist) of those points where
%0<percent<1


 if (percent>1.0 || percent < 0.0)
     error('percent must be between 0.0 and 1.0');
     return;
 end

n = size(StartList,1);

%disp(n)

%rand_list is a list of indexs to StartList
[S, rand_list] = sort(random('unif',1,10,[n,1]));

%disp(rand_list)

%taking percent*n at random from StartList using rand_list 
list = StartList(rand_list(1:floor(percent*n)),:);

end


