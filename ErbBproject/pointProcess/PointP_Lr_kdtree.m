function [Lr]=PointP_Lr_kdtree(pos,r)
%PointP_Lr_Cross takes two related point processes and evaulates
% the average number of points in posB that are with a distance of r(i) from a
% point in posA.
% Implemented using KDtreeBallQuery.
%
%Inputs: 
%       posA, reference postion list (nx2)
%       posB, queried list (mx2)
%          r, a vector of radii to consider
%
%Outputs:
%       Lr, the Lr statistic for 
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/30/2013
%

%edge = 1-(4/(3*pi))*((r/L)+(r/W))+((11/(3*pi))-1)*(r^2/(L*W));
edge=1;

%The number of total points
%n = size(posA,1);
%m = size(posB,1);

%pos = [posA ; posB]; %concatinating the two lists allows for processing with one dist matrix

%find boundaries
limX = [max(pos(:,1)), min(pos(:,1))];
limY = [max(pos(:,2)), min(pos(:,2))];
W = limX(1)-limX(2);
L = limY(1)-limY(2);

%to analyze only points > max(rA,rB) from the boundries of the rectangle
%removes edge effect problems

Rmax = max(r);
x = pos(:,1)-limX(2);
y = pos(:,2)-limY(2);

%Here we find all the points in pos A that are at least Rmax away from the
%bounderies. To help prevent edge effects.

%ind = find( (x > Rmax) & (x < W -Rmax) & (y > Rmax) & (y < L-Rmax));
ind = (x > Rmax) & (x < W -Rmax) & (y > Rmax) & (y < L-Rmax);
clear x y;

%keeps only points that meet this criteria as query points
posA= pos(ind,:);

%Uncomputed values are returned as NaNs
Lr = NaN(numel(r),1);

%Average density within the total area considered
lambda = numel(pos(:,1))/(W*L);

%c
tic
for i=1:numel(r)

    %for each r it finds the inPnts that are <r(i) from each query point
    [idx, dist] = KDTreeBallQuery(pos,posA,r(i));
    
    %idx is a cell array with indicies of all posB within r(i) of posA
    % we only need to know the number of points in idx for each point
    temp = cellfun(@numel,idx);
    temp = sqrt(sum(temp)/(lambda*pi*(numel(temp)-1)));
    
    
    %out(i,1)= sum(tempB)/(edge*pi*rB^2);
    %out(i,2)= sum(tempA)/(edge*pi*rA^2);
    Lr(i)= temp;
    
end
toc


end
    