function [bowtie]=testBowtie(boxesX, boxesY)
% TESTBOWTIE tests 1 or more quadrilaterals for self-intersection at a point
% (where opposite sides cross). It does not test if 2 sides are parallel.

% test quadrilaterals 1) convex, 2) sides 12 and 34 cross, 3) sides 23 and
% 41 cross, 4) sides 34 and 41 parallel (this case not bowtie, area is
% correct) boxesX=[1 1 1 2; 2 3 1 3; 3 1 2 3; 3 3 2 1; 1 1 1 2] boxesY=[1 1
% 1 1; 3 2 3 3; 2 3 1 1; 1 1 2 1; 1 1 1 1]

% this function is used in grid2boxCorners

sideCross=zeros(2,size(boxesX,2));

for i=1:2 % check 2 pairs of opposite edges

    if i==1 % check to see if sides 12 and 34 cross
        x=boxesX(1:4,:);
        y=boxesY(1:4,:);
    else % check to see if sides 23 and 41 cross
        x=boxesX(2:5,:);
        y=boxesY(2:5,:);
    end
    
    denom=(y(4,:)-y(3,:)).*(x(2,:)-x(1,:)) - (x(4,:)-x(3,:)).*(y(2,:)-y(1,:));
    % pairs of edges where denom is zero are parallel. to avoid divide by
    % zero warnings, make it small
    denom(denom==0)=10^-10;
    ua=((x(4,:)-x(3,:)).*(y(1,:)-y(3,:)) - (y(4,:)-y(3,:)).*(x(1,:)-x(3,:)))./denom;
    ub=((x(2,:)-x(1,:)).*(y(1,:)-y(3,:)) - (y(2,:)-y(1,:)).*(x(1,:)-x(3,:)))./denom;

    sideCross(i,ua>0 & ua<1 & ub>0 &ub<1)=1;
end
% bowtie gets 1 if box self-intersects at a point
bowtie=logical(sum(sideCross));
